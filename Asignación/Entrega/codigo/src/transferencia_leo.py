"""Transferencia tangencial LEO–Sol del material del curso."""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .constantes import DAY, GM, MU, RAD
from .hohmann import R1, R2

EARTH_RADIUS = 6378.137
EARTH_MU = 398600.4418
EARTH_SOI = R1*(EARTH_MU/MU)**.4
TRANSFER_AXIS = (R1+R2)/2
NOMINAL_TIME = np.pi*np.sqrt(TRANSFER_AXIS**3/MU)
EARTH_SOLAR_SPEED = np.sqrt(MU/R1)
NOMINAL_EXCESS = np.sqrt(MU*(2/R1-1/TRANSFER_AXIS))-EARTH_SOLAR_SPEED
JUPITER_RATE = np.sqrt(MU/R2**3)
JUPITER_PHASE = np.pi-JUPITER_RATE*NOMINAL_TIME


def hohmann_burn(altitude):
    circular = np.sqrt(EARTH_MU/(EARTH_RADIUS+altitude))
    return np.sqrt(NOMINAL_EXCESS**2+2*circular**2)-circular


def circle(radius, angles):
    return radius*np.column_stack((np.cos(angles), np.sin(angles)))


def tangent_arc(mu, radius, speed, end_angle):
    """Muestreo en anomalía verdadera desde un punto de retorno."""
    p = (radius*speed)**2/mu
    e = p/radius-1
    angles = np.linspace(0, end_angle, 2401)
    radii = p/(1+e*np.cos(angles))
    xy = radii[:,None]*circle(1,angles)
    return xy


@dataclass
class Transfer:
    circular_speed: float
    burn_speed: float
    earth_excess: float | None
    leo: np.ndarray
    escape: np.ndarray
    solar: np.ndarray | None
    arrival_days: float | None
    phase_error_deg: float | None
    nominal: bool


def transfer(altitude, delta_v):
    nominal_burn = hohmann_burn(altitude)
    nominal = abs(delta_v-nominal_burn) < 1e-8
    if nominal:
        # El control redondea el valor prescrito a ocho decimales.
        delta_v = nominal_burn
    radius = EARTH_RADIUS+altitude
    circular = np.sqrt(EARTH_MU/radius)
    burn = circular+delta_v
    energy_twice = burn**2-2*circular**2
    p = (radius*burn)**2/EARTH_MU
    e = p/radius-1
    if 0 < e < 1 and p/(1-e) >= EARTH_SOI:
        earth_angle = np.arccos(np.clip((p/EARTH_SOI-1)/e,-1,1))
    elif e >= 1:
        earth_angle = np.arccos(np.clip((p/EARTH_SOI-1)/e,-1,1))
    else:
        earth_angle = np.pi*2
    leo = circle(radius,np.linspace(0,2*np.pi,361))
    escape = tangent_arc(EARTH_MU,radius,burn,earth_angle)
    if energy_twice < 0:
        return Transfer(circular,burn,None,leo,escape,None,None,None,False)

    excess = np.sqrt(energy_twice)
    speed = EARTH_SOLAR_SPEED+excess
    p = (R1*speed)**2/MU
    e = p/R1-1
    if nominal:
        speed = EARTH_SOLAR_SPEED+NOMINAL_EXCESS
    reaches = nominal or p >= R2*(1-e)
    if reaches:
        angle = np.pi if nominal else np.arccos(np.clip((p/R2-1)/e,-1,1))
    else:
        angle = 2*np.pi
    solar = tangent_arc(MU,R1,speed,angle)
    if not reaches:
        return Transfer(circular,burn,excess,leo,escape,solar,None,None,False)

    # dt/dν = r²/L por la ley de áreas. La cuadratura sirve también en la
    # frontera parabólica y para excesos que producen una hipérbola solar.
    radii_sq = np.sum(solar**2,axis=1)
    elapsed = np.trapezoid(radii_sq,dx=angle/(len(solar)-1))/(R1*speed)
    phase = JUPITER_PHASE+JUPITER_RATE*elapsed
    error = np.degrees(np.arctan2(np.sin(angle-phase),np.cos(angle-phase)))
    return Transfer(circular,burn,excess,leo,escape,solar,elapsed/DAY,error,nominal)


@dataclass
class JovianScattering:
    periapsis_km: float
    eccentricity: float
    impact_parameter_km: float
    turn_rad: float
    incoming: np.ndarray
    outgoing: np.ndarray
    planet_velocity: np.ndarray
    heliocentric_exit: np.ndarray
    energy_change: float
    hyperbola: np.ndarray


@lru_cache(maxsize=1)
def nominal_scattering():
    """Dispersión prógrada hacia afuera en el afelio nominal de Hohmann."""
    jupiter_mu = GM['J']
    periapsis = 10*RAD['J']
    planet_velocity = np.array([0.,-np.sqrt(MU/R2)])
    approach = np.array([0.,-R1/R2*(EARTH_SOLAR_SPEED+NOMINAL_EXCESS)])
    incoming = approach-planet_velocity
    excess = np.linalg.norm(incoming)
    eccentricity = 1+periapsis*excess**2/jupiter_mu
    turn = 2*np.arcsin(1/eccentricity)
    asymptote = np.arccos(-1/eccentricity)
    periapsis_direction = np.arctan2(incoming[1],incoming[0])-(np.pi-asymptote)
    outgoing = excess*circle(1,[periapsis_direction+asymptote])[0]
    p = periapsis*(1+eccentricity)
    # Ventana local de ocho periapsis. Las velocidades del empalme son
    # asintóticas y se calculan aparte de esta ventana finita.
    edge = np.arccos((p/(8*periapsis)-1)/eccentricity)
    anomalies = np.linspace(-edge,edge,1201)
    hyperbola = (p/(1+eccentricity*np.cos(anomalies)))[:,None]*circle(1,anomalies+periapsis_direction)
    return JovianScattering(periapsis,eccentricity,p/np.sqrt(eccentricity**2-1),
                           turn,incoming,outgoing,planet_velocity,planet_velocity+outgoing,
                           float(planet_velocity@(outgoing-incoming)),hyperbola)


@dataclass
class SaturnLeg:
    approach: np.ndarray
    departure: np.ndarray
    relative_path: np.ndarray
    encounter_path: np.ndarray
    periapsis_position: np.ndarray
    saturn_at_periapsis: np.ndarray
    days_to_periapsis: float
    days_after_periapsis: float


@lru_cache(maxsize=1)
def nominal_saturn_leg():
    """Solución Sol–Saturno antes, durante y después del sobrevuelo."""
    from .optimizacion import design

    _, propagation, system, closest = design(return_refined=True)
    # La transferencia llega a Júpiter en -x y el diseño parte en +x.
    window = 5*DAY
    before, after = closest-window, closest+window
    local_times = np.linspace(before,after,601)
    saturn = system.planets['S']
    encounter_path = -propagation(local_times)[:2].T
    local = encounter_path+np.array([saturn(t)[:2] for t in local_times])
    if not np.allclose(-propagation(0)[3:5],nominal_scattering().heliocentric_exit,atol=1e-10,rtol=0):
        raise RuntimeError('La continuación solar no coincide con la salida del scattering joviano')
    arrival_times = np.linspace(0,before,601)
    end = propagation.ts[-1]
    departure_times = np.linspace(after,end,401)
    return SaturnLeg(-propagation(arrival_times)[:2].T,-propagation(departure_times)[:2].T,
                     local,encounter_path,-propagation(closest)[:2],-saturn(closest)[:2],
                     closest/DAY,(end-closest)/DAY)
