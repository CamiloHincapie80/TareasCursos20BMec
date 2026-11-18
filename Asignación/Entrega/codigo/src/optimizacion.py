"""Diseño circular coplanar con fase de Saturno libre, desde Hohmann.
Se conserva la hipótesis de asistencia instantánea asintótica del original.
El objetivo es un acercamiento finito de 500000 km en propagación solar;
no se identifica este diseño preliminar con un periapsis saturniano perturbado.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq,minimize_scalar
from .constantes import MU,GM,RAD,DAY
from .dinamica import grav,norm,System
from .hohmann import calculate,R2,R3
from .encuentros import minimum

def design(return_refined=False):
    h=calculate(); u=h['vinf_jupiter_kms']; V=h['v_jupiter_kms']
    rp=10*RAD['J']; delta=2*np.arcsin(1/(1+rp*u*u/GM['J']))
    incoming=np.array([0.,-u,0.]); outgoing=np.array([u*np.sin(delta),-u*np.cos(delta),0.])
    def rhs(t,y): return np.r_[y[3:],grav(MU,y[:3])]
    def at_saturn(t,y): return norm(y[:3])-R3
    at_saturn.terminal=True; at_saturn.direction=1
    y0=np.r_[[R2,0.,0.],outgoing+[0,V,0]]
    sol=solve_ivp(rhs,(0,15*365.25*DAY),y0,events=at_saturn,rtol=1e-12,atol=1e-8,dense_output=True,max_step=DAY*5)
    if not len(sol.t_events[0]): raise RuntimeError('El diseño no cruza la órbita de Saturno')
    tc=sol.t_events[0][0]
    # Extender después del cruce para medir un mínimo interior.
    sol=solve_ivp(rhs,(0,tc+300*DAY),y0,rtol=1e-12,atol=1e-8,dense_output=True,max_step=DAY*5)
    theta=np.arctan2(sol.sol(tc)[1],sol.sol(tc)[0]); n=np.sqrt(MU/R3**3)
    phase0=theta-n*tc
    def planet(t,phase): return R3*np.array([np.cos(phase+n*t),np.sin(phase+n*t),0.])
    def miss(phase): return minimum(sol.sol,lambda t:planet(t,phase),tc-250*DAY,tc+250*DAY)
    target=5e5
    phase=brentq(lambda p:miss(p)[1]-target,phase0,phase0+.1,xtol=1e-13)
    tm,dm,_=miss(phase)
    de=.5*(norm(y0[3:])**2-(V-u)**2)
    de_dot=np.dot([0,V,0],outgoing-incoming)
    # En el modelo de acercamiento de dos cuerpos, estimación del periapsis por el vector relativo.
    rrel=sol.sol(tm)[:3]-planet(tm,phase)
    vs=R3*n*np.array([-np.sin(phase+n*tm),np.cos(phase+n*tm),0.])
    vrel=sol.sol(tm)[3:]-vs
    E=.5*vrel@vrel-GM['S']/dm; ang=norm(np.cross(rrel,vrel))
    ecc=np.sqrt(1+2*E*ang*ang/GM['S']**2)
    estimated_rp=ang*ang/(GM['S']*(1+ecc))
    # Validación y ajuste adicional: integrar Sol+Saturno con fase prescrita.
    # No incluye Júpiter después del impulso asintótico idealizado.
    class Circular:
        def __init__(self,p): self.p=p
        def __call__(self,t,nu=0):
            angle=self.p+n*t+nu*np.pi/2
            return R3*n**nu*np.array([np.cos(angle),np.sin(angle),0.])
    def perturbed(p):
        sy=System({'S':Circular(p)})
        def collision(t,y): return norm(y[:3]-sy.planets['S'](t))-RAD['S']
        collision.terminal=True;collision.direction=-1
        integration=solve_ivp(lambda t,y:sy.rhs(t,y),(0,tc+250*DAY),y0,
                              events=collision,rtol=1e-12,atol=1e-8,
                              max_step=5*DAY,dense_output=True)
        if not integration.success: raise RuntimeError(integration.message)
        if len(integration.t_events[0]):
            return RAD['S'],integration.t[-1],integration.sol,sy
        mt,dist,_=minimum(integration.sol,sy.planets['S'],tc-250*DAY,tc+250*DAY)
        return dist,mt,integration.sol,sy
    # La rama de fase positiva evita colisión en las evaluaciones iniciales.
    center=minimize_scalar(lambda p:perturbed(p)[0],bounds=(phase0-.003,phase0+.003),method='bounded',options={'xatol':1e-9})
    phase_new=brentq(lambda p:perturbed(p)[0]-target,center.x,phase0+.1,xtol=1e-12)
    dpert,tpert,traj,sy=perturbed(phase_new)
    if dpert<RAD['S'] or rp<2*RAD['J']: raise RuntimeError('Diseño no admisible')
    result=dict(hypothesis='circular coplanar; Saturno fase libre; asistencia instantánea asintótica',
                rp_jupiter_km=rp,rp_min_allowed_km=2*RAD['J'],delta_deg=np.degrees(delta),
                vinf_in_kms=norm(incoming),vinf_out_kms=norm(outgoing),
                energy_gain_km2s2=de,energy_identity_residual=de-de_dot,
                saturn_phase_at_jupiter_deg=np.degrees(phase),time_to_saturn_days=tm/DAY,
                miss_solar_km=dm,estimated_saturn_rp_km=estimated_rp,
                saturn_radius_km=RAD['S'],initial_phase_miss_km=miss(phase+.2)[1],
                earth_launch_phase_saturn_deg=np.degrees(phase-n*h['tof_days']*DAY),
                earth_jupiter_tof_days=h['tof_days'],
                refined_phase_deg=np.degrees(phase_new),refined_saturn_rp_km=dpert,
                refined_time_days=tpert/DAY,
                refined_saturn_inside_soi=bool(dpert<sy.soi('S',tpert)),
                refined_model='Impulso J asintótico + integración Sol/Saturno prescrita')
    if return_refined:
        return result,traj,sy,tpert
    return result,sol,lambda t:planet(t,phase),lambda t:planet(t,phase+.2),tm
