"""Experimentos adicionales: errores de datos, origen y física omitida."""
import numpy as np
from .horizons import read
from .dinamica import System,grav,norm
from .constantes import MU,GM,DAY
from .encuentros import jd,minimum,crossings

def experiments(root):
    eph={k:read(root/f) for k,f in [('V','voyager_2.txt'),('J','jupiter.txt'),('S','saturno.txt')]}
    epoch=eph['V'].jd[0];T=lambda date:(jd(date)-epoch)*DAY
    planets={k:eph[k].spline(epoch) for k in ('J','S')};sys=System(planets)
    v=eph['V'].spline(epoch);state=lambda t:np.r_[v(t),v(t,1)]
    hourly=read(root/'Voyager_jupiter.txt');hq=(hourly.jd-epoch)*DAY
    local=(hq>=T('1979-07-08'))&(hq<=T('1979-07-12'))
    diff=np.linalg.norm(v(hq)-hourly.states[:,:3],axis=1)
    rows=[dict(experiment='interpolacion_diaria_contra_horaria_J',
               max_dr_km=float(diff[local].max()),note='Muestras horarias no usadas en Hermite diario')]
    # Aceleraciones descompuestas: solar, otro planeta, inconsistencia del origen.
    t0=T('1978-11-01');tf=T('1981-10-01');n=sys.run(state(t0),t0,tf)
    tin=crossings(n,sys,'J',t0,T('1979-07-10'))[0]
    tout=crossings(n,sys,'J',T('1979-07-10'),T('1980-01-01'))[-1]
    tp,_,_=minimum(n,planets['J'],T('1979-07-08'),T('1979-07-12'))
    for label,t in [('entrada',tin),('periapsis',tp),('salida',tout)]:
        R=planets['J'](t);r=n(t)[:3];rho=r-R
        tide=grav(MU,r)-grav(MU,R)
        other=grav(GM['S'],r-planets['S'](t))-grav(GM['S'],R-planets['S'](t))
        modeled=grav(MU+GM['J'],R)+grav(GM['S'],R-planets['S'](t))+grav(GM['S'],planets['S'](t))
        residual=modeled-planets['J'](t,2)
        denom=norm(grav(GM['J'],rho))
        rows.append(dict(experiment='descomposicion_'+label,eta_solar=norm(tide)/denom,
                         eta_saturn=norm(other)/denom,eta_origin=norm(residual)/denom))
    # Separar omisión de Saturno de la transformación de coordenadas.
    a=T('1979-06-01');b=T('1979-08-01');y=state(a)
    base=sys.run(y,a,b)
    for label,sy in [('sin_S',System({'J':planets['J']})),
                     ('planetas_cada_dos_dias',System({k:eph[k].spline(epoch,stride=2) for k in ('J','S')}))]:
        sol=sy.run(y,a,b)
        rows.append(dict(experiment=label,dr_km=norm(sol(b)[:3]-base(b)[:3]),dv_kms=norm(sol(b)[3:]-base(b)[3:])))
    # Convergencia global P-N con mismas CI en f=1 y control del paso solar.
    for tol,step in [(1e-11,DAY),(1e-12,DAY/2),(1e-13,DAY/4),(2.3e-14,DAY/8)]:
        ref=sys.run(state(t0),t0,tf,rtol=tol,atol=tol*1e4,max_step=step)
        p=sys.run(ref(tin)-sys.state('J',tin),tin,tf,'P',frame='J',rtol=tol,atol=tol*1e4,max_step=step)
        rows.append(dict(experiment='convergencia_global_P_N',rtol=tol,max_step_s=step,
                         dr_km=norm(p(tf)[:3]-ref(tf)[:3]),dv_kms=norm(p(tf)[3:]-ref(tf)[3:])))
    # Invariantes K solo mientras permanece en región joviana.
    k=sys.run(n(tin)-sys.state('J',tin),tin,tout,'K',frame='J')
    fr,segment=k.segments[0]
    from .dinamica import invariants
    vals=[invariants(segment.sol(t),GM['J']) for t in np.linspace(segment.t[0],segment.t[-1],500)]
    rows.append(dict(experiment='invariantes_K_J',
                     relative_energy_range=np.ptp([x['energy_km2s2'] for x in vals])/abs(vals[0]['energy_km2s2']),
                     relative_h_range=np.ptp([x['h_km2s'] for x in vals])/vals[0]['h_km2s']))
    # Medidas de velocidades finitas, sin llamarlas asintóticas.
    yi=n(tin);yf=n(tout);ui=yi[3:]-planets['J'](tin,1);uf=yf[3:]-planets['J'](tout,1)
    rows.append(dict(experiment='velocidades_frontera_N',
                     t_in_jd=epoch+tin/DAY,t_out_jd=epoch+tout/DAY,
                     u_in_kms=ui.tolist(),u_out_kms=uf.tolist(),
                     finite_angle_deg=np.degrees(np.arccos(ui@uf/(norm(ui)*norm(uf)))),
                     heliocentric_speed_in=norm(yi[3:]),heliocentric_speed_out=norm(yf[3:]),
                     energy_in=.5*yi[3:]@yi[3:]-MU/norm(yi[:3]),
                     energy_out=.5*yf[3:]@yf[3:]-MU/norm(yf[:3])))
    return rows

def saturn_exit_experiment(root):
    """Continuación independiente hasta diciembre para registrar también la salida S."""
    eph={k:read(root/f) for k,f in [('V','voyager_2.txt'),('J','jupiter.txt'),('S','saturno.txt')]}
    epoch=eph['V'].jd[0];T=lambda date:(jd(date)-epoch)*DAY
    sy=System({k:eph[k].spline(epoch) for k in ('J','S')})
    v=eph['V'].spline(epoch);t0=T('1978-11-01');tf=T('1981-12-01')
    n=sy.run(np.r_[v(t0),v(t0,1)],t0,tf)
    tin=crossings(n,sy,'J',t0,T('1979-07-10'))[0]
    rows=[]
    from .encuentros import calendar
    for model in ('P','P0','K'):
        sol=sy.run(n(tin)-sy.state('J',tin),tin,tf,model,frame='J')
        for j in sol.jumps:
            j['date_tdb']=calendar(epoch+j['t']/DAY)
        rows.append(dict(model=model,events=sol.jumps,final_frame=sol.segments[-1][0],
                         final_dr_vs_N_km=norm(sol(tf)[:3]-n(tf)[:3])))
    return rows
