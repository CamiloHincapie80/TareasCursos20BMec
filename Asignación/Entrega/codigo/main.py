"""Punto de entrada reproducible. Ejecutar con Python del entorno local."""
from pathlib import Path
import argparse
import json
import csv
import numpy as np
from src.horizons import read,compatible
from src.constantes import MU,GM,DAY
from src.dinamica import System,invariants,norm
from src.encuentros import jd,encounter,crossings
from src.hohmann import calculate
from src.optimizacion import design

ROOT=Path(__file__).resolve().parent
HORIZONS=ROOT/'datos'/'horizons'
OUT=ROOT/'resultados'/'tablas'
OUT.mkdir(parents=True,exist_ok=True)

def save(name,data):
    (OUT/(name+'.json')).write_text(json.dumps(data,indent=2,ensure_ascii=False,default=lambda x:x.item())+'\n')

def table(name,rows):
    if not rows:return
    with (OUT/(name+'.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def data_audit():
    data=[read(p) for p in sorted(HORIZONS.glob('*.txt'))]
    compatible(*data)
    audits=[e.audit() for e in data]
    save('efemerides',audits)
    processed=ROOT/'datos'/'procesados'
    processed.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(processed/'estados.npz',
                        **{e.path.stem:np.column_stack((e.jd,e.states)) for e in data})
    print('Datos:',[(e.path.name,len(e.jd)) for e in data],flush=True)

def load():
    eph={k:read(HORIZONS/f) for k,f in [('V','voyager_2.txt'),('J','jupiter.txt'),('S','saturno.txt')]}
    compatible(*eph.values()); epoch=eph['V'].jd[0]
    system=System({k:eph[k].spline(epoch) for k in ('J','S')})
    v=eph['V'].spline(epoch)
    state=lambda t:np.r_[v(t),v(t,1)]
    return eph,epoch,system,state

def history(compare=True,figures=False):
    eph,epoch,sys,reference=load()
    T=lambda date:(jd(date)-epoch)*DAY
    t0=T('1978-11-01'); tf=T('1981-10-01')
    y0=reference(t0)
    n=sys.run(y0,t0,tf)
    windows={'J':(T('1979-07-01'),T('1979-07-20')),'S':(T('1981-08-01'),T('1981-09-20'))}
    rows=[]
    hourly=read(HORIZONS/'Voyager_jupiter.txt'); hj=read(HORIZONS/'Jupiter1979.txt');compatible(hourly,hj)
    hv=hourly.spline(epoch); hsys=System({'J':hj.spline(epoch)})
    href=lambda t:np.r_[hv(t),hv(t,1)]
    for k,(a,b) in windows.items():
        er=encounter(href if k=='J' else reference,hsys if k=='J' else sys,k,a,b,epoch)
        er['model']='Efemerides';rows.append(er)
        er=encounter(n,sys,k,a,b,epoch);er['model']='N';rows.append(er)
    save('encuentros',rows)
    runs={'N':n}; comparison=[]
    if compare:
        for fac in (.1,.25,.5,1.):
            cr=crossings(n,sys,'J',t0,T('1979-07-10'),fac)
            if not cr: raise RuntimeError(f'N no entra a esfera J f={fac}')
            tin=cr[0]; yin=n(tin)-sys.state('J',tin)
            refJ=encounter(n,sys,'J',*windows['J'],epoch)
            refS=encounter(n,sys,'S',T('1981-01-01'),tf,epoch)
            for model in ('P','P0','K'):
                sol=sys.run(yin,tin,tf,model,fac,frame='J')
                y=sol(tf); delta=y-n(tf)
                eJ=encounter(sol,sys,'J',max(tin,windows['J'][0]),windows['J'][1],epoch)
                eS=encounter(sol,sys,'S',T('1981-01-01'),tf,epoch)
                comparison.append(dict(f=fac,model=model,dr_km=norm(delta[:3]),dv_kms=norm(delta[3:]),
                    dE_km2s2=invariants(y,MU)['energy_km2s2']-invariants(n(tf),MU)['energy_km2s2'],
                    drp_J_km=eJ['rp_km']-refJ['rp_km'],ddelta_deg=eJ['delta_deg']-refJ['delta_deg'],
                    rp_S_km=eS['rp_km'],drp_S_km=eS['rp_km']-refS['rp_km'],
                    dt_S_s=(eS['jd_tdb']-refS['jd_tdb'])*DAY,
                    minimum_at_boundary=eS['boundary_minimum'],
                    drp_S_relative=(eS['rp_km']-refS['rp_km'])/refS['rp_km'],
                    drp_J_relative=(eJ['rp_km']-refJ['rp_km'])/refJ['rp_km'],
                    dE_pre_saturn=invariants(sol(T('1981-06-01')),MU)['energy_km2s2']-invariants(n(T('1981-06-01')),MU)['energy_km2s2'],
                    saturn_inside_soi=eS['rp_km']<sys.soi('S',(eS['jd_tdb']-epoch)*DAY,fac),
                    equivalent_1000km=abs(eS['rp_km']-refS['rp_km'])<1000,
                    max_jump_km=max((a['dr_km'] for a in sol.jumps),default=0.)))
                print('Comparación',fac,model,comparison[-1]['dr_km'],flush=True)
                if fac==1.:
                    runs[model]=sol
                    save('eventos_'+model,dict(jumps=sol.jumps,segments=[dict(frame=fr,start_jd=epoch+ss.t[0]/DAY,end_jd=epoch+ss.t[-1]/DAY) for fr,ss in sol.segments]))
                    for k in ('J','S'):
                        er=encounter(sol,sys,k,*windows[k],epoch);er['model']=model;rows.append(er)
                save('comparacion',comparison);table('comparacion',comparison)
        save('encuentros',rows);table('encuentros',rows)
    # Comparación histórica: error de posición en la fecha del mínimo de la efeméride.
    for row in rows:
        k=row['planet'];truth=next(r for r in rows if r['planet']==k and r['model']=='Efemerides')
        qt=(truth['jd_tdb']-epoch)*DAY
        model_state=(href if k=='J' else reference) if row['model']=='Efemerides' else runs[row['model']]
        true_state=href if k=='J' else reference
        row['position_error_at_reference_periapsis_km']=norm(model_state(qt)[:3]-true_state(qt)[:3])
        row['periapsis_distance_error_km']=row['rp_km']-truth['rp_km']
        row['periapsis_distance_relative_error']=(row['rp_km']-truth['rp_km'])/truth['rp_km']
        row['periapsis_time_error_s']=(row['jd_tdb']-truth['jd_tdb'])*DAY
        row['reference_quality']='horaria aproximada' if k=='J' else 'diaria insuficiente: no certificar periapsis'
    save('encuentros',rows);table('encuentros',rows)
    table('comparacion_saturno',[r for r in rows if r['planet']=='S'])
    # Diagnóstico físico sobre la misma trayectoria N en tres instantes del encuentro.
    cs=crossings(n,sys,'J',t0,T('1980-01-01'))
    per=(rows[1]['jd_tdb']-epoch)*DAY
    diagnostics=[]
    for label,t in zip(('entrada','periapsis','salida'),(cs[0],per,cs[-1])):
        y=n(t)-sys.state('J',t)
        central=-GM['J']*y[:3]/norm(y[:3])**3
        perturb=sys.rhs(t,y,'P','J')[3:]-central
        diagnostics.append(dict(point=label,jd_tdb=epoch+t/DAY,eta=norm(perturb)/norm(central),
                                **invariants(y,GM['J'])))
    save('diagnostico_jupiter',diagnostics)
    if figures:
        from src.graficas import historical
        historical(sys,epoch,reference,n,runs,comparison,windows,cs)
    return sys,epoch,reference,n,runs,windows

def convergence():
    _,epoch,sys,reference=load(); T=lambda date:(jd(date)-epoch)*DAY
    # Mismo IC de N en entrada f=0.25; prueba más exigente que una órbita circular.
    t0=T('1979-06-01'); tf=T('1979-08-01');y0=reference(t0)
    sols=[]; rows=[]
    for tol,step,method in [(1e-9,DAY,'DOP853'),(1e-11,DAY/4,'DOP853'),(1e-12,DAY/8,'DOP853'),(1e-11,DAY/4,'Radau')]:
        sol=sys.run(y0,t0,tf,rtol=tol,atol=tol*100, max_step=step,method=method)
        sols.append(sol); rows.append(dict(rtol=tol,max_step_s=step,method=method))
    ref=sols[2](tf)
    for row,sol in zip(rows,sols):
        row['dr_km']=norm(sol(tf)[:3]-ref[:3]);row['dv_kms']=norm(sol(tf)[3:]-ref[3:])
    save('convergencia',rows);table('convergencia',rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('task',choices=['datos','hohmann','historia','comparacion','saturno','diseno','convergencia','figuras','validacion','todo'])
    a=p.parse_args().task
    if a in ('datos','todo'):data_audit()
    if a in ('hohmann','todo'):save('hohmann',calculate());print(calculate())
    if a in ('historia','comparacion','saturno','figuras','todo'):
        history(compare=a!='historia',figures=a in ('figuras','todo'))
    if a in ('diseno','todo','figuras'):
        d,*_=design();save('diseno',d);print(d)
        if a in ('todo','figuras'):
            from src.graficas import idealized
            idealized()
    if a in ('convergencia','todo'):convergence()

    if a in ('validacion','todo'):
        from src.validacion import experiments
        save('validacion_adicional',experiments(HORIZONS))
    if a in ('saturno','todo'):
        from src.validacion import saturn_exit_experiment
        save('salida_saturno',saturn_exit_experiment(HORIZONS))
