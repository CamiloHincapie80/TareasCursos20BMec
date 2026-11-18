"""Tres trayectorias sobre un marco solar común, dirigidas a Saturno."""
import hashlib
import json
from pathlib import Path
from tempfile import NamedTemporaryFile
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from .animacion_datos import CACHE,load_ephemerides,sample_grid,voyager_positions,signature
from .ajuste_saturno import aimed_models
from .constantes import DAY,GM,MU,RAD
from .dinamica import System, grav
from .encuentros import calendar, encounter, jd


# GM de la Tierra en la cabecera de datos/horizons/tierra.txt, en km³/s².
EARTH_SOI_RATIO = (398600.435436 / MU)**.4


def earth_departures(a, epoch, earth, trajectories):
    """Primer cruce saliente en la ventana inicial de tres días."""
    departures = {}
    grid = np.linspace(a, a + 3*DAY, 289)
    for key, position in trajectories.items():
        def gap(t):
            center = earth(t)
            return np.linalg.norm(position(t) - center) - EARTH_SOI_RATIO*np.linalg.norm(center)

        distance = float(np.linalg.norm(position(a) - earth(a)))
        if not np.isfinite(distance):
            departures[key] = {'status':'no_data','initial_distance_km':None,
                               'elapsed_seconds':None,'date_tdb':None}
            continue
        entry = {'initial_distance_km': distance, 'elapsed_seconds': None, 'date_tdb': None}
        departures[key] = entry
        if gap(a) >= 0:
            entry['status'] = 'outside_at_start'
            continue
        entry['status'] = 'no_crossing'
        values = [gap(t) for t in grid]
        for lo, hi, left, right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if left < 0 <= right:
                t = brentq(gap, lo, hi, xtol=.001)
                entry.update(status='crossing', elapsed_seconds=t-a, date_tdb=calendar(epoch+t/DAY))
                break
    return departures


# 14:29 UTC (cabecera Voyager de Horizons) corresponde aproximadamente a 14:29:48 TDB.
LAUNCH_JD = jd('1977-08-20T14:29:48')


def extend_planets_to_launch(spl, launch):
    """Prolongación solar de dos cuerpos hasta el primer estado tabulado."""
    for key in ('T','J','S'):
        original = spl[key]
        start = original.x[0]
        initial = np.r_[original(start), original(start,1)]
        prefix = np.linspace(launch,start,65)[:-1]
        arc = solve_ivp(lambda t,y: np.r_[y[3:],grav(MU,y[:3])],(start,launch),
                        initial,t_eval=prefix[::-1],method='DOP853',rtol=1e-12,atol=1e-8)
        if not arc.success:
            raise RuntimeError(arc.message)
        times = np.r_[prefix,original.x]
        states = arc.y[:,::-1].T
        spl[key] = CubicHermiteSpline(times,np.vstack((states[:,:3],original(original.x))),
                                     np.vstack((states[:,3:],original(original.x,1))),extrapolate=False)


def comparison(f=1.,rebuild=False,*,launch=False):
    if f not in (.1,.25,.5,1.):raise ValueError('Factor permitido: 0.1, 0.25, 0.5 o 1')
    stamp={'source':signature('comparacion',f),
           'comparison_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'historical_launch':launch,
           'targeting_sha256':hashlib.sha256(Path(__file__).with_name('ajuste_saturno.py').read_bytes()).hexdigest()}
    CACHE.mkdir(parents=True,exist_ok=True);suffix='_lanzamiento' if launch else ''
    path=CACHE/f'comparacion_f{f:g}{suffix}.npz'
    if path.exists() and not rebuild:
        with np.load(path,allow_pickle=False) as saved:
            if json.loads(str(saved['signature']))==stamp:
                return {key:saved[key] for key in saved.files}
    epoch,spl=load_ephemerides();a=float(spl['T'].x[0])
    if launch:
        a=(LAUNCH_JD-epoch)*DAY
        extend_planets_to_launch(spl,a)
    print('Ajustando salida terrestre, encuentro joviano y llegada a Saturno…',flush=True)
    h,k,system,design=aimed_models(spl,a,f);b=design['end'];duration=b-a
    summaries={};saturn={}
    for key,trajectory in [('Hohmann',h),('K',k)]:
        summaries[key]=encounter(trajectory,system,'J',a,design['saturn_target_time'],epoch)
        saturn[key]=encounter(trajectory,system,'S',design['jupiter_time'],b,epoch)
        rp=saturn[key]['rp_km'];t=(saturn[key]['jd_tdb']-epoch)*DAY
        inside=rp<system.soi('S',t,f)
        saturn[key]['inside_soi']=bool(inside)
        if saturn[key]['boundary_minimum'] or not inside or rp<RAD['S']:
            raise RuntimeError(f'{key}: el ajuste no produce un sobrevuelo exterior dentro de la SOI de Saturno')
        if summaries[key]['rp_km']<RAD['J']:
            raise RuntimeError(f'{key}: la trayectoria ajustada intersecta Júpiter')
    entered_jupiter=any(event['frame_after']=='J' for event in k.jumps)
    if not entered_jupiter:raise RuntimeError('Cónicas empalmadas no alcanza la esfera joviana')
    summaries['K']['description']='Salida terrestre ajustada · cónicas empalmadas · δ jovicéntrico osculante.'
    summaries['Hohmann']['description']='Transferencia ajustada con dispersión gravitatoria continua Sol+Júpiter+Saturno. δ jovicéntrico osculante en periapsis.'
    vs=encounter(lambda t:np.r_[spl['VH'](t),spl['VH'](t,1)],System({'J':spl['JH']}),'J',spl['VH'].x[0],spl['VH'].x[-1],epoch)
    vs['description']='Datos originales · δ jovicéntrico osculante. La rapidez en Saturno está limitada por el muestreo diario.'
    summaries['Voyager']=vs
    # La referencia diaria de Saturno se etiqueta como estimada.
    vstate=lambda t:np.r_[spl['V'](t),spl['V'](t,1)]
    saturn['Voyager']=encounter(vstate,system,'S',max(a,spl['V'].x[0]),b,epoch)
    for summary in saturn.values():summary['elapsed_days']=(summary['jd_tdb']-epoch-a/DAY)
    departures=earth_departures(a,epoch,spl['T'],{
        'Hohmann':lambda t:h(t)[:3], 'K':lambda t:k(t)[:3], 'Voyager':spl['V']})
    events=[a+x['elapsed_seconds'] for x in departures.values() if x['elapsed_seconds'] is not None]
    events += [event['t'] for event in k.jumps]+[design['jupiter_time'],design['Hohmann_jupiter_time']]
    minima=[(summary['jd_tdb']-epoch)*DAY for summary in list(summaries.values())+list(saturn.values())]
    dense=[np.arange(a,a+3*DAY,600.)]
    dense += [np.arange(max(a,t-20*DAY),min(b,t+20*DAY),3600.) for t in minima]
    # Cinco minutos cerca del periapsis: evita facetas al ampliar el marco joviano.
    ht=design['Hohmann_jupiter_time']
    dense.append(np.arange(max(a,ht-2*DAY),min(b,ht+2*DAY),300.))
    absolute=np.unique(np.concatenate([sample_grid(a,b,epoch),np.array([a,b,spl['V'].x[0]]),
                       np.array(events),np.array(events)+1.,np.array(minima)]+dense))
    absolute=absolute[(absolute>=a)&(absolute<=b)];times=absolute-a
    hs=h(absolute);ks=k(absolute)
    result={'t':times,'Hohmann':hs[:,:3],'v_Hohmann':hs[:,3:],
            'K':ks[:,:3],'v_K':ks[:,3:],
            'Voyager':voyager_positions(spl,absolute),'v_Voyager':voyager_positions(spl,absolute,1)}
    for key in ('T','J','S'):result[key]=spl[key](absolute)
    result['soi_T']=EARTH_SOI_RATIO*np.linalg.norm(result['T'],axis=1)
    for key in ('J','S'):result['soi_'+key]=f*np.linalg.norm(result[key],axis=1)*(GM[key]/MU)**.4
    metadata={'epoch_jd':epoch+a/DAY,'factor':f,'initial_state_K':design['initial_state_K'],
              'initial_state_Hohmann':design['initial_state_Hohmann'],
              'launch_reference':'Ambos modelos salen de la posición terrestre; sus velocidades se ajustan por separado para alcanzar Saturno histórico.',
              'earth_departures':departures,'targeting':design,'K_events':k.jumps,'K_enters_J':entered_jupiter,
              'jupiter_encounters':summaries,'saturn_encounters':saturn,
              'comparison_note':('Salida: 20 ago. 1977. Posiciones planetarias prolongadas con dinámica solar hasta el 5 sep.; Voyager se dibuja desde el 5 sep. según Horizons. ' if launch else 'Simulación de 10 años desde 05 sep 1977, con planetas históricos. ')+
                  'Ambos modelos se ajustan para sobrevolar Saturno: '
                  'esferas usa cónicas empalmadas; Hohmann ajustada propaga la gravedad del Sol, Júpiter y Saturno sin impulsos (no es Hohmann circular pura). '
                  'Voyager conserva sus datos originales. Velocidad y energía son heliocéntricas; mínimos al centro y ángulos corresponden al encuentro completo. '
                  'Saturno en Voyager: estimación con datos diarios. ε = v²/2 − μ☉/r.'}
    result['metadata']=np.array(json.dumps(metadata,ensure_ascii=False));result['signature']=np.array(json.dumps(stamp,sort_keys=True))
    with NamedTemporaryFile(dir=CACHE,suffix='.npz',delete=False) as tmp:
        temporary=Path(tmp.name)
    try:
        np.savez_compressed(temporary,**result);temporary.replace(path)
    finally:temporary.unlink(missing_ok=True)
    return result
