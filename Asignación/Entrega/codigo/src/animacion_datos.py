"""Muestreo de las soluciones existentes para VisPy; no altera sus ecuaciones.
La caché valida hashes de datos/código físico y nunca extrapola efemérides.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp
from .horizons import read,compatible
from .constantes import AU,DAY,MU,GM
from .dinamica import System,grav
from .encuentros import jd,crossings,minimum,encounter
from .hohmann import calculate,R1,R2,R3
from .optimizacion import design

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'resultados'/'animaciones'
SCHEMA=3
HORIZONS=Path('datos')/'horizons'
DATA_FILES=['voyager_2.txt','Voyager_jupiter.txt','Jupiter1979.txt','tierra.txt','jupiter.txt','saturno.txt']
FILES=[*(HORIZONS/f for f in DATA_FILES),
       'src/horizons.py','src/dinamica.py','src/hohmann.py','src/optimizacion.py','src/constantes.py',
       'src/animacion_datos.py','src/encuentros.py']

def signature(mode,f):
    hashes={str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
    return {'schema':SCHEMA,'mode':mode,'factor':f,'sha256':hashes}

def load_ephemerides():
    data={k:read(ROOT/HORIZONS/f) for k,f in [('V','voyager_2.txt'),('T','tierra.txt'),('J','jupiter.txt'),('S','saturno.txt'),
                                    ('VH','Voyager_jupiter.txt'),('JH','Jupiter1979.txt')]}
    compatible(*data.values());epoch=data['V'].jd[0]
    return epoch,{k:v.spline(epoch) for k,v in data.items()}

def sample_grid(a,b,epoch):
    parts=[np.linspace(a,b,2400)]
    for start,end in [('1979-07-01','1979-07-20'),('1981-08-10','1981-09-10')]:
        lo=max(a,(jd(start)-epoch)*DAY);hi=min(b,(jd(end)-epoch)*DAY)
        if lo<hi:parts.append(np.arange(lo,hi,3600.))
    return np.unique(np.concatenate(parts+[np.array([a,b])]))

def voyager_positions(spl,times,derivative=0):
    """Utiliza la serie horaria donde existe; el resto sigue siendo diario."""
    y=spl['V'](times,derivative)
    hourly=spl['VH']; mask=(times>=hourly.x[0])&(times<=hourly.x[-1])
    y[mask]=hourly(times[mask],derivative)
    return y

def history(f):
    epoch,spl=load_ephemerides()
    a=(jd('1978-11-01')-epoch)*DAY;b=(jd('1981-12-01')-epoch)*DAY
    sys=System({k:spl[k] for k in ('J','S')})
    times=sample_grid(a,b,epoch); y0=np.r_[spl['V'](a),spl['V'](a,1)]
    print('Integrando N hasta diciembre de 1981…',flush=True)
    n=sys.run(y0,a,b)
    tin=crossings(n,sys,'J',a,(jd('1979-07-10')-epoch)*DAY,f)[0]
    # El punto exacto de entrada queda en la rejilla para mostrar la continuidad.
    peri_J=minimum(lambda t:spl['VH'](t),spl['JH'],spl['VH'].x[0],spl['VH'].x[-1])[0]
    peri_S=minimum(lambda t:spl['V'](t),spl['S'],(jd('1981-08-01')-epoch)*DAY,(jd('1981-09-20')-epoch)*DAY)[0]
    times=np.unique(np.r_[times,tin,peri_J,peri_S]);mask=times>=tin
    result={'t':times-times[0],'jd':epoch+times/DAY,'Voyager':voyager_positions(spl,times),
            'N':n(times)[:,:3], 'v_N':n(times)[:,3:],
            'v_Voyager':voyager_positions(spl,times,1)}
    for k in ('T','J','S'):result[k]=spl[k](times)
    metadata={'epoch_jd':float(epoch+times[0]/DAY),'factor':f,'entry_jd':float(epoch+tin/DAY),
              'shared_prefix':'P/P0/K comparten N antes de entrada J; no son integraciones propias en ese tramo.',
              'voyager':'Horizons original aproximado: horario en J (1–20 julio 1979), diario fuera.',
              'events':{},'frames':{},
              'jupiter_encounters':{'N':encounter(n,sys,'J',a,b,epoch)}}
    for model in ('P','P0','K'):
        print(f'Integrando {model}, f={f}…',flush=True)
        sol=sys.run(n(tin)-sys.state('J',tin),tin,b,model,f,frame='J')
        positions=result['N'].copy();positions[mask]=sol(times[mask])[:,:3]
        result[model]=positions
        velocities=result['v_N'].copy();velocities[mask]=sol(times[mask])[:,3:]
        result['v_'+model]=velocities
        def full_state(t):
            return n(t) if t<tin else sol(t)
        metadata['jupiter_encounters'][model]=encounter(full_state,sys,'J',a,b,epoch)
        metadata['events'][model]=[{**event,'jd':epoch+event['t']/DAY} for event in sol.jumps]
        # Estados de región sobre la rejilla; N siempre H, P cambia solo representación.
        frame=np.zeros(len(times),dtype=np.int8)
        for fr,segment in sol.segments:
            flag=(times>=segment.t[0])&(times<=segment.t[-1]);frame[flag]={'H':0,'J':1,'S':2}[fr]
        result['region_'+model]=frame
    for k in ('J','S'):
        result['soi_'+k]=f*np.linalg.norm(result[k],axis=1)*(GM[k]/MU)**.4
    # Saltos rápidos usan los mínimos de Voyager; S diario se marca como insuficiente.
    metadata['peri_J']=float(peri_J-times[0]);metadata['peri_S']=float(peri_S-times[0])
    result['metadata']=np.array(json.dumps(metadata,ensure_ascii=False))
    return result

def hohmann():
    epoch,spl=load_ephemerides();h=calculate();th=h['tof_days']*DAY
    print('Calculando Hohmann y extensión idealizada hacia Saturno…',flush=True)
    design_info,post,system,tperi=design(return_refined=True)
    end=th+tperi+100*DAY
    times=np.unique(np.r_[np.linspace(0,end,3000),th,th+tperi])
    initial=np.array([R1,0,0,0,h['vt1_kms'],0])
    transfer=solve_ivp(lambda t,y:np.r_[y[3:],grav(MU,y[:3])],(0,th),initial,
                       method='DOP853',rtol=1e-12,atol=1e-8,dense_output=True)
    if not transfer.success:raise RuntimeError(transfer.message)
    before=times<=th;pos=np.empty((len(times),3));vel=np.empty_like(pos)
    pos[before]=transfer.sol(times[before])[:3].T
    vel[before]=transfer.sol(times[before])[3:].T
    # El diseño calcula J en +x; girarlo π lo conecta con el afelio Hohmann en −x.
    # Es rotación de la geometría idealizada, nunca de los datos históricos.
    pos[~before]=np.asarray([post(t-th)[:3] for t in times[~before]])*[-1,-1,1]
    vel[~before]=np.asarray([post(t-th)[3:] for t in times[~before]])*[-1,-1,1]
    phaseJ=np.radians(h['phase_deg'])
    phaseS=np.pi+np.radians(design_info['refined_phase_deg'])-np.sqrt(MU/R3**3)*th
    def circle(R,phase):
        angle=phase+np.sqrt(MU/R**3)*times
        return R*np.column_stack((np.cos(angle),np.sin(angle),np.zeros_like(angle)))
    hist_a=(jd('1977-09-05')-epoch)*DAY;hist_b=(jd('1981-12-01')-epoch)*DAY
    ht=sample_grid(hist_a,hist_b,epoch)
    result={'t':times,'Hohmann':pos,'v_Hohmann':vel,'T':circle(R1,0),'J':circle(R2,phaseJ),'S':circle(R3,phaseS),
            'hist_t':ht-ht[0],'hist_jd':epoch+ht/DAY,'hist_Voyager':voyager_positions(spl,ht),'hist_v_Voyager':voyager_positions(spl,ht,1)}
    for k in ('T','J','S'):result['hist_'+k]=spl[k](ht)
    metadata={'transfer_days':h['tof_days'],'duration_days':end/DAY,'peri_J':th,'peri_S':th+tperi,
              'historical_start_jd':epoch+hist_a/DAY,'historical_end_jd':epoch+hist_b/DAY,
              'design':design_info,'hohmann':h,
              'jupiter_encounters':{
                  'Hohmann':{'rp_km':float(design_info['rp_jupiter_km']),
                             'delta_deg':float(design_info['delta_deg']),
                             'boundary_minimum':False},
                  'Voyager':encounter(lambda t:np.r_[spl['VH'](t),spl['VH'](t,1)],
                                      System({'J':spl['JH']}),'J',
                                      spl['VH'].x[0],spl['VH'].x[-1],epoch)},
              'clock':'Los paneles comparten avance porcentual; no comparten fecha ni geometría.',
              'initial':'La vista solar empieza tras el escape asintótico de LEO 500 km; no integra la órbita LEO.',
              'extension':'Tras Hohmann: giro joviano asintótico instantáneo y Sol+Saturno circular prescrito.'}
    result['metadata']=np.array(json.dumps(metadata,ensure_ascii=False,default=lambda x:x.item()))
    return result

def get_data(mode,f=1.,rebuild=False):
    if mode not in ('hohmann','esferas'):raise ValueError('Modo no reconocido')
    if f not in (.1,.25,.5,1.):raise ValueError('Factor permitido: 0.1, 0.25, 0.5 o 1')
    CACHE.mkdir(parents=True,exist_ok=True)
    path=CACHE/(f'{mode}_f{f:g}.npz' if mode=='esferas' else 'hohmann.npz')
    stamp=signature(mode,f if mode=='esferas' else 1.)
    if path.exists() and not rebuild:
        with np.load(path,allow_pickle=False) as archive:
            if json.loads(str(archive['signature']))==stamp:
                return {k:archive[k] for k in archive.files}
    data=history(f) if mode=='esferas' else hohmann()
    data['signature']=np.array(json.dumps(stamp,sort_keys=True))
    # No dejar una caché parcial si se interrumpe el cálculo.
    tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**data);tmp.replace(path)
    return data
