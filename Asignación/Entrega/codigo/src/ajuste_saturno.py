"""Diseño dirigido a las efemérides: cónicas empalmadas y dispersión gravitatoria continua.
La segunda es una extensión ajustada de la referencia Hohmann, no una
transferencia Hohmann circular pura. No modifica las efemérides.
"""
import numpy as np
from scipy.optimize import brentq,least_squares
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
from .constantes import MU,GM,RAD,DAY
from .dinamica import System,invariants
from .encuentros import minimum


def stumpff(z):
    if abs(z)<1e-7:return .5-z/24+z*z/720,1/6-z/120+z*z/5040
    q=np.sqrt(abs(z))
    if z>0:return (1-np.cos(q))/z,(q-np.sin(q))/q**3
    return (np.cosh(q)-1)/(-z),(np.sinh(q)-q)/q**3


def lambert(r1,r2,dt):
    """Solución prógrada de cero revoluciones, km y segundos."""
    if dt<=0:raise ValueError('Tiempo de transferencia no positivo')
    n1=np.linalg.norm(r1);n2=np.linalg.norm(r2)
    cosine=np.clip(r1@r2/(n1*n2),-1,1)
    sine=np.sqrt(1-cosine*cosine)*np.sign(np.cross(r1,r2)[2])
    if abs(sine)<1e-10:raise ValueError('Geometría Lambert colineal')
    A=sine*np.sqrt(n1*n2/(1-cosine))
    def residual(z):
        C,S=stumpff(z);y=n1+n2+A*(z*S-1)/np.sqrt(C)
        if y<0:return np.nan
        return (y/C)**1.5*S+A*np.sqrt(y)-np.sqrt(MU)*dt
    previous=None
    for z in np.linspace(-4*np.pi**2,4*np.pi**2-1e-5,400):
        value=residual(z)
        if not np.isfinite(value):continue
        if previous is not None and previous[1]*value<=0:
            root=brentq(residual,previous[0],z,xtol=1e-12);break
        previous=(z,value)
    else:raise RuntimeError('No se encontró una transferencia Lambert de cero revoluciones')
    C,S=stumpff(root);y=n1+n2+A*(root*S-1)/np.sqrt(C)
    ff=1-y/n1;gg=A*np.sqrt(y/MU);gd=1-y/n2
    return (r2-ff*r1)/gg,(gd*r2-r1)/gg


def aimed_models(spl,a,f):
    earth=spl['T'](a)
    def legs(x):
        tj=a+x[0]*DAY;ts=tj+x[1]*DAY
        v0,vj=lambert(earth,spl['J'](tj),tj-a)
        vo,_=lambert(spl['J'](tj),spl['S'](ts),ts-tj)
        u=vj-spl['J'](tj,1);w=vo-spl['J'](tj,1)
        un=np.linalg.norm(u);wn=np.linalg.norm(w)
        delta=np.arccos(np.clip(u@w/un/wn,-1,1))
        rp=GM['J']/un**2*(1/np.sin(delta/2)-1)
        return v0,u,w,rp,un-wn
    def objective(x):
        *_,rp,difference=legs(x)
        return [difference,(rp-10*RAD['J'])/1e6]
    dates=least_squares(objective,[670,780],bounds=([450,400],[1100,1800]),
                        xtol=1e-11,ftol=1e-11,gtol=1e-11)
    if np.linalg.norm(dates.fun)>1e-5:raise RuntimeError('No convergieron las fechas del diseño joviano')
    tj=a+dates.x[0]*DAY;ts=tj+dates.x[1]*DAY
    # Dos años después del último encuentro previsto; ambos estados continúan integrándose.
    end=max(ts+2*365.25*DAY,a+10*365.25*DAY)
    if any(end>spl[key].x[-1] for key in ('T','J','S','V')):
        raise ValueError('El diseño excede la cobertura de las efemérides')
    _,u,w,rp,_=legs(dates.x);un=np.linalg.norm(u)
    perpendicular=w-u*(w@u)/(u@u)
    impact=-perpendicular/np.linalg.norm(perpendicular)*rp*np.sqrt(1+2*GM['J']/rp/un**2)
    seed,_=lambert(earth,spl['J'](tj)+impact,tj-a)
    # Un punto a 1 millón de km sobre el plano de referencia saturniano evita apuntar al centro.
    target=spl['S'](ts)+np.array([0.,0.,1e6])
    system=System({key:spl[key] for key in ('J','S')})
    settings=dict(rtol=1e-11,atol=1e-7,max_step=3*DAY)
    def shoot_k(velocity):
        trajectory=system.run(np.r_[earth,velocity],a,ts,model='K',f=f,**settings)
        return (trajectory(ts)[:3]-target)/1e6
    fit=least_squares(shoot_k,seed,diff_step=1e-6,x_scale='jac',max_nfev=70,
                       xtol=1e-11,ftol=1e-10,gtol=1e-8)
    if np.linalg.norm(shoot_k(fit.x))>.01:raise RuntimeError('No convergió la salida de cónicas empalmadas hacia Saturno')
    initial=np.r_[earth,fit.x]
    k=system.run(initial,a,end,model='K',f=f,**settings)

    # Disparo desde periapsis en ambos sentidos: evita la sensibilidad extrema
    # de apuntar a Saturno variando sólo la velocidad desde la Tierra.
    kt,_,_=minimum(k,spl['J'],tj-80*DAY,tj+80*DAY)
    relative=k(kt)-system.state('J',kt)
    radial=relative[:3]/np.linalg.norm(relative[:3])
    tangent=relative[3:]-radial*(relative[3:]@radial)
    tangent/=np.linalg.norm(tangent)
    rp_target=10*RAD['J']
    def peri_state(parameters):
        rotation=Rotation.from_rotvec(parameters[:3])
        peri=kt+parameters[4]*DAY
        speed=np.sqrt(parameters[3]**2+2*GM['J']/rp_target)
        return peri,system.state('J',peri)+np.r_[rp_target*rotation.apply(radial),
                                                         speed*rotation.apply(tangent)]
    def integrate_peri(parameters,finish):
        peri,state=peri_state(parameters)
        sol=solve_ivp(lambda t,y:system.rhs(t,y), (peri,finish),state,
                      method='DOP853',dense_output=True,**settings)
        if not sol.success:raise RuntimeError(sol.message)
        return sol
    def residual_h(parameters):
        arrival=ts+parameters[5]*DAY
        past=integrate_peri(parameters,a)
        future=integrate_peri(parameters,arrival)
        destination=spl['S'](arrival)+np.array([0.,0.,1e6])
        return np.r_[(past.y[:3,-1]-earth)/1e6,
                     (future.y[:3,-1]-destination)/1e6]
    def jacobian_h(parameters):
        # Diferencias centrales con pasos físicos absolutos. El paso relativo
        # automático se anula en los offsets de fecha inicialmente nulos.
        steps=np.array([1e-5,1e-5,1e-5,1e-4,.001,.001])
        columns=[]
        for i,step in enumerate(steps):
            offset=np.zeros(6);offset[i]=step
            columns.append((residual_h(parameters+offset)-residual_h(parameters-offset))/(2*step))
        return np.column_stack(columns)
    seed_vinf=invariants(relative,GM['J'])['vinf_kms']
    corrected=least_squares(residual_h,[0.,0.,0.,seed_vinf,0.,0.],
                            bounds=([-1.,-1.,-1.,.1,-60.,-60.],[1.,1.,1.,30.,60.,60.]),
                            jac=jacobian_h,x_scale=[.1,.1,.1,5.,10.,10.],
                            max_nfev=100,xtol=1e-12,ftol=1e-12,gtol=1e-9)
    if np.linalg.norm(residual_h(corrected.x))>1e-5:
        raise RuntimeError(f'No convergió la dispersión continua hacia Saturno: {corrected.fun}; {corrected.x}; {corrected.message}; {corrected.optimality}')
    ht,_=peri_state(corrected.x)
    arrival=ts+corrected.x[5]*DAY
    past=integrate_peri(corrected.x,a)
    # Repropagar desde la posición terrestre exacta con la velocidad convergida.
    # La trayectoria entregada es una única solución, sin empalmes ni impulsos.
    def forward_residual(parameters):
        time=ts+parameters[3]*DAY
        trajectory=system.run(np.r_[earth,parameters[:3]],a,time,**settings)
        _,radius,_=minimum(trajectory,spl['J'],ht-DAY,ht+DAY)
        destination=spl['S'](time)+np.array([0.,0.,1e6])
        return np.r_[(trajectory(time)[:3]-destination)/1e6,(radius-rp_target)/1e6]
    def forward_jacobian(parameters):
        columns=[]
        for i,step in enumerate([1e-7,1e-7,1e-7,.001]):
            offset=np.zeros(4);offset[i]=step
            columns.append((forward_residual(parameters+offset)-forward_residual(parameters-offset))/(2*step))
        return np.column_stack(columns)
    polished=least_squares(forward_residual,np.r_[past.y[3:,-1],corrected.x[5]],
                           jac=forward_jacobian,x_scale=[10.,10.,10.,10.],
                           max_nfev=20,xtol=1e-12,ftol=1e-12,gtol=1e-8)
    if np.linalg.norm(forward_residual(polished.x))>1e-5:
        raise RuntimeError('No convergió la verificación heliocéntrica del sobrevuelo')
    arrival=ts+polished.x[3]*DAY
    h=system.run(np.r_[earth,polished.x[:3]],a,end,**settings)
    ht,_,_=minimum(h,spl['J'],ht-DAY,ht+DAY)
    hrp=np.linalg.norm(h(ht)[:3]-spl['J'](ht))
    relative=h(ht)-system.state('J',ht)
    scattering=invariants(relative,GM['J'])
    vinf=scattering['vinf_kms'];ecc=scattering['e']
    semimajor=-GM['J']/vinf**2
    destination=spl['S'](arrival)+np.array([0.,0.,1e6])
    design={'jupiter_time':float(tj),'saturn_target_time':float(ts),'end':float(end),
            'initial_state_K':initial.tolist(),'initial_state_Hohmann':h(a).tolist(),
            'K_target_error_km':float(np.linalg.norm(k(ts)[:3]-target)),
            'Hohmann_target_error_km':float(np.linalg.norm(h(arrival)[:3]-destination)),
            'Hohmann_saturn_target_time':float(arrival),
            'Hohmann_jupiter_time':float(ht),'Hohmann_rp_km':float(hrp),
            'Hohmann_turn_deg':scattering['delta_deg'],
            'Hohmann_hyperbola':{'a_km':float(semimajor),'e':ecc,
                'rp_km':float(semimajor*(1-ecc)),
                'b_km':float(abs(semimajor)*np.sqrt(ecc**2-1)),
                'vinf_kms':vinf,'interpretation':'Elementos jovicéntricos osculantes en periapsis; incluye perturbaciones solares.'},
            'Hohmann_propagation':'Gravedad continua Sol+Júpiter+Saturno, sin impulsos.',
            'K_events':k.jumps}
    return h,k,system,design
