"""N: partícula de prueba heliocéntrica; P: misma física en marcos trasladados.
P0 conserva la ecuación planetocéntrica original; K conserva solo Kepler.
"""
import numpy as np
from scipy.integrate import solve_ivp
from .constantes import MU, GM, DAY
norm=np.linalg.norm

def grav(mu,r):
    return -mu*r/norm(r)**3

class System:
    def __init__(self, planets): self.planets=planets
    def state(self,k,t):
        s=self.planets[k]
        return np.r_[s(t),s(t,1)]
    def soi(self,k,t,f=1.):
        return f*norm(self.planets[k](t))*(GM[k]/MU)**.4
    def acceleration(self,t,r,pert=True):
        a=grav(MU,r)
        if pert:
            for k,s in self.planets.items():
                R=s(t); a+=grav(GM[k],r-R)+grav(GM[k],R)
        return a
    def rhs(self,t,y,model='N',frame='H'):
        if frame=='H': a=self.acceleration(t,y[:3],model!='K')
        elif model=='P':
            s=self.planets[frame]
            a=self.acceleration(t,y[:3]+s(t))-s(t,2)
        else:
            a=grav(GM[frame],y[:3])
            if model=='P0':
                R=self.planets[frame](t)
                a+=grav(MU,R+y[:3])-grav(MU,R)
        return np.r_[y[3:],a]
    def run(self,y0,t0,tf,model='N',f=1.,frame='H',rtol=1e-12,
            atol=1e-8,max_step=DAY/2,method='DOP853'):
        if tf<=t0: raise ValueError('El tiempo final debe superar el inicial')
        if np.shape(y0)!=(6,) or not np.all(np.isfinite(y0)):
            raise ValueError('Se requiere un estado inicial finito de seis componentes')
        for k,spline in self.planets.items():
            if hasattr(spline,'x') and (t0<spline.x[0] or tf>spline.x[-1]):
                raise ValueError(f'Intervalo fuera de cobertura de efeméride {k}')
        segments=[]; jumps=[]; t=t0; y=np.array(y0,float)
        for _ in range(100000):
            events=[]; keys=list(self.planets)
            if model!='N':
                if frame=='H':
                    for k in keys:
                        def event(q,z,k=k): return norm(z[:3]-self.planets[k](q))-self.soi(k,q,f)
                        event.direction=-1; event.terminal=True; events.append(event)
                else:
                    def event(q,z): return norm(z[:3])-self.soi(frame,q,f)
                    event.direction=1; event.terminal=True; events=[event]
            end=tf
            if model=='P' and frame!='H':
                knots=self.planets[frame].x
                idx=np.searchsorted(knots,t+1e-5)
                if idx<len(knots): end=min(tf,float(knots[idx]))
            sol=solve_ivp(lambda q,z:self.rhs(q,z,model,frame),(t,end),y,
                          events=events or None,rtol=rtol,atol=atol,
                          max_step=max_step,method=method,dense_output=True)
            if not sol.success: raise RuntimeError(sol.message)
            segments.append((frame,sol)); old_t=t; t=sol.t[-1]; y=sol.y[:,-1]
            if t>=tf: break
            if t>=end-1e-6 and not any(len(a) for a in sol.t_events): continue
            if t<=old_t+1e-5: raise RuntimeError('Evento repetido sin avance temporal')
            before=y if frame=='H' else y+self.state(frame,t)
            if frame=='H':
                hits=[i for i,a in enumerate(sol.t_events) if len(a)]
                frame=keys[hits[0]]; y=y-self.state(frame,t)
            else:
                y=y+self.state(frame,t); frame='H'
            after=y if frame=='H' else y+self.state(frame,t)
            jumps.append(dict(t=float(t),frame_after=frame,dr_km=float(norm(after[:3]-before[:3])),
                              dv_kms=float(norm(after[3:]-before[3:]))))
        else: raise RuntimeError('Demasiados cambios de marco')
        return Trajectory(self,segments,jumps)

class Trajectory:
    def __init__(self,system,segments,jumps): self.system,self.segments,self.jumps=system,segments,jumps
    def __call__(self,t):
        if np.ndim(t): return np.array([self(q) for q in t])
        for fr,s in self.segments:
            if s.t[0]-1e-6<=t<=s.t[-1]+1e-6:
                y=s.sol(t)
                return y if fr=='H' else y+self.system.state(fr,t)
        raise ValueError(f'Tiempo fuera de trayectoria: {t}')

def invariants(y,mu):
    r,v=y[:3],y[3:]
    E=.5*v@v-mu/norm(r); h=np.cross(r,v)
    evec=np.cross(v,h)/mu-r/norm(r); e=norm(evec)
    return dict(energy_km2s2=float(E),h_km2s=float(norm(h)),e=float(e),
                delta_deg=float(np.degrees(2*np.arcsin(1/e))) if e>1 else None,
                vinf_kms=float(np.sqrt(2*E)) if E>0 else None)
