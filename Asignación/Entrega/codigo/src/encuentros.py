import datetime
import numpy as np
from scipy.optimize import minimize_scalar,brentq
from .constantes import DAY,GM
from .dinamica import invariants,norm

def jd(date):
    return 2451545.+(datetime.datetime.fromisoformat(date)-datetime.datetime(2000,1,1,12)).total_seconds()/DAY

def calendar(value):
    return (datetime.datetime(2000,1,1,12)+datetime.timedelta(days=value-2451545.)).isoformat(timespec='seconds')+' TDB'

def minimum(state,planet,t0,tf):
    fun=lambda t:norm(state(t)[:3]-planet(t))
    ts=np.linspace(t0,tf,max(101,int((tf-t0)/DAY)*4))
    ds=np.array([fun(t) for t in ts]); i=int(np.argmin(ds))
    if i in (0,len(ts)-1): return float(ts[i]),float(ds[i]),True
    s=minimize_scalar(fun,bounds=(ts[i-1],ts[i+1]),method='bounded',options={'xatol':.001})
    return float(s.x),float(s.fun),False

def encounter(state,system,k,t0,tf,epoch):
    t,r,bound=minimum(state,system.planets[k],t0,tf)
    relative=state(t)-system.state(k,t)
    inv=invariants(relative,GM[k])
    return dict(planet=k,jd_tdb=epoch+t/DAY,date_tdb=calendar(epoch+t/DAY),
                rp_km=r,vp_kms=float(norm(relative[3:])),boundary_minimum=bound,**inv)

def crossings(state,system,k,t0,tf,f=1.):
    fun=lambda t:norm(state(t)[:3]-system.planets[k](t))-system.soi(k,t,f)
    ts=np.arange(t0,tf,DAY/4); ts=np.r_[ts,tf]
    return [brentq(fun,a,b,xtol=.001) for a,b in zip(ts[:-1],ts[1:]) if fun(a)*fun(b)<0]
