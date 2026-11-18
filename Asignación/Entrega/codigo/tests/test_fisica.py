import tempfile
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import solve_ivp
from src.horizons import read,compatible
from src.dinamica import System,grav,invariants,norm
from src.constantes import MU,GM,DAY
from src.hohmann import calculate,R1,R2

HORIZONS=Path(__file__).resolve().parents[1]/'datos'/'horizons'
class Physics(unittest.TestCase):
    def test_hohmann_independent_integration(self):
        h=calculate();y=np.array([R1,0,0,0,h['vt1_kms'],0.])
        sol=solve_ivp(lambda t,z:np.r_[z[3:],grav(MU,z[:3])],(0,h['tof_days']*DAY),y,
                      rtol=1e-12,atol=1e-9,method='DOP853',dense_output=True)
        self.assertLess(norm(sol.y[:3,-1]-[-R2,0,0]),.1)
        vals=[invariants(sol.sol(t),MU) for t in np.linspace(0,sol.t[-1],100)]
        self.assertLess(np.ptp([a['energy_km2s2'] for a in vals])/abs(vals[0]['energy_km2s2']),1e-10)
        self.assertLess(np.ptp([a['h_km2s'] for a in vals])/vals[0]['h_km2s'],1e-10)

    def test_hyperbola_analytical(self):
        # Periapsis -> fronteras grandes; desviación finita converge a asíntota.
        rp=7e5;u=6.;mu=GM['J'];vp=np.sqrt(u*u+2*mu/rp)
        y=np.array([rp,0,0,0,vp,0.]);T=3000*DAY
        rhs=lambda t,z:np.r_[z[3:],grav(mu,z[:3])]
        sols=[solve_ivp(rhs,(0,t),y,rtol=1e-12,atol=1e-8,method='DOP853',dense_output=True) for t in (-T,T)]
        vi,vf=[s.y[3:,-1] for s in sols]
        angle=np.arccos(np.clip(vi@vf/(norm(vi)*norm(vf)),-1,1))
        delta=2*np.arcsin(1/(1+rp*u*u/mu))
        self.assertLess(abs(angle-delta),1e-4)
        for s in sols:
            inv=[invariants(s.sol(t),mu) for t in np.linspace(0,s.t[-1],200)]
            self.assertLess(np.ptp([a['energy_km2s2'] for a in inv])/(u*u/2),1e-9)
            self.assertLess(np.ptp([a['h_km2s'] for a in inv])/(rp*vp),1e-8)
        self.assertLess(2*np.arcsin(1/(1+1e16*u*u/mu)),1e-8)
        self.assertGreater(2*np.arcsin(1/(1+1e-6*u*u/mu)),3.14)

    def test_exact_frame_acceleration(self):
        e=read(HORIZONS/'jupiter.txt');ep=e.jd[0];s=e.spline(ep);sys=System({'J':s})
        for t in (600*DAY+123,650*DAY+567):
            rho=np.array([1e6,2e6,3e5]);z=np.r_[rho,[1,2,3]]
            actual=sys.rhs(t,z,'P','J')[3:]+s(t,2)
            self.assertTrue(np.allclose(actual,sys.acceleration(t,rho+s(t)),atol=1e-18,rtol=1e-12))

    def test_parser_and_derivative(self):
        e=read(HORIZONS/'Jupiter1979.txt');self.assertEqual(len(e.jd),457)
        self.assertAlmostEqual(e.states[0,0],-5.819361044472382e8,places=4)
        s=e.spline();q=4321.;eps=.1
        self.assertLess(norm((s(q+eps)-s(q-eps))/(2*eps)-s(q,1)),1e-5)
        self.assertTrue(np.all(np.isnan(s(-1))))

    def test_parser_labeled_and_duplicate(self):
        head='Target body name: Test\nCenter body name: Sun (10)\nOutput units: KM-S\nReference frame: Ecliptic of J2000.0\nOutput type: GEOMETRIC cartesian states\nStart time : A.D. 1979-Jul-01 TDB\n'
        blocks='\n'.join(f'{2451545+i}.0 = test\n X = 1D+08 Y = 2E+08 Z = 3E+08\n VX = 1 VY = 2 VZ = 3' for i in range(3))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.txt';p.write_text(head+'$$SOE\n'+blocks+'\n$$EOE')
            e=read(p);self.assertEqual(e.states.shape,(3,6))
            p.write_text(p.read_text().replace('2451546.0','2451545.0'))
            with self.assertRaisesRegex(ValueError,'duplicados'):read(p)
        bad=read(HORIZONS/'jupiter.txt');bad.meta['time_scale']='TT'
        with self.assertRaises(ValueError):compatible(read(HORIZONS/'saturno.txt'),bad)


class FrameIntegration(unittest.TestCase):
    def test_integrated_P_N(self):
        e=read(HORIZONS/'jupiter.txt');ep=e.jd[0];s=e.spline(ep);sys=System({'J':s})
        a=(2444063.5-ep)*DAY;b=a+2*DAY
        relative=np.array([1e6,2e6,3e5,-10.,-10.,1.])
        y=relative+sys.state('J',a)
        n=sys.run(y,a,b,rtol=2.3e-14,atol=1e-10,max_step=3600);p=sys.run(relative,a,b,'P',frame='J',rtol=2.3e-14,atol=1e-10,max_step=3600)
        self.assertLess(norm(p(b)[:3]-n(b)[:3]),.05)
        self.assertLess(norm(p(b)[3:]-n(b)[3:]),1e-6)

if __name__=='__main__': unittest.main()
