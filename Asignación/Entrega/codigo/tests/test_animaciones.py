"""Comprobaciones de trazabilidad del muestreo, independientes de OpenGL."""
import json
import unittest
import numpy as np
from src.animacion_datos import get_data
from src.constantes import DAY

class AnimationSampling(unittest.TestCase):
    def test_historical_prefix_and_frames(self):
        d=get_data('esferas');m=json.loads(str(d['metadata']))
        self.assertTrue(np.all(np.diff(d['t'])>0))
        entry=(m['entry_jd']-m['epoch_jd'])*DAY
        # La conversión JD introduce decenas de microsegundos: excluir el punto frontera.
        shared=d['t']<entry-.1
        for name in ('N','P','P0','K','Voyager'):
            self.assertEqual(d[name].shape,(len(d['t']),3))
            self.assertTrue(np.all(np.isfinite(d[name])))
        for name in ('P','P0','K'):
            np.testing.assert_array_equal(d[name][shared],d['N'][shared])
            self.assertTrue(set(np.unique(d['region_'+name])).issubset({0,1,2}))
        # K f=1 falla Saturno; la animación no inventa una entrada.
        self.assertNotIn(2,d['region_K'])
        self.assertIn(2,d['region_P'])
        for key in ('J','S'):
            self.assertLess(np.min(abs(d['t']-m['peri_'+key])),1e-6)

    def test_hohmann_keeps_its_initial_conditions(self):
        from src.hohmann import R1,R2
        d=get_data('hohmann');m=json.loads(str(d['metadata']))
        np.testing.assert_allclose(d['Hohmann'][0],[R1,0,0],atol=1e-7)
        index=np.argmin(abs(d['t']-m['peri_J']))
        np.testing.assert_allclose(d['Hohmann'][index],[-R2,0,0],atol=.1)
        self.assertTrue(np.all(np.isfinite(d['Hohmann'])))
        self.assertNotEqual(m['historical_end_jd']-m['historical_start_jd'],d['t'][-1]/DAY)
        self.assertLess(abs(m['design']['refined_saturn_rp_km']-500000),1.)

class AnimationTelemetry(unittest.TestCase):
    def test_initial_hohmann_speed_and_energy(self):
        from src.hohmann import calculate,R1,R2
        from src.constantes import MU
        d=get_data('hohmann');h=calculate()
        np.testing.assert_allclose(d['v_Hohmann'][0],[0,h['vt1_kms'],0],atol=1e-10)
        energy=.5*np.dot(d['v_Hohmann'][0],d['v_Hohmann'][0])-MU/R1
        self.assertAlmostEqual(energy,-MU/(R1+R2),places=8)
        for key in ('v_Hohmann','hist_v_Voyager'):
            self.assertTrue(np.all(np.isfinite(d[key])))

    def test_voyager_velocity_uses_ephemeris_derivative(self):
        from src.animacion_datos import load_ephemerides,voyager_positions
        epoch,spl=load_ephemerides();d=get_data('esferas')
        times=(d['jd']-epoch)*DAY
        np.testing.assert_allclose(d['v_Voyager'],voyager_positions(spl,times,1),atol=1e-5)
        for key in ('N','P','P0','K'):
            self.assertEqual(d['v_'+key].shape,d[key].shape)
            self.assertTrue(np.all(np.isfinite(d['v_'+key])))

class JupiterEncounterTelemetry(unittest.TestCase):
    def test_hyperbolic_periapsis_and_scattering(self):
        from scipy.integrate import solve_ivp
        from src.constantes import GM
        from src.dinamica import System,grav
        from src.encuentros import encounter
        mu=GM['J'];rp=714920.;vinf=6.
        vp=np.sqrt(vinf**2+2*mu/rp)
        initial=np.array([rp,0,0,0,vp,0])
        rhs=lambda t,y:np.r_[y[3:],grav(mu,y[:3])]
        past=solve_ivp(rhs,(0,-86400),initial,dense_output=True,rtol=1e-11,atol=1e-8)
        future=solve_ivp(rhs,(0,86400),initial,dense_output=True,rtol=1e-11,atol=1e-8)
        class Planet:
            def __call__(self,t,nu=0):return np.zeros(3)
        result=encounter(lambda t:(past if t<0 else future).sol(t),System({'J':Planet()}),'J',-86400,86400,2451545.)
        self.assertFalse(result['boundary_minimum'])
        self.assertAlmostEqual(result['rp_km'],rp,places=3)
        expected=np.degrees(2*np.arcsin(1/(1+rp*vinf**2/mu)))
        self.assertAlmostEqual(result['delta_deg'],expected,places=7)

    def test_hohmann_uses_design_periapsis(self):
        d=get_data('hohmann');m=json.loads(str(d['metadata']))
        result=m['jupiter_encounters']['Hohmann']
        self.assertEqual(result['rp_km'],m['design']['rp_jupiter_km'])
        self.assertEqual(result['delta_deg'],m['design']['delta_deg'])
        voyager=m['jupiter_encounters']['Voyager']
        self.assertFalse(voyager['boundary_minimum'])
        self.assertGreater(voyager['rp_km'],0)
        self.assertTrue(0<voyager['delta_deg']<180)


class SuperposedComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from src.animacion_comparacion import comparison
        cls.data=comparison();cls.meta=json.loads(str(cls.data['metadata']))

    def test_esferas_starts_at_earth_without_shared_prefix(self):
        d=self.data
        np.testing.assert_array_equal(d['K'][0],d['T'][0])
        np.testing.assert_array_equal(d['Hohmann'][0],d['T'][0])
        np.testing.assert_allclose(d['v_K'][0],self.meta['initial_state_K'][3:],atol=1e-12)
        self.assertNotIn('N',d);self.assertNotIn('P',d);self.assertNotIn('P0',d)
        self.assertEqual(set(self.meta['jupiter_encounters']),{'Hohmann','K','Voyager'})

    def test_voyager_remains_in_original_solar_frame(self):
        from src.animacion_datos import load_ephemerides,voyager_positions
        epoch,spl=load_ephemerides();d=self.data
        times=(self.meta['epoch_jd']-epoch)*DAY+d['t']
        np.testing.assert_allclose(d['Voyager'],voyager_positions(spl,times),atol=1e-6)
        np.testing.assert_allclose(d['v_Voyager'],voyager_positions(spl,times,1),atol=1e-10)

    def test_common_clock_and_scattering_validity(self):
        d=self.data;self.assertTrue(np.all(np.diff(d['t'])>0))
        self.assertEqual(d['t'][0],0)
        for key in ('Hohmann','K','Voyager','T','J','S'):
            self.assertEqual(d[key].shape,(len(d['t']),3))
            self.assertTrue(np.all(np.isfinite(d[key])))
        if not self.meta['K_enters_J']:
            self.assertIsNone(self.meta['jupiter_encounters']['K']['delta_deg'])
        self.assertTrue(np.all(np.isfinite(d['v_K'])))

    def test_earth_sphere_departures_match_trajectory_samples(self):
        d=self.data
        for key in ('Hohmann','K'):
            departure=self.meta['earth_departures'][key]
            self.assertEqual(departure['status'],'crossing')
            crossing=departure['elapsed_seconds']
            self.assertGreater(crossing,0)
            self.assertLess(crossing,3*DAY)
            index=int(np.argmin(abs(d['t']-crossing)))
            distance=np.linalg.norm(d[key]-d['T'],axis=1)
            self.assertAlmostEqual(d['t'][index],crossing,places=6)
            self.assertAlmostEqual(distance[index],d['soi_T'][index],delta=.02)
            self.assertLess(distance[0],d['soi_T'][0])
        voyager=self.meta['earth_departures']['Voyager']
        self.assertEqual(voyager['status'],'outside_at_start')
        self.assertIsNone(voyager['elapsed_seconds'])
        self.assertGreater(voyager['initial_distance_km'],d['soi_T'][0])

    def test_continuous_jupiter_scattering_and_osculating_metrics(self):
        from src.constantes import GM, MU, RAD
        d=self.data;m=self.meta
        j=m['jupiter_encounters']['Hohmann']
        self.assertAlmostEqual(j['rp_km'],10*RAD['J'],delta=1.)
        hyper=m['targeting']['Hohmann_hyperbola']
        self.assertAlmostEqual(hyper['rp_km'],j['rp_km'],delta=1.)
        self.assertAlmostEqual(j['delta_deg'],np.degrees(2*np.arcsin(1/hyper['e'])),places=6)
        # The velocity derivative must match gravity through periapsis.
        # A hidden impulse or rounded/interpolated corner violates this bound.
        index=np.argmin(np.linalg.norm(d['Hohmann']-d['J'],axis=1))
        indices=np.arange(max(0,index-20),min(len(d['t'])-1,index+20))
        for i in indices:
            dt=d['t'][i+1]-d['t'][i]
            if dt<.1:continue
            r=(d['Hohmann'][i]+d['Hohmann'][i+1])/2
            acceleration=-MU*r/np.linalg.norm(r)**3
            for key in ('J','S'):
                R=(d[key][i]+d[key][i+1])/2
                rho=r-R
                acceleration-=GM[key]*(rho/np.linalg.norm(rho)**3+R/np.linalg.norm(R)**3)
            measured=(d['v_Hohmann'][i+1]-d['v_Hohmann'][i])/dt
            np.testing.assert_allclose(measured,acceleration,rtol=.01,atol=1e-8)
        # Gravitational assist increases heliocentric energy across Jupiter.
        before=index;after=index
        while before>0 and d['t'][index]-d['t'][before]<20*DAY:before-=1
        while after<len(d['t'])-1 and d['t'][after]-d['t'][index]<20*DAY:after+=1
        energy=lambda i: .5*np.dot(d['v_Hohmann'][i],d['v_Hohmann'][i])-MU/np.linalg.norm(d['Hohmann'][i])
        self.assertGreater(energy(after),energy(before))

    def test_both_models_encounter_saturn_and_continue(self):
        from src.constantes import RAD
        d=self.data;m=self.meta
        self.assertGreaterEqual(d['t'][-1]/DAY,3652.5)
        for key in ('Hohmann','K'):
            saturn=m['saturn_encounters'][key]
            self.assertFalse(saturn['boundary_minimum'])
            self.assertTrue(saturn['inside_soi'])
            self.assertGreater(saturn['rp_km'],RAD['S'])
            self.assertLess(saturn['rp_km'],2e6)
            self.assertGreater(d['t'][-1]/DAY-saturn['elapsed_days'],2*365.25)
            self.assertLess(m['targeting'][key+'_target_error_km'],1e4)
            self.assertGreater(m['jupiter_encounters'][key]['rp_km'],RAD['J'])
            self.assertGreater(np.linalg.norm(d[key][-1]-d[key][-100]),1e6)
        self.assertTrue(any(event['frame_after']=='S' for event in m['K_events']))

class LambertTargeting(unittest.TestCase):
    def test_transfer_reaches_endpoint_under_solar_gravity(self):
        from scipy.integrate import solve_ivp
        from src.ajuste_saturno import lambert
        from src.constantes import AU,MU
        from src.dinamica import grav
        r1=np.array([AU,0.,0.]);r2=5.2*AU*np.array([np.cos(2.),np.sin(2.),.01])
        duration=700*DAY;v0,v1=lambert(r1,r2,duration)
        integration=solve_ivp(lambda t,y:np.r_[y[3:],grav(MU,y[:3])],(0,duration),np.r_[r1,v0],rtol=1e-12,atol=1e-8)
        np.testing.assert_allclose(integration.y[:3,-1],r2,atol=.02,rtol=0)
        np.testing.assert_allclose(integration.y[3:,-1],v1,atol=1e-8,rtol=0)

if __name__=='__main__':unittest.main()
