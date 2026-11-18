import unittest

import numpy as np

from src.constantes import DAY
from src.hohmann import R1, R2, calculate
from src.transferencia_leo import (
    EARTH_RADIUS, EARTH_SOI, NOMINAL_TIME, hohmann_burn, nominal_saturn_leg,
    nominal_scattering, transfer,
)


class ParkingOrbitTransfer(unittest.TestCase):
    def test_nominal_is_tangent_at_both_heliocentric_orbits(self):
        for altitude in (200, 500, 1500):
            with self.subTest(altitude=altitude):
                burn = hohmann_burn(altitude)
                result = transfer(altitude, burn)
                self.assertTrue(result.nominal)
                np.testing.assert_allclose(
                    np.linalg.norm(result.leo,axis=1),EARTH_RADIUS+altitude,
                )
                self.assertAlmostEqual(np.linalg.norm(result.escape[-1]), EARTH_SOI, delta=.01)
                np.testing.assert_allclose(result.solar[[0,-1]],[[R1,0],[-R2,0]],atol=.01)
                self.assertAlmostEqual(result.arrival_days, NOMINAL_TIME/DAY, places=7)
                self.assertAlmostEqual(result.phase_error_deg, 0., places=7)
        self.assertAlmostEqual(hohmann_burn(500), calculate()['dv_leo_kms'], places=9)

    def test_independent_burn_can_miss_or_fail_to_escape(self):
        burn = hohmann_burn(500)
        self.assertIsNone(transfer(500, 1).solar)
        self.assertIsNone(transfer(500, burn-.05).arrival_days)
        faster = transfer(500, burn+.05)
        self.assertFalse(faster.nominal)
        self.assertLess(faster.arrival_days, NOMINAL_TIME/DAY)
        self.assertGreater(abs(faster.phase_error_deg), 1.)

    def test_scattering_is_an_elastic_hyperbola_and_gains_solar_energy(self):
        from src.constantes import GM, RAD

        result = nominal_scattering()
        incoming, outgoing = result.incoming, result.outgoing
        np.testing.assert_allclose(np.linalg.norm(incoming),np.linalg.norm(outgoing))
        angle = np.arctan2(incoming[0]*outgoing[1]-incoming[1]*outgoing[0],incoming@outgoing)
        self.assertAlmostEqual(angle,result.turn_rad)
        self.assertAlmostEqual(np.min(np.linalg.norm(result.hyperbola,axis=1)),10*RAD['J'])
        v_inf = np.linalg.norm(incoming)
        self.assertAlmostEqual(result.eccentricity,1+result.periapsis_km*v_inf**2/GM['J'])
        solar_in = result.planet_velocity+incoming
        from_velocities = .5*(result.heliocentric_exit@result.heliocentric_exit-solar_in@solar_in)
        self.assertAlmostEqual(result.energy_change,from_velocities)
        self.assertGreater(result.energy_change,0)
        self.assertLess(result.heliocentric_exit[0],0)  # Radial hacia afuera en -x.

    def test_continuation_starts_at_hohmann_aphelion_and_reaches_saturn(self):
        leg = nominal_saturn_leg()
        np.testing.assert_allclose(leg.approach[0],[-R2,0],atol=.01)
        np.testing.assert_allclose(leg.approach[-1],leg.encounter_path[0],atol=.01)
        np.testing.assert_allclose(leg.departure[0],leg.encounter_path[-1],atol=.01)
        self.assertAlmostEqual(np.linalg.norm(leg.periapsis_position-leg.saturn_at_periapsis),500000,delta=1)
        self.assertAlmostEqual(leg.days_to_periapsis,1163.123,places=2)
        self.assertGreater(leg.days_after_periapsis,200)
        radii = np.linalg.norm(leg.relative_path,axis=1)
        self.assertEqual(np.argmin(radii),len(radii)//2)
        self.assertGreater(radii[-1],3*radii.min())
        incoming = leg.relative_path[1]-leg.relative_path[0]
        outgoing = leg.relative_path[-1]-leg.relative_path[-2]
        turn = np.arccos((incoming@outgoing)/(np.linalg.norm(incoming)*np.linalg.norm(outgoing)))
        self.assertGreater(np.degrees(turn),10)
