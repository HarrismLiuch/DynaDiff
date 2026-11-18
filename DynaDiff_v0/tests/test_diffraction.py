"""Numerical regression checks for the shared crystal, view and diffraction pipeline."""

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from dynadiff.orientation import Orientation
from dynadiff.physics import SimulationSettings, simulate
from dynadiff.structure import Crystal, load_structure


EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "silicon.cif"


def frame(view):
    return np.stack([view.right, view.up, view.beam])


class DiffractionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.si = load_structure(EXAMPLE)
        cls.view = Orientation.from_zone(cls.si.cell, (0, 0, 1))
        cls.dynamic = SimulationSettings(mode="electron_dynamical")

    def test_zero_thickness_and_elastic_flux(self):
        for thickness in (0.0, 8.0, 20.0):
            with self.subTest(thickness=thickness):
                result = simulate(self.si, frame(self.view),
                                  replace(self.dynamic, thickness_nm=thickness))
                direct = np.all(result.hkl == 0, axis=1)
                self.assertEqual(np.count_nonzero(direct), 1)
                self.assertAlmostEqual(float(result.intensity.sum()), 1.0, places=9)
                if thickness == 0:
                    self.assertAlmostEqual(float(result.intensity[direct][0]), 1.0)
                else:
                    self.assertGreater(float(result.intensity[~direct].sum()), 1e-4)

    def test_silicon_extinctions_and_spacings(self):
        settings = replace(self.dynamic, mode="electron_kinematic")
        result = simulate(self.si, frame(self.view), settings)
        indices = {tuple(hkl): i for i, hkl in enumerate(result.hkl)}
        allowed = indices[(2, 2, 0)]
        self.assertGreater(result.intensity[allowed], 0)
        radius = np.linalg.norm(result.xy[allowed])
        self.assertAlmostEqual(radius, np.sqrt(8) / np.linalg.norm(self.si.cell[0]))
        for forbidden in [(1, 0, 0), (2, 0, 0)]:
            if forbidden in indices:
                self.assertLess(result.intensity[indices[forbidden]],
                                1e-12 * result.intensity[allowed])

    def test_two_beam_bragg_oscillation(self):
        # Only 000 and 100 enter this basis. At exact Bragg incidence the
        # analytical diffracted intensity is sin²(t λ γ f_100 / V).
        crystal = Crystal("Two-beam test", np.diag([4., 1., 1.]), [[0, 0, 0]],
                          ("Si",), [1.0], [0.0])
        wavelength, gamma = 0.0250793404, 1.39139024  # 200 kV
        g = 0.25
        x = -wavelength * g / 2
        view = Orientation.from_cartesian([x, 0, np.sqrt(1 - x*x)])
        f_100 = 3.9377072  # Si elastic electron factor in Å at sin(θ)/λ = 0.125 Å⁻¹
        rate = wavelength * gamma * f_100 / 4.0
        for phase in (np.pi / 4, np.pi / 2, np.pi):
            settings = replace(self.dynamic, thickness_nm=phase / (10 * rate),
                               max_beams=2, g_max=0.30)
            pattern = simulate(crystal, frame(view), settings)
            np.testing.assert_array_equal(pattern.hkl, [[0, 0, 0], [1, 0, 0]])
            self.assertAlmostEqual(pattern.intensity[1], np.sin(phase)**2, places=6)

    def test_roll_rotates_the_detector_with_the_preview(self):
        tilted = Orientation.from_zone(self.si.cell, (1, 1, 2))
        rolled = tilted.roll(34.0)
        np.testing.assert_allclose(rolled.beam, tilted.beam, atol=1e-12)
        for mode in ("xray_kinematic", "electron_kinematic", "electron_dynamical"):
            with self.subTest(mode=mode):
                result = simulate(self.si, frame(rolled), replace(self.dynamic, mode=mode))
                reciprocal = np.linalg.inv(self.si.cell).T
                expected = (result.hkl @ reciprocal) @ frame(rolled)[:2].T
                np.testing.assert_allclose(result.xy, expected, atol=1e-10)

    def test_nonorthogonal_poscar_direction(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "POSCAR"
            path.write_text("Hexagonal cell\n1.0\n3 0 0\n-1.5 2.598076 0\n0 0 5\n"
                            "Mg\n1\nDirect\n0 0 0\n", encoding="utf-8")
            crystal = load_structure(path)
        direction = np.array([2.0, 1.0, -1.0]) @ crystal.cell
        view = Orientation.from_zone(crystal.cell, (2, 1, -1))
        np.testing.assert_allclose(view.beam, direction / np.linalg.norm(direction))
        np.testing.assert_allclose(frame(view) @ frame(view).T, np.eye(3), atol=1e-12)


if __name__ == "__main__":
    unittest.main()
