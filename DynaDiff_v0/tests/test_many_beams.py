"""The large-basis propagator must preserve the dense elastic Bloch model."""

from concurrent.futures import CancelledError
from dataclasses import replace
import unittest
from unittest.mock import patch

import numpy as np
from scipy import fft

from dynadiff.orientation import Orientation
from dynadiff.physics import SimulationSettings, _compact_lattice, simulate
from dynadiff.structure import Crystal


class ManyBeamTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # An oblique lattice and a noncentrosymmetric, partially occupied motif
        # exercise nontrivial complex structure factors and skew-zone projections.
        cls.crystal = Crystal("Oblique crystal", [[7, 0, 0], [1, 6, 0], [1, 1.5, 6.5]],
                              [[.12, .23, .34], [.28, .61, .78]], ("Al", "O"),
                              [1.0, .71], [.3, .5])
        cls.frame = Orientation.from_cartesian([.83, .41, .22])
        cls.opts = SimulationSettings(mode="electron_dynamical", max_beams=419,
                                      thickness_nm=7)

    def test_matrix_free_converges_to_dense_exit_intensities(self):
        actual = simulate(self.crystal, self.frame, self.opts)
        self.assertGreater(len(actual.hkl), 401)
        with patch("dynadiff.physics._DENSE_BEAM_LIMIT", 500):
            reference = simulate(self.crystal, self.frame, self.opts)
        np.testing.assert_array_equal(actual.hkl, reference.hkl)
        np.testing.assert_allclose(actual.intensity, reference.intensity,
                                   atol=2e-12, rtol=1e-9)

    def test_integer_compaction_of_sparse_zone(self):
        h, k = np.meshgrid(np.arange(-18, 19), np.arange(-14, 15))
        hkl = np.column_stack((2*h.ravel(), k.ravel(), (177*h - 93*k).ravel()))
        compact, transform = _compact_lattice(hkl)
        np.testing.assert_array_equal(compact @ transform, hkl)
        self.assertLess(np.prod(2*np.ptp(compact, axis=0) + 1), 10_000)

    def test_maximum_requested_basis_and_zero_thickness(self):
        crystal = Crystal("Beam limit", np.diag([24., 22., 6.]),
                          [[.21, .16, .1]], ("Si",), [1.0], [.5])
        frame = Orientation.from_cartesian([.057, .081, 1])
        for thickness in (0., .5):
            with self.subTest(thickness=thickness):
                pattern = simulate(crystal, frame, replace(self.opts, max_beams=9999,
                                                           thickness_nm=thickness))
                self.assertEqual(len(pattern.hkl), 9999)
                self.assertAlmostEqual(float(pattern.intensity.sum()), 1., places=10)
                if thickness == 0:
                    self.assertEqual(float(pattern.intensity[0]), 1.)
                else:
                    self.assertGreater(float(pattern.intensity[1:].sum()), 0.)
        with self.assertRaisesRegex(ValueError, "9999"):
            simulate(crystal, frame, replace(self.opts, max_beams=10_000))

    def test_cancellation_while_propagating_and_workspace_limit(self):
        with patch("dynadiff.physics.fft.ifftn", wraps=fft.ifftn) as inverse:
            with self.assertRaises(CancelledError):
                simulate(self.crystal, self.frame, self.opts, cancelled=lambda: inverse.call_count >= 2)
        with patch("dynadiff.physics._MAX_FFT_POINTS", 10):
            with self.assertRaisesRegex(ValueError, "workspace limit"):
                simulate(self.crystal, self.frame, self.opts)
