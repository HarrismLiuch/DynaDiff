"""Round-trip the expanded unit cell and its diffraction-relevant metadata."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from dynadiff.structure import Crystal, crystal_cif, load_structure


def round_trip(crystal):
    with TemporaryDirectory() as folder:
        path = Path(folder) / "export.cif"
        path.write_text(crystal_cif(crystal), encoding="utf-8")
        return load_structure(path)


class CrystalExportTests(unittest.TestCase):
    def test_expanded_silicon(self):
        example = Path(__file__).resolve().parents[1] / "examples" / "silicon.cif"
        silicon = load_structure(example)
        result = round_trip(silicon)
        self.assertEqual(result.spacegroup, "P 1")
        self.assertEqual(len(result.symbols), 8)
        self.assertEqual(result.symbols, silicon.symbols)
        np.testing.assert_allclose(result.frac, silicon.frac)
        np.testing.assert_allclose(result.cell, silicon.cell)

    def test_oblique_cell_and_chirality(self):
        geometry = np.array([[0.4, 4.1, 0.1], [-3.6, 0.7, 0.2], [0.8, -1.2, 3.2]])
        positions = [[0.12, 0.22, 0.04], [0.19, 0.33, 0.48], [0.71, 0.92, 0.66],
                     [0.44, 0.87, 0.31]]
        for handedness in (1, -1):
            with self.subTest(handedness=handedness):
                cell = geometry.copy()
                cell[2] *= handedness
                original = Crystal("A 'mixed' motif", cell, positions, ("O", "Cu", "Fe", "O"),
                                   [0.75, 0.5, 1, 0.25], [np.nan, 0, 1.24, 0.31], 150)
                result = round_trip(original)
                # Cartesian coordinates can undergo a rigid rotation when CIF
                # chooses its conventional frame. Reflections invert chirality.
                rotation = np.linalg.lstsq(original.cart, result.cart, rcond=None)[0]
                np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-10)
                self.assertAlmostEqual(np.linalg.det(rotation), 1.0, places=10)
                np.testing.assert_allclose(original.cart @ rotation, result.cart)
                self.assertEqual(result.name, original.name)
                self.assertEqual(result.reference_temperature_k, 150)
                self.assertEqual(result.symbols, original.symbols)
                np.testing.assert_allclose(result.occupancy, original.occupancy)
                np.testing.assert_allclose(result.b_iso, original.b_iso, equal_nan=True)
