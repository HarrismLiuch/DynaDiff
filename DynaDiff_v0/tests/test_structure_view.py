"""Check camera scale and element colors through the structure widget."""

from pathlib import Path
import unittest

from PySide6.QtGui import QColor
from PySide6.QtWidgets import QApplication

from dynadiff.orientation import Orientation
from dynadiff.structure import load_structure
from dynadiff.views import StructureView


class StructureViewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication(["DynaDiff tests", "-platform", "offscreen"])
        cls.crystal = load_structure(Path(__file__).resolve().parents[1] / "examples/silicon.cif")

    def setUp(self):
        self.view = StructureView()
        self.addCleanup(self.view.close)
        self.view.resize(500, 400)
        self.view.set_crystal(self.crystal)
        self.front = Orientation.from_zone(self.crystal.cell, (0, 0, 1))
        self.tilted = Orientation.from_zone(self.crystal.cell, (1, 1, 1))
        self.view.set_frame(self.front)

    def test_rotation_and_roll_preserve_manual_scale(self):
        self.view.zoom = 1.7
        original = self.view.pixels_per_angstrom()
        self.view.set_frame(self.tilted.roll(29))
        self.assertEqual(self.view.pixels_per_angstrom(), original)
        self.view.fit_view()
        self.assertNotEqual(self.view.pixels_per_angstrom(), original)

    def test_optional_auto_fit_and_lock(self):
        self.view.set_auto_fit(True)
        original = self.view.pixels_per_angstrom()
        self.view.set_frame(self.tilted)
        fitted = self.view.pixels_per_angstrom()
        self.assertLess(fitted, original)
        self.view.set_auto_fit(False)
        self.view.set_frame(self.front)
        self.assertEqual(self.view.pixels_per_angstrom(), fitted)

    def test_presets_and_element_override(self):
        mineral = self.view.colors[14]
        self.view.set_atom_palette("Jmol")
        self.assertNotEqual(self.view.colors[14], mineral)
        self.view.color_overrides["Si"] = QColor("#a13f89")
        self.view.refresh_colors()
        self.assertEqual(self.view.colors[14].name(), "#a13f89")
        self.view.set_atom_palette("Mineral")
        self.assertEqual(self.view.colors[14], mineral)
