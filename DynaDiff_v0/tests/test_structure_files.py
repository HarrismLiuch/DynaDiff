"""File actions in the desktop workspace."""

from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from PySide6.QtCore import QSettings
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication, QFileDialog

from dynadiff.gui import EXAMPLE, MainWindow
from dynadiff.structure import load_structure


class StructureFileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication(["DynaDiff tests", "-platform", "offscreen"])

    def setUp(self):
        self.tempdir = TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)
        store = self.root / "preferences.ini"
        factory = patch("dynadiff.gui.QSettings", side_effect=lambda *_:
                        QSettings(str(store), QSettings.Format.IniFormat))
        factory.start()
        self.addCleanup(factory.stop)

    def window(self):
        window = MainWindow()
        self.addCleanup(window.close)
        return window

    def test_open_folder_survives_restart_example_cancel_and_failed_open(self):
        folder = self.root / "Measured crystals"
        folder.mkdir()
        path = folder / "Si.cif"
        copyfile(EXAMPLE, path)
        first = self.window()
        with patch.object(QFileDialog, "getOpenFileName", return_value=(str(path), "")):
            first.choose_structure()
        first.close()

        restarted = self.window()  # Opens bundled silicon automatically.
        with patch.object(restarted, "_error"):
            self.assertFalse(restarted.open_structure(self.root / "missing.cif"))
        with patch.object(QFileDialog, "getOpenFileName", return_value=("", "")) as picker:
            restarted.choose_structure()
            restarted.choose_structure()
            for call in picker.call_args_list:
                self.assertEqual(Path(call.args[2]), folder)

    def test_exports_available_while_diffraction_is_pending(self):
        window = self.window()
        window.show()
        window.set_zone([1, 1, 0])
        window.structure_view.set_repeat(2)
        self.assertIsNone(window.current_pattern)
        self.assertTrue(window.atomic_png_action.isEnabled())
        expected = window.structure_view.grab().toImage()
        dest = self.root / "Rotated atomic view.png"
        with patch.object(QFileDialog, "getSaveFileName", return_value=(str(dest), "")):
            window.atomic_png_action.trigger()
        image = QImage(str(dest))
        self.assertFalse(image.isNull())
        self.assertEqual(image, expected.convertToFormat(image.format()))

        cif = self.root / "Exported cell.cif"
        with patch.object(QFileDialog, "getSaveFileName", return_value=(str(cif), "")):
            window.atomic_cif_action.trigger()
        exported = load_structure(cif)
        self.assertEqual(len(exported.symbols), 8)
        self.assertEqual(exported.symbols, window.crystal.symbols)
