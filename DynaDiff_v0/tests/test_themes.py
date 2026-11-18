"""Application theme selection and its boundaries with the scientific views."""

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import QApplication

from dynadiff.gui import MainWindow
from dynadiff.themes import DEFAULT_THEME, THEMES


class ThemeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication(["DynaDiff tests", "-platform", "offscreen"])

    def setUp(self):
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        path = str(Path(temp.name) / "preferences.ini")
        self.store = QSettings(path, QSettings.Format.IniFormat)
        factory = patch("dynadiff.gui.QSettings", side_effect=lambda *_: QSettings(path, QSettings.Format.IniFormat))
        factory.start()
        self.addCleanup(factory.stop)
        self.addCleanup(self.app.setPalette, self.app.palette())
        self.addCleanup(self.app.setStyleSheet, self.app.styleSheet())

    def window(self):
        window = MainWindow()
        self.addCleanup(window.close)
        return window

    def test_theme_actions_change_widget_colors_and_preserve_scientific_views(self):
        window = self.window()
        window.show()
        window.set_zone([1, 1, 0])
        structure, pattern = window.structure_view, window.pattern_view
        orientation = window.orientation
        structure.set_atom_palette("Jmol")
        atom_colors = structure.colors.copy()
        window.detector_mode.setCurrentIndex(window.detector_mode.findData("detector"))
        settings = replace(pattern.detector_options, colormap="magma", reverse=True)
        pattern.set_detector_options(settings)
        with patch("dynadiff.detector_view.colormap_groups", return_value=[("Sequential", ("gray", "magma"))]):
            window.detector_settings_button.click()
            dialog = window.detector_dialog
            dialog._map_future.result(timeout=5)
            dialog._maps_loaded()
        pattern.zoom_in()
        zoom = pattern.zoom

        self.assertEqual(window.theme_id, DEFAULT_THEME)
        self.assertEqual(set(window.theme_actions), set(THEMES))
        appearance_menu = next(item.menu() for item in window.menuBar().actions()
                               if item.text().replace("&", "") == "Appearance")
        theme_menu = next(item.menu() for item in appearance_menu.actions() if item.text().replace("&", "") == "Theme")
        choices = list(window.theme_actions.values())
        self.assertEqual(theme_menu.actions(), choices)
        self.assertEqual(window.theme_button.menu().actions(), choices)
        palette_colors = set()
        heading_colors = set()
        for theme_id, action in window.theme_actions.items():
            with self.subTest(theme=theme_id):
                action.trigger()
                self.assertEqual(window.theme_id, theme_id)
                self.assertEqual([key for key, item in window.theme_actions.items() if item.isChecked()], [theme_id])
                palette = self.app.palette()
                palette_colors.add((palette.window().color().name(), palette.highlight().color().name()))
                ink = dialog.cmap.itemData(0, Qt.ItemDataRole.ForegroundRole).color().name()
                paper = dialog.cmap.itemData(0, Qt.ItemDataRole.BackgroundRole).color().name()
                self.assertNotEqual(ink, paper)
                heading_colors.add((ink, paper))
                self.assertIs(window.orientation, orientation)
                self.assertEqual(structure.colors, atom_colors)
                self.assertEqual(structure.grab().toImage().pixelColor(0, 0).name(), "#ffffff")
                self.assertEqual(pattern.detector_options, settings)
                self.assertEqual(pattern.zoom, zoom)
        self.assertEqual(len(palette_colors), len(THEMES))
        self.assertGreater(len(heading_colors), 1)

    def test_choice_survives_restart_and_unknown_preference_uses_default(self):
        self.store.setValue("appearance/theme", "removed-theme")
        self.store.sync()
        first = self.window()
        self.assertEqual(first.theme_id, DEFAULT_THEME)
        self.assertTrue(first.theme_actions[DEFAULT_THEME].isChecked())
        choice = next(key for key in THEMES if key != DEFAULT_THEME)
        first.theme_actions[choice].trigger()
        surface = self.app.palette().color(QPalette.ColorRole.Window)
        first.close()

        restarted = self.window()
        self.assertEqual(restarted.theme_id, choice)
        self.assertTrue(restarted.theme_actions[choice].isChecked())
        self.assertEqual(self.app.palette().color(QPalette.ColorRole.Window), surface)
