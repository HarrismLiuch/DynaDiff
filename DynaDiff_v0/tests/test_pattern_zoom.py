"""Exercise diffraction navigation through the rendered Qt canvas."""

import unittest
from unittest.mock import patch

import numpy as np
from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import QWheelEvent
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from dynadiff.physics import Pattern
from dynadiff.views import PatternView, scale_bar


class PatternZoomTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication(["DynaDiff tests", "-platform", "offscreen"])

    def setUp(self):
        self.view = PatternView()
        self.addCleanup(self.view.close)
        self.view.resize(480, 360)
        self.pattern = Pattern(np.array([[0, 0, 0], [1, 0, 0]]),
                               np.array([[0., 0.], [.375, .25]]),
                               np.array([1., .5]), np.zeros(2), .025, ())
        self.view.set_pattern(self.pattern, 1.0)
        self.view.show()

    def wheel(self, pos, pixels=0, angle=0):
        event = QWheelEvent(pos, pos, QPoint(0, pixels), QPoint(0, angle),
                            Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
                            Qt.ScrollPhase.ScrollUpdate, False)
        self.app.sendEvent(self.view, event)

    def test_zoom_anchors_spot_and_double_click_resets_in_both_modes(self):
        for mode in ("spots", "detector"):
            with self.subTest(mode=mode):
                self.view.set_mode(mode)
                self.view.grab()
                original = self.view._screen.copy()
                anchor = QPointF(*original[1])
                self.wheel(anchor, angle=120)
                self.view.grab()
                self.assertAlmostEqual(self.view.zoom, 1.25)
                np.testing.assert_allclose(self.view._screen[1], original[1])
                np.testing.assert_allclose(self.view._screen[1] - self.view._screen[0],
                                           1.25 * (original[1] - original[0]))
                self.assertEqual(self.view.nearest(anchor), 1)
                self.wheel(anchor, pixels=20)
                self.assertGreater(self.view.zoom, 1.25)
                QTest.mouseDClick(self.view, Qt.MouseButton.LeftButton, pos=anchor.toPoint())
                self.view.grab()
                np.testing.assert_allclose(self.view._screen, original)
                self.assertEqual(self.view.zoom, 1.0)

    def test_detector_zoom_preserves_colorbar_exposure_and_updates_reciprocal_scale(self):
        self.view.set_mode("detector")
        scales = []

        def draw_ruler(painter, width, height, scale, unit, color):
            scales.append(scale)
            scale_bar(painter, width, height, scale, unit, color)

        with patch("dynadiff.views.scale_bar", draw_ruler):
            before = self.view.grab().toImage()
            old_scale = scales[-1]
            exposure = self.view.raster.exposure
            self.view.zoom_in()
            after = self.view.grab().toImage()
            self.assertAlmostEqual(scales[-1], 1.25 * old_scale)
        self.assertIs(self.view.raster.exposure, exposure)
        height = round(45 * before.devicePixelRatio())
        self.assertEqual(before.copy(0, before.height()-height, before.width(), height),
                         after.copy(0, after.height()-height, after.width(), height))
        self.assertNotEqual(before, after)
        self.assertEqual(self.view.nearest(QPointF(240, 340)), -1)
        self.wheel(QPointF(240, 340), angle=120)
        self.assertEqual(self.view.zoom, 1.25)
        self.view.set_message("Calculating…")
        self.view.set_pattern(self.pattern, 1.0)
        self.assertEqual(self.view.zoom, 1.25)
        self.view.zoom_out()
        self.assertEqual(self.view.zoom, 1.)
