"""The detector picker exposes the same swatches in its full and filtered lists."""

import unittest

from PySide6.QtCore import Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from dynadiff.detector import color_lut, colormap_groups
from dynadiff.detector_view import DetectorControls, DetectorOptions, PREVIEW_SIZE


class ColormapPickerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication(["DynaDiff tests", "-platform", "offscreen"])

    def setUp(self):
        self.dialog = DetectorControls(DetectorOptions(), None)
        self.addCleanup(self.dialog.close)
        self.dialog._map_future.result(timeout=60)
        self.dialog._maps_loaded()

    def test_group_headings_order_and_keyboard_selection(self):
        picker = self.dialog.cmap
        groups = colormap_groups()
        expected_rows = [entry for label, names in groups for entry in (label, *names)]
        self.assertEqual([picker.itemText(i) for i in range(picker.count())], expected_rows)
        for label, _ in groups:
            row = picker.findText(label)
            self.assertFalse(picker.model().index(row, 0).flags() & Qt.ItemFlag.ItemIsSelectable)
            self.assertIsNone(picker.itemData(row))

        self.dialog.show()
        picker.setCurrentIndex(picker.findData(groups[0][1][-1]))
        QTest.keyClick(picker, Qt.Key.Key_Down)
        self.assertEqual(picker.currentData(), groups[1][1][0])
        picker.lineEdit().selectAll()
        QTest.keyClicks(picker, "batlow")
        popup = picker.completer().popup()
        matches = []
        for label, names in groups:
            filtered = [n for n in names if "batlow" in n.casefold()]
            if filtered:
                matches.extend([label, *filtered])
        model = popup.model()
        self.assertEqual([model.index(i, 0).data() for i in range(model.rowCount())], matches)
        QTest.keyClick(popup, Qt.Key.Key_Down)
        self.assertEqual(popup.currentIndex().data(), matches[1])
        QTest.keyClick(popup, Qt.Key.Key_Return)
        self.assertEqual(picker.currentData(), matches[1])

    def test_swatches_match_detector_and_reverse_in_search(self):
        dialog = self.dialog
        picker = dialog.cmap
        self.assertTrue(all(not picker.itemIcon(i).isNull() for i in range(picker.count())
                            if picker.itemData(i)))
        completer = picker.completer()
        for reversed_map in (False, True):
            dialog.reverse.setChecked(reversed_map)
            for name in ("gray", "cmc.batlow", "cmr.eclipse"):
                picker.setCurrentIndex(picker.findText(name))
                completer.setCompletionPrefix(name)
                model = completer.completionModel()
                row = next(i for i in range(model.rowCount()) if model.index(i, 0).data() == name)
                current = picker.itemIcon(picker.currentIndex())
                filtered = model.index(row, 0).data(Qt.ItemDataRole.DecorationRole)
                self.assertEqual(current.cacheKey(), filtered.cacheKey())
                image = current.pixmap(PREVIEW_SIZE).toImage()
                table = color_lut(name, reversed_map)
                for x, i in ((0, 0), (image.width() - 1, -1)):
                    rgba = image.pixelColor(x, image.height() // 2).getRgb()
                    self.assertTrue(all(abs(int(a) - int(b)) <= 5 for a, b in zip(rgba, table[i])))
