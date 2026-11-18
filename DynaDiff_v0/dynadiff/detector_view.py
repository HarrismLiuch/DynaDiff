"""Display controls and cached Qt rasters for the simulated diffraction detector."""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import QEvent, QSize, QSortFilterProxyModel, Qt, QTimer, Signal
from PySide6.QtGui import (
    QFontMetrics, QGuiApplication, QIcon, QImage, QPalette, QPixmap, QStandardItem, QTransform, qRgba,
)
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QCompleter, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
    QLabel, QVBoxLayout,
)

from .detector import PROFILES, color_lut, colormap_groups, detector_image, scale_image


@dataclass(frozen=True)
class DetectorOptions:
    profile: str = "Gaussian"
    fwhm: float = 0.025
    colormap: str = "gray"
    reverse: bool = False
    clip_percent: float = 0.1
    logarithmic: bool = True


PREVIEW_SIZE = QSize(68, 14)
_SEARCH_TEXT = Qt.ItemDataRole.UserRole + 1


class _ColormapHeading(QStandardItem):
    def data(self, role=Qt.ItemDataRole.UserRole + 1):
        # Both the main picker and the search popup share this item, including
        # when a user changes the theme while detector settings are open.
        if role == Qt.ItemDataRole.ForegroundRole:
            return QGuiApplication.palette().brush(QPalette.ColorGroup.Active, QPalette.ColorRole.Text)
        if role == Qt.ItemDataRole.BackgroundRole:
            return QGuiApplication.palette().brush(QPalette.ColorGroup.Active, QPalette.ColorRole.AlternateBase)
        return super().data(role)


class _ColormapCompleter(QCompleter):
    """Keep matching palettes under their published headings in source order."""

    def __init__(self, picker):
        proxy = QSortFilterProxyModel(picker)
        proxy.setSourceModel(picker.model())
        proxy.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        proxy.setFilterRole(_SEARCH_TEXT)
        super().__init__(proxy, picker)

    def splitPath(self, text):
        self.model().setFilterFixedString(text.strip())
        # The proxy includes matching group headers. A second name filter
        # in QCompleter would hide them along with their category context.
        return [""]

    def eventFilter(self, widget, event):
        if event.type() == QEvent.Type.KeyPress and widget is self.popup():
            key = event.key()
            if key in (Qt.Key.Key_Down, Qt.Key.Key_Up):
                # QCompleter selects the first row even when it is a disabled
                # heading. Navigate the enabled rows of this grouped popup.
                model = widget.model()
                step = 1 if key == Qt.Key.Key_Down else -1
                row = widget.currentIndex().row()
                if row < 0 and step < 0:
                    row = model.rowCount()
                for i in range(model.rowCount()):
                    candidate = model.index((row + (i + 1) * step) % model.rowCount(), 0)
                    if candidate.flags() & Qt.ItemFlag.ItemIsSelectable:
                        widget.setCurrentIndex(candidate)
                        break
                return True
        return super().eventFilter(widget, event)


def _colormap_icons(name):
    """Use the detector's exact RGBA table, including steps in discrete maps."""
    lut = color_lut(name)
    strip = QImage(lut.data, len(lut), 1, lut.nbytes, QImage.Format.Format_RGBA8888)
    pixmap = QPixmap.fromImage(strip.scaled(PREVIEW_SIZE * 2))
    pixmap.setDevicePixelRatio(2)
    return QIcon(pixmap), QIcon(pixmap.transformed(QTransform().scale(-1, 1)))


class DetectorControls(QDialog):
    changed = Signal(object)

    def __init__(self, options, parent):
        super().__init__(parent)
        self.setWindowTitle("Detector settings")
        self.setModal(False)
        self.setMinimumWidth(420)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 16, 18, 14)
        layout.setSpacing(12)
        form = QFormLayout()
        form.setVerticalSpacing(11)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        layout.addLayout(form)

        self.profile = QComboBox()
        self.profile.setObjectName("detectorProfile")
        self.profile.addItems(PROFILES)
        self.profile.setCurrentText(options.profile)
        self.profile.setToolTip("Gaussian models instrumental blur. Airy models an ideal circular aperture or domain.")
        form.addRow("Peak profile", self.profile)
        self.width_spin = QDoubleSpinBox()
        self.width_spin.setObjectName("detectorFwhm")
        self.width_spin.setDecimals(3)
        self.width_spin.setRange(0.001, 0.250)
        self.width_spin.setSingleStep(0.005)
        self.width_spin.setSuffix(" Å⁻¹")
        self.width_spin.setValue(options.fwhm)
        self.width_spin.setKeyboardTracking(False)
        self.width_spin.setToolTip("Full width at half maximum across a peak in reciprocal space.")
        form.addRow("Peak FWHM", self.width_spin)

        self.cmap = QComboBox()
        self.cmap.setObjectName("detectorColormap")
        self.cmap.setAccessibleName("Detector colormap")
        self.cmap.setEditable(True)
        self.cmap.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.cmap.setIconSize(PREVIEW_SIZE)
        self._active_colormap = options.colormap
        self.cmap.addItem(options.colormap, options.colormap)
        self.cmap.setEnabled(False)
        self.cmap.setMinimumContentsLength(12)
        self.cmap.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.cmap.setCompleter(None)
        self.cmap.lineEdit().setPlaceholderText("Search colormaps…")
        self.cmap.setToolTip("Loading Matplotlib and scientific colormaps…")
        form.addRow("Colormap", self.cmap)
        self.reverse = QCheckBox("Reverse colormap")
        self.reverse.setObjectName("detectorReverse")
        self.reverse.setChecked(options.reverse)
        self._preview_icons = {}
        form.addRow("", self.reverse)
        self.clip = QDoubleSpinBox()
        self.clip.setObjectName("detectorOutliers")
        self.clip.setDecimals(2)
        self.clip.setRange(0.0, 99.99)
        self.clip.setSingleStep(0.1)
        self.clip.setSuffix(" %")
        self.clip.setValue(options.clip_percent)
        self.clip.setKeyboardTracking(False)
        self.clip.setToolTip("Saturate this percentage of the brightest detector pixels. Set 0 for the full range.")
        form.addRow("Clip brightest pixels", self.clip)
        self.log_scale = QCheckBox("Log intensity")
        self.log_scale.setObjectName("detectorLog")
        self.log_scale.setChecked(options.logarithmic)
        self.log_scale.setToolTip("Use a logarithmic color scale spanning four decades below the upper limit.")
        form.addRow("", self.log_scale)
        note = QLabel("Peak width is shared by all reflections. Each beam contributes its integrated intensity.")
        note.setObjectName("subtle")
        note.setWordWrap(True)
        layout.addWidget(note)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)

        for combo in (self.profile, self.cmap):
            combo.currentIndexChanged.connect(self._changed)
        self.cmap.editTextChanged.connect(self._match_colormap)
        self.cmap.lineEdit().editingFinished.connect(self._finish_colormap)
        for switch in (self.reverse, self.log_scale):
            switch.toggled.connect(self._changed)
        self.reverse.toggled.connect(self._refresh_previews)
        for spin in (self.width_spin, self.clip):
            spin.valueChanged.connect(self._changed)
        # On first use the palette packages can build Matplotlib's font cache.
        # Keep the instrument controls and both viewports responsive meanwhile.
        pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="colormaps")
        self._map_future = pool.submit(colormap_groups)
        pool.shutdown(wait=False)
        self._map_timer = QTimer(self)
        self._map_timer.timeout.connect(self._maps_loaded)
        self._map_timer.start(75)

    def _maps_loaded(self):
        if not self._map_future.done():
            return
        self._map_timer.stop()
        try:
            groups = self._map_future.result()
        except Exception as error:
            from matplotlib import colormaps
            from .mpl_groups import matplotlib_groups

            groups = [(f"Matplotlib · {category}", names)
                      for category, names in matplotlib_groups(colormaps)]
            self.cmap.setToolTip(f"Scientific colormaps could not load. {error}")
        else:
            self.cmap.setToolTip("Search by colormap name, collection or category.")
        self.cmap.blockSignals(True)
        self.cmap.clear()
        model = self.cmap.model()
        heading_font = self.cmap.font()
        heading_font.setBold(True)
        for label, names in groups:
            heading = _ColormapHeading(label)
            heading.setFlags(Qt.ItemFlag.NoItemFlags)
            heading.setFont(heading_font)
            heading.setData("\x1f".join((label, *names)), _SEARCH_TEXT)
            model.appendRow(heading)
            for name in names:
                self.cmap.addItem(name, name)
                self.cmap.setItemData(self.cmap.count() - 1, f"{label}\x1f{name}", _SEARCH_TEXT)
        self.cmap.setCurrentIndex(self.cmap.findData(self._active_colormap))
        self._refresh_previews()
        self.cmap.setCompleter(_ColormapCompleter(self.cmap))
        search_popup = self.cmap.completer().popup()
        search_popup.setIconSize(PREVIEW_SIZE)
        font_metrics = QFontMetrics(heading_font)
        heading_width = max(font_metrics.horizontalAdvance(label) for label, _ in groups)
        for view in (self.cmap.view(), search_popup):
            view.setMinimumWidth(heading_width + view.verticalScrollBar().sizeHint().width() + 12)
        self.cmap.blockSignals(False)
        self.cmap.setEnabled(True)

    def _refresh_previews(self):
        reversed_map = self.reverse.isChecked()
        for i in range(self.cmap.count()):
            name = self.cmap.itemData(i)
            if name is None:
                continue
            if name not in self._preview_icons:
                self._preview_icons[name] = _colormap_icons(name)
            self.cmap.setItemIcon(i, self._preview_icons[name][reversed_map])

    def _match_colormap(self, text):
        index = self.cmap.findText(text, Qt.MatchFlag.MatchFixedString)
        if index >= 0 and self.cmap.itemData(index) and index != self.cmap.currentIndex():
            self.cmap.setCurrentIndex(index)

    def _finish_colormap(self):
        self.cmap.setEditText(self.cmap.itemText(self.cmap.currentIndex()))

    def _changed(self, *_):
        name = self.cmap.currentData()
        if name is None:
            # An editable combo can try to commit an exactly typed heading
            # on Return even though that row is disabled in the popup.
            self.cmap.setCurrentIndex(self.cmap.findData(self._active_colormap))
            return
        self._active_colormap = name
        self.changed.emit(DetectorOptions(
            profile=self.profile.currentText(), fwhm=self.width_spin.value(),
            colormap=name, reverse=self.reverse.isChecked(),
            clip_percent=self.clip.value(), logarithmic=self.log_scale.isChecked(),
        ))


class DetectorRaster:
    """Keep the detector exposure, display levels and color table separately."""

    def __init__(self):
        self.exposure_key = None
        self.level_key = None
        self.color_key = None
        self.image = QImage()
        self.colorbar = QImage()

    def prepare(self, pattern, g_max, options):
        exposure_key = (id(pattern), g_max, options.profile, options.fwhm)
        if exposure_key != self.exposure_key:
            self.exposure = detector_image(pattern.xy, pattern.intensity, g_max,
                                           profile=options.profile, fwhm=options.fwhm)
            self.exposure_key = exposure_key
            self.level_key = None
        level_key = (options.clip_percent, options.logarithmic)
        levels_changed = level_key != self.level_key
        if levels_changed:
            levels, self.lower, self.upper = scale_image(
                self.exposure, clip_percent=options.clip_percent, logarithmic=options.logarithmic)
            indices = np.ascontiguousarray(np.rint(levels * 255), dtype=np.uint8)
            height, width = indices.shape
            self.image = QImage(indices.data, width, height, indices.strides[0],
                                QImage.Format.Format_Indexed8).copy()
            self.level_key = level_key
        color_key = (options.colormap, options.reverse)
        color_changed = color_key != self.color_key
        if color_changed:
            colors = color_lut(*color_key)
            self.table = [qRgba(*map(int, color)) for color in colors]
            indices = np.arange(256, dtype=np.uint8)
            self.colorbar = QImage(indices.data, len(indices), 1, len(indices),
                                   QImage.Format.Format_Indexed8).copy()
            self.colorbar.setColorTable(self.table)
            self.color_key = color_key
        if levels_changed or color_changed:
            self.image.setColorTable(self.table)
