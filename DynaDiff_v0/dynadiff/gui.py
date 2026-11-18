"""DynaDiff desktop workspace."""

import csv
from pathlib import Path
import re
import sys
from threading import Event

import numpy as np
from PySide6.QtCore import (
    QAbstractTableModel, QIODevice, QModelIndex, QObject, QPointF, QRunnable, QSaveFile,
    QSettings, Qt, QThreadPool, QTimer, Signal,
)
from PySide6.QtGui import (
    QAction, QActionGroup, QFont, QFontDatabase, QIcon, QKeySequence, QPainter, QPixmap,
)
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QColorDialog, QComboBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QFrame, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow, QMenu,
    QMessageBox, QPushButton, QScrollArea, QSpinBox, QSplitter, QTableView,
    QToolButton, QVBoxLayout, QWidget,
)

from .orientation import Orientation
from .physics import MAX_DYNAMICAL_BEAMS, SimulationSettings, simulate
from .structure import crystal_cif, load_structure
from .themes import DEFAULT_THEME, THEMES, apply_theme
from .views import PatternView, StructureView


EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "silicon.cif"
MODES = (
    ("Dynamical electron", "electron_dynamical"),
    ("Kinematic electron", "electron_kinematic"),
    ("Kinematic X-ray", "xray_kinematic"),
)


def label(text, role=None):
    widget = QLabel(text)
    if role:
        widget.setObjectName(role)
    return widget


def action(menu, text, function, shortcut=None):
    item = QAction(text, menu)
    item.triggered.connect(function)
    if shortcut:
        item.setShortcut(shortcut)
    menu.addAction(item)
    return item


def line(vertical=False):
    rule = QFrame()
    rule.setObjectName("rule")
    if vertical:
        rule.setFixedSize(1, 30)
    else:
        rule.setFixedHeight(1)
    return rule


def popup_button(text):
    button = QToolButton()
    button.setText(text + "  ▾")
    button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
    button.setMenu(QMenu(button))
    return button


class Logo(QWidget):
    def __init__(self):
        super().__init__()
        self.setFixedSize(35, 35)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        banner_palette = self.parentWidget().palette()
        p.setBrush(banner_palette.windowText())
        p.setPen(Qt.PenStyle.NoPen)
        p.drawRoundedRect(self.rect(), 8, 8)
        p.setBrush(banner_palette.window())
        for i, j, r in ((0, 0, 2.9), (-1, 0, 2.1), (1, 0, 2.1), (0, 1, 2.1),
                        (0, -1, 2.1), (1, 1, 1.25), (-1, -1, 1.25),
                        (1, -1, 1.25), (-1, 1, 1.25)):
            p.drawEllipse(QPointF(17.5 + i * 6.8 + j * 2.0, 17.5 + j * 7 - i * 1.6), r, r)


class ReflectionModel(QAbstractTableModel):
    columns = ("h", "k", "l", "|g|  (Å⁻¹)", "s  (Å⁻¹)", "I / I₀")

    def __init__(self):
        super().__init__()
        self.pattern = None
        self.order = np.array([], dtype=int)

    def replace(self, pattern, cell=None):
        self.beginResetModel()
        self.pattern = pattern
        if pattern is None:
            self.order = np.array([], dtype=int)
        else:
            self.order = np.argsort(-np.asarray(pattern.intensity), kind="stable")
            self.inverse = np.empty_like(self.order)
            self.inverse[self.order] = np.arange(len(self.order))
            g = pattern.hkl @ np.linalg.inv(cell).T
            self.g_length = np.linalg.norm(g, axis=1)
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.order)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.columns)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.DisplayRole:
            return self.columns[section]
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.ToolTipRole:
            return ("Miller index h", "Miller index k", "Miller index l",
                    "Reciprocal vector magnitude", "Excitation error",
                    "Linear intensity relative to the incident beam")[section]
        return None

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or self.pattern is None:
            return None
        if role == Qt.ItemDataRole.TextAlignmentRole:
            return Qt.AlignmentFlag.AlignCenter
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        i, column = self.order[index.row()], index.column()
        if column < 3:
            return str(self.pattern.hkl[i, column])
        if column == 3:
            return f"{self.g_length[i]:.5f}"
        if column == 4:
            return f"{self.pattern.excitation_error[i]:.5f}"
        return f"{self.pattern.intensity[i]:.6g}"


class WorkerSignals(QObject):
    finished = Signal(int, object, str)


class SimulationJob(QRunnable):
    def __init__(self, version, crystal, frame, settings, signals, cancel_event):
        super().__init__()
        self.version = version
        self.crystal = crystal
        self.frame = frame
        self.settings = settings
        self.signals = signals
        self.cancel_event = cancel_event

    def run(self):
        try:
            pattern = simulate(self.crystal, self.frame, self.settings, cancelled=self.cancel_event.is_set)
        except Exception as exc:
            self.signals.finished.emit(self.version, None, str(exc) or type(exc).__name__)
        else:
            self.signals.finished.emit(self.version, pattern, "")


class MainWindow(QMainWindow):
    def __init__(self, initial_path=None):
        super().__init__()
        self.setWindowTitle("DynaDiff")
        self.setMinimumSize(890, 620)
        self.setAcceptDrops(True)
        self.preferences = QSettings("DynaDiff", "DynaDiff")
        saved_theme = self.preferences.value("appearance/theme", DEFAULT_THEME, type=str)
        self.theme_id = apply_theme(QApplication.instance(), saved_theme)
        self.crystal = None
        self.source = None
        self.orientation = None
        self.current_pattern = None
        self._generation = 0
        self._busy = False
        self._cancel_event = Event()
        self._closed = False
        self._signals = WorkerSignals(self)
        self._signals.finished.connect(self._completed)
        self._pool = QThreadPool(self)
        self._pool.setMaxThreadCount(1)
        self._debounce = QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.timeout.connect(self._start_job)

        container = QWidget()
        container.setObjectName("workspace")
        self.setCentralWidget(container)
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._banner())
        content = QHBoxLayout()
        content.setSpacing(0)
        layout.addLayout(content, 1)
        content.addWidget(self._controls())
        content.addWidget(self._workspace(), 1)
        self._menus()
        self.statusBar().setSizeGripEnabled(False)
        self._enable_controls()
        area = QApplication.primaryScreen().availableGeometry()
        self.resize(min(area.width(), 1340), min(area.height(), 840))
        self.open_structure(initial_path if initial_path is not None else EXAMPLE,
                            remember_folder=initial_path is not None)

    def _banner(self):
        banner = QFrame()
        banner.setObjectName("banner")
        row = QHBoxLayout(banner)
        row.setContentsMargins(19, 12, 19, 12)
        row.setSpacing(12)
        row.addWidget(Logo())
        row.addWidget(label("DynaDiff", "brand"))
        row.addSpacing(8)
        row.addWidget(line(vertical=True))
        details = QVBoxLayout()
        details.setSpacing(0)
        self.file_label = label("Crystal structure", "file")
        self.cell_label = label("", "subtle")
        details.addWidget(self.file_label)
        details.addWidget(self.cell_label)
        row.addLayout(details)
        row.addStretch()
        self.theme_button = popup_button("Theme")
        self.theme_button.setObjectName("themeButton")
        self.theme_button.setAccessibleName("Interface theme")
        row.addWidget(self.theme_button)
        open_file = QPushButton("Open structure…")
        open_file.setObjectName("openButton")
        open_file.clicked.connect(self.choose_structure)
        row.addWidget(open_file)
        return banner

    def _menus(self):
        files = self.menuBar().addMenu("&File")
        action(files, "Open structure…", self.choose_structure, QKeySequence.StandardKey.Open)
        action(files, "Open example silicon", lambda: self.open_structure(EXAMPLE, remember_folder=False))
        files.addSeparator()
        files.addAction(self.atomic_png_action)
        files.addAction(self.atomic_cif_action)
        files.addSeparator()
        self.export_csv_action = action(files, "Export reflections as CSV…", self.export_csv)
        self.export_png_action = action(files, "Export pattern as PNG…", self.export_png)
        files.addSeparator()
        action(files, "Quit", self.close, QKeySequence.StandardKey.Quit)
        self._set_export_enabled(False)

        theme_menu = self.menuBar().addMenu("&Appearance").addMenu("&Theme")
        theme_group = QActionGroup(self)
        self.theme_actions = {}
        for theme_id, name in THEMES.items():
            item = theme_group.addAction(name)
            item.setData(theme_id)
            item.setCheckable(True)
            item.setChecked(theme_id == self.theme_id)
            theme_menu.addAction(item)
            self.theme_button.menu().addAction(item)
            self.theme_actions[theme_id] = item
        theme_group.triggered.connect(self._select_theme)
        self.theme_button.setToolTip(f"Interface theme · {THEMES[self.theme_id]}")

    def _select_theme(self, item):
        theme_id = item.data()
        if theme_id == self.theme_id:
            return
        self.theme_id = apply_theme(QApplication.instance(), theme_id)
        self.theme_button.setToolTip(f"Interface theme · {THEMES[theme_id]}")
        self.preferences.setValue("appearance/theme", theme_id)
        self.preferences.sync()

    def _controls(self):
        scroll = QScrollArea()
        scroll.setObjectName("controls")
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(265)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        side = QWidget()
        side.setObjectName("sidebar")
        scroll.setWidget(side)
        box = QVBoxLayout(side)
        box.setContentsMargins(17, 17, 17, 17)
        box.setSpacing(10)
        box.addWidget(label("Simulation", "section"))
        self.mode = QComboBox()
        self.mode.setObjectName("modeCombo")
        self.mode.setAccessibleName("Simulation mode")
        for name, value in MODES:
            self.mode.addItem(name, value)
        box.addWidget(self.mode)
        form = QFormLayout()
        form.setHorizontalSpacing(5)
        form.setVerticalSpacing(9)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        box.addLayout(form)
        self.voltage = self._number(form, "Voltage", "voltageSpin", 1, 2000, 200, 1, " kV", 10)
        self.wavelength = self._number(form, "Wavelength", "wavelengthSpin", 0.05, 10, 1.5406,
                                       4, " Å", 0.05)
        self.thickness = self._number(form, "Thickness", "thicknessSpin", 0, 1000, 20, 2, " nm", 1)
        self.temperature = self._number(form, "Temperature", "temperatureSpin", 0, 3000, 293.15,
                                        2, " K", 10)
        self.g_limit = self._number(form, "Max. |g|", "gMaxSpin", 0.01, 4, 1.25, 2, " Å⁻¹", 0.1)
        self.beams = QSpinBox()
        self.beams.setObjectName("beamSpin")
        self.beams.setRange(2, MAX_DYNAMICAL_BEAMS)
        self.beams.setValue(101)
        self.beams.setSingleStep(10)
        self.beams.setKeyboardTracking(False)
        self.beams.setToolTip("Maximum beams in the dynamical scattering calculation, including 000.")
        form.addRow("Beam limit", self.beams)
        self.form = form
        box.addSpacing(6)
        box.addWidget(line())
        box.addSpacing(2)
        box.addWidget(label("Scattering effects", "group"))
        self.dw = self._effect(box, "Debye–Waller attenuation", "debyeWallerCheck",
                               "Attenuate coherent scattering using thermal displacement factors.")
        self.multiple = self._effect(box, "Multiple scattering", "multipleScatteringCheck",
                                     "Allow amplitudes to scatter between diffracted beams.")
        self.relativistic = self._effect(box, "Relativistic electrons", "relativisticCheck",
                                        "Apply voltage-dependent relativistic corrections.")
        self.ewald = self._effect(box, "Ewald-sphere curvature", "ewaldCheck",
                                  "Include curvature when computing the excitation error.")
        box.addSpacing(7)
        calculate = QPushButton("Calculate pattern")
        calculate.setObjectName("simulateButton")
        calculate.setToolTip("Recalculate using the current settings and view direction.")
        calculate.clicked.connect(lambda: self._schedule(0))
        box.addWidget(calculate)
        self.auto_hint = label("Updates automatically when the view or settings change.", "subtle")
        self.auto_hint.setWordWrap(True)
        box.addWidget(self.auto_hint)
        box.addStretch()
        self.mode.currentIndexChanged.connect(self._settings_changed)
        for number in (self.voltage, self.wavelength, self.thickness, self.temperature,
                       self.g_limit, self.beams):
            number.valueChanged.connect(self._settings_changed)
        for effect in (self.dw, self.multiple, self.relativistic, self.ewald):
            effect.toggled.connect(self._settings_changed)
        return scroll

    @staticmethod
    def _number(form, title, name, minimum, maximum, value, digits, suffix, step):
        spin = QDoubleSpinBox()
        spin.setObjectName(name)
        spin.setRange(minimum, maximum)
        spin.setDecimals(digits)
        spin.setValue(value)
        spin.setSuffix(suffix)
        spin.setSingleStep(step)
        spin.setKeyboardTracking(False)
        form.addRow(title, spin)
        return spin

    @staticmethod
    def _effect(box, title, name, hint):
        item = QCheckBox(title)
        item.setObjectName(name)
        item.setChecked(True)
        item.setToolTip(hint)
        box.addWidget(item)
        return item

    def _workspace(self):
        area = QWidget()
        area.setObjectName("workspace")
        layout = QVBoxLayout(area)
        layout.setContentsMargins(13, 13, 13, 12)
        self.vertical_splitter = QSplitter(Qt.Orientation.Vertical)
        self.vertical_splitter.setHandleWidth(10)
        self.vertical_splitter.setChildrenCollapsible(False)
        layout.addWidget(self.vertical_splitter)
        self.views_splitter = QSplitter(Qt.Orientation.Horizontal)
        self.views_splitter.setChildrenCollapsible(False)
        self.views_splitter.setHandleWidth(10)
        self.views_splitter.addWidget(self._structure_card())
        self.views_splitter.addWidget(self._detector_card())
        self.views_splitter.setSizes([500, 500])
        self.vertical_splitter.addWidget(self.views_splitter)
        self.vertical_splitter.addWidget(self._reflections_card())
        self.vertical_splitter.setStretchFactor(0, 3)
        self.vertical_splitter.setStretchFactor(1, 1)
        self.vertical_splitter.setSizes([490, 212])
        return area

    @staticmethod
    def _card(title):
        card = QFrame()
        card.setObjectName("panel")
        column = QVBoxLayout(card)
        column.setContentsMargins(1, 1, 1, 1)
        column.setSpacing(0)
        header = QHBoxLayout()
        header.setContentsMargins(12, 9, 10, 9)
        header.setSpacing(8)
        header.addWidget(label(title, "section"))
        header.addStretch()
        column.addLayout(header)
        return card, column, header

    def _structure_card(self):
        card, layout, header = self._card("Atomic structure")
        self.structure_view = StructureView()
        menu_button = popup_button("View")
        header.addWidget(menu_button)
        menu = menu_button.menu()
        menu.addSection("Zone axis")
        for text, uvw in (("[001]", [0, 0, 1]), ("[100]", [1, 0, 0]), ("[110]", [1, 1, 0]),
                          ("[111]", [1, 1, 1])):
            action(menu, text, lambda checked=False, v=uvw: self.set_zone(v))
        action(menu, "Roll 90°", lambda: self._roll(90))
        menu.addSeparator()
        action(menu, "Fit current view", self.structure_view.fit_view)
        self.auto_fit_action = menu.addAction("Auto-fit on rotation")
        self.auto_fit_action.setCheckable(True)
        self.auto_fit_action.toggled.connect(self.structure_view.set_auto_fit)
        menu.addSection("Displayed cells")
        repeat_group = QActionGroup(menu)
        for n in (1, 2, 3):
            item = action(menu, f"{n} × {n} × {n}", lambda checked=False, v=n: self.structure_view.set_repeat(v))
            item.setCheckable(True)
            item.setChecked(n == 1)
            repeat_group.addAction(item)
        menu.addSeparator()
        bonds = menu.addAction("Show bonds")
        bonds.setCheckable(True)
        bonds.setChecked(True)
        bonds.toggled.connect(self._show_bonds)
        colors = menu.addMenu("Atom colors")
        palettes = QActionGroup(colors)
        for name in ("Mineral", "Jmol", "CPK"):
            item = action(colors, name, lambda checked=False, p=name: self.structure_view.set_atom_palette(p))
            item.setCheckable(True)
            item.setChecked(name == self.structure_view.atom_palette)
            palettes.addAction(item)
        colors.addSeparator()
        elements = colors.addMenu("Customize element…")
        elements.aboutToShow.connect(lambda: self._populate_element_colors(elements))
        menu.addSeparator()
        self.atomic_png_action = action(menu, "Export atomic view as PNG…", self.export_atomic_png)
        self.atomic_cif_action = action(menu, "Export unit cell as CIF…", self.export_atomic_cif)
        self.atomic_png_action.setEnabled(False)
        self.atomic_cif_action.setEnabled(False)
        self.structure_view.orbitRequested.connect(self._orbit)
        self.structure_view.rollRequested.connect(self._roll)
        layout.addWidget(self.structure_view, 1)
        footer = QVBoxLayout()
        footer.setContentsMargins(11, 8, 11, 9)
        footer.setSpacing(5)
        self.axis_readout = label("Beam direction", "readout")
        footer.addWidget(self.axis_readout)
        row = QHBoxLayout()
        row.setSpacing(5)
        self.direction_kind = QComboBox()
        self.direction_kind.setObjectName("directionKind")
        self.direction_kind.setAccessibleName("Direction coordinate system")
        self.direction_kind.addItems(["Zone [uvw]", "Cartesian"])
        self.direction_kind.currentIndexChanged.connect(self._sync_direction)
        self.direction_edit = QLineEdit("0 0 1")
        self.direction_edit.setObjectName("directionEdit")
        self.direction_edit.setAccessibleName("Beam direction components")
        self.direction_edit.setToolTip("Enter three direct-lattice indices or Cartesian components.")
        self.direction_edit.returnPressed.connect(self.apply_direction)
        self.direction_edit.setMinimumWidth(55)
        apply = QPushButton("Set")
        apply.setObjectName("directionApply")
        apply.clicked.connect(self.apply_direction)
        row.addWidget(self.direction_kind)
        row.addWidget(self.direction_edit, 1)
        row.addWidget(apply)
        footer.addLayout(row)
        layout.addLayout(footer)
        return card

    def _detector_card(self):
        card, layout, header = self._card("Diffraction")
        view_button = popup_button("Display")
        header.addWidget(view_button)
        menu = view_button.menu()
        self.grid_action = menu.addAction("Reciprocal-space grid")
        self.grid_action.setCheckable(True)
        self.grid_action.setChecked(True)
        self.label_action = menu.addAction("Label strongest reflections")
        self.label_action.setCheckable(True)
        self.pattern_view = PatternView()
        menu.addSeparator()
        self.addActions([
            action(menu, "Zoom in", self.pattern_view.zoom_in, QKeySequence.StandardKey.ZoomIn),
            action(menu, "Zoom out", self.pattern_view.zoom_out, QKeySequence.StandardKey.ZoomOut),
            action(menu, "Reset zoom", self.pattern_view.reset_zoom, QKeySequence("Ctrl+0")),
        ])
        self.pattern_view.zoomChanged.connect(self._update_detector_readout)
        self.grid_action.toggled.connect(self._display_settings)
        self.label_action.toggled.connect(self._display_settings)
        self.pattern_view.reflectionSelected.connect(self._select_reflection)
        layout.addWidget(self.pattern_view, 1)
        footer = QVBoxLayout()
        footer.setContentsMargins(12, 9, 12, 10)
        footer.setSpacing(5)
        self.detector_readout = label("λ    —", "readout")
        footer.addWidget(self.detector_readout)
        controls = QHBoxLayout()
        controls.setSpacing(5)
        self.detector_mode = QComboBox()
        self.detector_mode.setObjectName("detectorMode")
        self.detector_mode.setAccessibleName("Diffraction display mode")
        self.detector_mode.addItem("Indexed spots", "spots")
        self.detector_mode.addItem("Simulated detector", "detector")
        self.detector_mode.setToolTip("Indexed spot sizes use fourth-root intensity. The simulated detector broadens each reflection into a peak.")
        self.detector_mode.currentIndexChanged.connect(self._detector_mode_changed)
        controls.addWidget(self.detector_mode, 1)
        self.detector_settings_button = QPushButton("Settings…")
        self.detector_settings_button.setObjectName("detectorSettingsButton")
        self.detector_settings_button.setToolTip("Adjust the simulated detector peak shape, colormap and intensity range.")
        self.detector_settings_button.setEnabled(False)
        self.detector_settings_button.clicked.connect(self._open_detector_settings)
        controls.addWidget(self.detector_settings_button)
        self.detector_dialog = None
        footer.addLayout(controls)
        layout.addLayout(footer)
        return card

    def _detector_mode_changed(self):
        mode = self.detector_mode.currentData()
        self.pattern_view.set_mode(mode)
        self.grid_action.setChecked(self.pattern_view.show_grid)
        self.detector_settings_button.setEnabled(mode == "detector")
        if mode != "detector" and self.detector_dialog is not None:
            self.detector_dialog.hide()

    def _open_detector_settings(self):
        if self.detector_dialog is None:
            from .detector_view import DetectorControls

            self.detector_dialog = DetectorControls(self.pattern_view.detector_options, self)
            self.detector_dialog.changed.connect(self.pattern_view.set_detector_options)
        self.detector_dialog.show()
        self.detector_dialog.raise_()
        self.detector_dialog.activateWindow()

    def _reflections_card(self):
        card, layout, header = self._card("Reflections")
        self.count = label("", "count")
        header.insertWidget(1, self.count)
        self.notes_button = QPushButton("Model notes")
        self.notes_button.setObjectName("notesButton")
        self.notes_button.clicked.connect(self._notes)
        self.notes_button.setEnabled(False)
        header.addWidget(self.notes_button)
        self.export_button = popup_button("Export")
        self.export_button.setObjectName("exportButton")
        action(self.export_button.menu(), "Reflections as CSV…", self.export_csv)
        action(self.export_button.menu(), "Pattern as PNG…", self.export_png)
        header.addWidget(self.export_button)
        self.table = QTableView()
        self.table.setObjectName("reflectionTable")
        self.model = ReflectionModel()
        self.table.setModel(self.model)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().hide()
        self.table.verticalHeader().setDefaultSectionSize(24)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)
        self.table.selectionModel().currentRowChanged.connect(self._row_selected)
        self.table.setMinimumHeight(64)
        layout.addWidget(self.table, 1)
        return card

    def _enable_controls(self):
        mode = self.mode.currentData()
        dynamic = mode == "electron_dynamical"
        xray = mode == "xray_kinematic"
        self.form.setRowVisible(self.voltage, not xray)
        self.form.setRowVisible(self.wavelength, xray)
        self.form.setRowVisible(self.beams, dynamic)
        self.multiple.setVisible(dynamic)
        self.relativistic.setVisible(not xray)
        self.temperature.setEnabled(self.dw.isChecked())
        self.form.labelForField(self.temperature).setEnabled(self.dw.isChecked())

    def _settings_changed(self, *_):
        self._enable_controls()
        self._schedule()

    def settings(self):
        return SimulationSettings(
            mode=self.mode.currentData(), voltage_kv=self.voltage.value(),
            wavelength_angstrom=self.wavelength.value(), thickness_nm=self.thickness.value(),
            temperature_k=self.temperature.value(), g_max=self.g_limit.value(),
            max_beams=self.beams.value(), debye_waller=self.dw.isChecked(),
            multiple_scattering=self.multiple.isChecked(), relativistic=self.relativistic.isChecked(),
            ewald_curvature=self.ewald.isChecked(),
        )

    def choose_structure(self):
        folder = self.preferences.value("files/lastStructureFolder", "", type=str)
        if not folder or not Path(folder).is_dir():
            folder = str(self.source.parent if self.source and self.source.parent.is_dir() else Path.home())
        path, _ = QFileDialog.getOpenFileName(
            self, "Open crystal structure", folder,
            "Crystal structures (*.cif *.mcif *.mmcif *.pdb *.ent *.vasp *.poscar *.contcar "
            "POSCAR POSCAR.* CONTCAR CONTCAR.* *.xsf *.xyz *.extxyz *.pwi *.gz);;All files (*)")
        if path:
            self.open_structure(path)

    def open_structure(self, path, *, remember_folder=True):
        path = Path(path).expanduser().absolute()
        try:
            crystal = load_structure(path)
            frame = Orientation.from_zone(crystal.cell, (0, 0, 1))
        except (ValueError, OSError, RuntimeError) as exc:
            if self.crystal is None:
                self.pattern_view.set_message("Open a crystal structure to begin")
            self._error("Could not open structure", str(exc))
            return False
        self.source = path
        self.crystal = crystal
        self.orientation = frame
        self.structure_view.set_crystal(crystal)
        self.structure_view.set_frame(frame)
        self.atomic_png_action.setEnabled(True)
        self.atomic_cif_action.setEnabled(True)
        if remember_folder:
            self.preferences.setValue("files/lastStructureFolder", str(path.parent))
            self.preferences.sync()
        self.file_label.setText(path.name)
        self.file_label.setToolTip(str(path))
        elements = ", ".join(dict.fromkeys(crystal.symbols))
        n = len(crystal.frac)
        metadata = [f"{n:,} site{'s' if n != 1 else ''}"]
        if crystal.spacegroup:
            metadata.append(crystal.spacegroup)
        metadata.append(f"{abs(np.linalg.det(crystal.cell)):.2f} Å³")
        self.cell_label.setText("   ·   ".join(metadata))
        self.cell_label.setToolTip(f"{crystal.name}\nElements: {elements}\n"
                                  + "\n".join(f"{axis} = {np.linalg.norm(row):.4f} Å"
                                              for axis, row in zip("abc", crystal.cell)))
        self.setWindowTitle(f"{path.name} — DynaDiff")
        self.temperature.blockSignals(True)
        self.temperature.setValue(crystal.reference_temperature_k)
        self.temperature.blockSignals(False)
        self._sync_direction()
        self._schedule(0)
        return True

    def _sync_direction(self, *_):
        if self.orientation is None:
            return
        zone = self.orientation.label(self.crystal.cell)
        self.axis_readout.setText(f"Beam {zone}    ·    drag to rotate")
        if self.direction_kind.currentIndex() == 0 and not zone.startswith("≈"):
            text = zone.strip("[]")
        else:
            vec = (self.orientation.zone(self.crystal.cell) if self.direction_kind.currentIndex() == 0
                   else self.orientation.beam)
            text = " ".join(f"{v:.6g}" for v in vec)
        self.direction_edit.setText(text)
        self._direction_error(False)

    def _direction_error(self, state):
        self.direction_edit.setProperty("invalid", state)
        self.direction_edit.style().unpolish(self.direction_edit)
        self.direction_edit.style().polish(self.direction_edit)

    def apply_direction(self):
        if self.crystal is None:
            return
        raw = self.direction_edit.text().strip().strip("[]()").replace("−", "-")
        parts = raw.replace(",", " ").split()
        if self.direction_kind.currentIndex() == 0 and len(parts) == 1 and re.fullmatch(r"([+-]?\d){3}", raw):
            parts = re.findall(r"[+-]?\d", raw)
        try:
            values = [float(s) for s in parts]
            if self.direction_kind.currentIndex() == 0:
                frame = Orientation.from_zone(self.crystal.cell, values)
            else:
                frame = Orientation.from_cartesian(values)
        except (ValueError, OverflowError):
            self._direction_error(True)
            self.statusBar().showMessage("Enter three finite components with at least one nonzero value.")
            return
        self._change_frame(frame)

    def set_zone(self, uvw):
        if self.crystal is not None:
            self._change_frame(Orientation.from_zone(self.crystal.cell, uvw))

    def _orbit(self, yaw, pitch):
        if self.orientation is not None:
            self._change_frame(self.orientation.orbit(-yaw, pitch))

    def _roll(self, angle):
        if self.orientation is not None:
            self._change_frame(self.orientation.roll(angle))

    def _change_frame(self, frame):
        self.orientation = frame
        self.structure_view.set_frame(frame)
        self._sync_direction()
        self._schedule()

    def _show_bonds(self, checked):
        self.structure_view.show_bonds = checked
        self.structure_view.update()

    def _populate_element_colors(self, menu):
        menu.clear()
        if self.crystal is None:
            return
        elements = sorted(set(zip(self.crystal.numbers, self.crystal.symbols)))
        for number, symbol in elements:
            swatch = QPixmap(16, 16)
            swatch.fill(self.structure_view.colors[number])
            item = menu.addAction(QIcon(swatch), symbol)
            item.triggered.connect(lambda checked=False, z=number, s=symbol: self._choose_element_color(z, s))

    def _choose_element_color(self, number, symbol):
        view = self.structure_view
        color = QColorDialog.getColor(view.colors[number], self, f"{symbol} atom color")
        if color.isValid():
            view.color_overrides[symbol] = color
            view.refresh_colors()

    def _display_settings(self, *_):
        self.pattern_view.show_grid = self.grid_action.isChecked()
        self.pattern_view.show_labels = self.label_action.isChecked()
        self.pattern_view.update()

    def _schedule(self, delay=180):
        if self.crystal is None or self._closed:
            return
        self._generation += 1
        self._cancel_event.set()
        self.current_pattern = None
        self.pattern_view.g_max = self.g_limit.value()
        self.pattern_view.set_message("Calculating…")
        if self.model.pattern is not None:
            self.model.replace(None)
        self.count.clear()
        self.detector_readout.setText("λ    —")
        self._set_export_enabled(False)
        self.notes_button.setEnabled(False)
        self.statusBar().showMessage("Calculating diffraction pattern…")
        self._debounce.start(delay)

    def _start_job(self):
        if self._busy or self._closed:
            return
        self._busy = True
        self._cancel_event = Event()
        self._pool.start(SimulationJob(self._generation, self.crystal, self.orientation,
                                       self.settings(), self._signals, self._cancel_event))

    def _completed(self, version, pattern, error):
        self._busy = False
        if self._closed:
            return
        if version != self._generation:
            if not self._debounce.isActive():
                self._start_job()
            return
        if error:
            self.statusBar().showMessage(error)
            self.pattern_view.set_message("Calculation unavailable")
            self.pattern_view.setToolTip(error)
            return
        self.current_pattern = pattern
        self.pattern_view.setToolTip("")
        self.pattern_view.set_pattern(pattern, self.g_limit.value())
        self.model.replace(pattern, self.crystal.cell)
        self.count.setText(f"{len(pattern.hkl):,} beams")
        self._update_detector_readout()
        self._set_export_enabled(True)
        self.notes_button.setEnabled(True)
        self.statusBar().showMessage(f"{self.mode.currentText()}    ·    {self.thickness.value():g} nm    ·    "
                                     f"Beam {self.orientation.label(self.crystal.cell)}")

    def _update_detector_readout(self, *_):
        if self.current_pattern is not None:
            self.detector_readout.setText(
                f"λ {self.current_pattern.wavelength_angstrom:.5f} Å  ·  "
                f"|g| ≤ {self.g_limit.value():g} Å⁻¹  ·  {self.pattern_view.zoom:.0%}")

    def _set_export_enabled(self, enabled):
        self.export_button.setEnabled(enabled)
        self.export_csv_action.setEnabled(enabled)
        self.export_png_action.setEnabled(enabled)

    def _row_selected(self, current, previous):
        if current.isValid() and current.row() < len(self.model.order):
            self.pattern_view.set_selected(int(self.model.order[current.row()]))

    def _select_reflection(self, index):
        if self.model.pattern is not None:
            row = int(self.model.inverse[index])
            self.table.selectRow(row)
            self.table.scrollTo(self.model.index(row, 0))

    def _notes(self):
        if self.current_pattern is None:
            return
        mode = self.mode.currentData()
        if mode == "electron_dynamical" and self.multiple.isChecked():
            summary = ("Elastic Bloch-wave calculation for a perfect, parallel-sided crystal. "
                       "The table gives each beam’s fraction of incident intensity.")
        else:
            summary = ("Single-scattering slab calculation. Intensities are relative to the "
                       "undepleted incident beam. Strong reflections can exceed unity.")
        notes = [summary]
        if self.dw.isChecked():
            notes.append("Isotropic B factors scale linearly with temperature.")
        notes.extend(self.current_pattern.notes)
        self._dialog("Model notes", "\n\n".join(notes), QMessageBox.Icon.Information)

    def _dialog(self, title, text, icon):
        dialog = QMessageBox(icon, title, text, QMessageBox.StandardButton.Ok, self)
        dialog.setTextFormat(Qt.TextFormat.PlainText)
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        dialog.open()

    def _error(self, title, message):
        self.statusBar().showMessage(message)
        self._dialog(title, message, QMessageBox.Icon.Warning)

    def export_csv(self):
        if self.current_pattern is None:
            return
        destination, _ = QFileDialog.getSaveFileName(
            self, "Export reflection intensities", str(self.source.with_suffix(".reflections.csv")),
            "CSV (*.csv)")
        if not destination:
            return
        if not Path(destination).suffix:
            destination += ".csv"
        # The modal file picker processes events, so recheck the live result.
        pattern = self.current_pattern
        if pattern is None:
            return
        try:
            with open(destination, "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["h", "k", "l", "qx_A^-1", "qy_A^-1", "g_A^-1",
                                 "excitation_error_A^-1", "I_over_I0"])
                for i in self.model.order:
                    writer.writerow([*pattern.hkl[i], *pattern.xy[i], self.model.g_length[i],
                                     pattern.excitation_error[i], pattern.intensity[i]])
        except OSError as exc:
            self._error("Could not export CSV", str(exc))
            return
        self.statusBar().showMessage(f"Saved {Path(destination).name}")

    def export_png(self):
        if self.current_pattern is not None:
            self._export_canvas_png(self.pattern_view, ".diffraction.png", "Export diffraction pattern")

    def export_atomic_png(self):
        if self.crystal is not None:
            self._export_canvas_png(self.structure_view, ".structure.png", "Export atomic view")

    def _export_canvas_png(self, canvas, suffix, title):
        # Save the view that the user saw when they requested the export.
        snapshot = canvas.grab()
        destination, _ = QFileDialog.getSaveFileName(
            self, title, str(self.source.with_suffix(suffix)), "PNG image (*.png)")
        if destination:
            self._write_export(destination, ".png", lambda device: snapshot.save(device, "PNG"))

    def export_atomic_cif(self):
        crystal = self.crystal
        if crystal is None:
            return
        destination, _ = QFileDialog.getSaveFileName(
            self, "Export unit cell", str(self.source.with_suffix(".unit-cell.cif")),
            "Crystallographic information file (*.cif)")
        if destination:
            def write_cif(device):
                payload = crystal_cif(crystal).encode("utf-8")
                return device.write(payload) == len(payload)
            self._write_export(destination, ".cif", write_cif)

    def _write_export(self, destination, suffix, write):
        if not Path(destination).suffix:
            destination += suffix
        output = QSaveFile(destination)
        try:
            if not output.open(QIODevice.OpenModeFlag.WriteOnly):
                raise OSError(output.errorString())
            if not write(output) or not output.commit():
                raise OSError(output.errorString())
        except (OSError, ValueError, RuntimeError) as exc:
            if output.isOpen():
                output.cancelWriting()
            self._error("Could not export file", str(exc))
            return
        self.statusBar().showMessage(f"Saved {Path(destination).name}")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls() and any(u.isLocalFile() for u in event.mimeData().urls()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            if url.isLocalFile():
                self.open_structure(url.toLocalFile())
                event.acceptProposedAction()
                return

    def closeEvent(self, event):
        self._closed = True
        self._cancel_event.set()
        self._debounce.stop()
        self._pool.clear()
        self._pool.waitForDone()
        super().closeEvent(event)


def main(argv=None):
    arguments = sys.argv if argv is None else argv
    app = QApplication(arguments)
    app.setApplicationName("DynaDiff")
    app.setStyle("Fusion")
    typeface = app.font()
    installed = set(QFontDatabase.families())
    for family in ("Inter", "Helvetica Neue", "Arial", "DejaVu Sans"):
        if family in installed:
            typeface.setFamily(family)
            break
    typeface.setStyleHint(QFont.StyleHint.SansSerif)
    app.setFont(typeface)
    paths = [arg for arg in arguments[1:] if not arg.startswith("-psn_")]
    window = MainWindow(paths[0] if paths else None)
    window.show()
    return app.exec()
