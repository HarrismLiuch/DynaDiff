"""Qt canvases for the direct-space structure and the reciprocal-space detector."""

from itertools import product
from math import floor, log10

import numpy as np
from ase.data import chemical_symbols, covalent_radii
from ase.data.colors import cpk_colors, jmol_colors
from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QColor, QFont, QPainter, QPen, QRadialGradient,
)
from PySide6.QtWidgets import QApplication, QSizePolicy, QToolTip, QWidget
from scipy.spatial import cKDTree


def font(size, bold=False):
    result = QApplication.font()
    result.setPixelSize(size)
    if bold:
        result.setWeight(QFont.Weight.DemiBold)
    return result


def pen(color, width=1):
    return QPen(QColor(color), width)


def point(xy):
    return QPointF(float(xy[0]), float(xy[1]))


def scale_length(maximum):
    decade = 10 ** floor(log10(max(maximum, 1e-12)))
    return max(factor * decade for factor in (0.5, 1, 2, 5, 10)
               if factor * decade <= maximum)


def scale_bar(painter, width, height, px_per_unit, unit, color):
    length = scale_length(min(85, width * 0.22) / px_per_unit)
    x, y = width - 18, height - 17
    left = x - length * px_per_unit
    painter.setPen(pen(color, 1.4))
    painter.drawLine(QPointF(left, y), QPointF(x, y))
    for end in (left, x):
        painter.drawLine(QPointF(end, y - 3), QPointF(end, y + 2))
    painter.setFont(font(11))
    painter.drawText(QRectF(x - 98, y - 21, 98, 17),
                     Qt.AlignmentFlag.AlignRight, f"{length:g} {unit}")


# Mineral pigments for common elements. The fallback hue is stable by atomic
# number so opening another structure never reassigns its element colors.
MINERAL_COLORS = {
    "H": "#eee4d1", "Li": "#ae80b4", "B": "#bf8062", "C": "#51495a",
    "N": "#416fb0", "O": "#be5746", "F": "#83aa80", "Na": "#9674b5",
    "Mg": "#83a76b", "Al": "#9aa9ba", "Si": "#4684a1", "P": "#d38e3e",
    "S": "#d4af4e", "Cl": "#559480", "K": "#8061a5", "Ca": "#a6a751",
    "Ti": "#7b91a8", "Cr": "#829863", "Mn": "#b4738b", "Fe": "#a45540",
    "Co": "#45638c", "Ni": "#749489", "Cu": "#bd7848", "Zn": "#839ca3",
    "Ga": "#a482a1", "Ge": "#8293ae", "As": "#99799b", "Se": "#b68d47",
    "Br": "#965543", "Ag": "#a9abb1", "I": "#77569a", "W": "#667995",
    "Pt": "#a5a4a8", "Au": "#c5a14b", "Pb": "#73718b",
}


class StructureView(QWidget):
    orbitRequested = Signal(float, float)
    rollRequested = Signal(float)

    def __init__(self):
        super().__init__()
        self.setObjectName("structureView")
        self.setAccessibleName("Atomic structure preview")
        self.setToolTip("Drag to rotate. Shift-drag to roll. Scroll to zoom. Double-click to fit.")
        self.setMinimumSize(220, 180)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        self.crystal = None
        self.frame = None
        self.repeat = 1
        self.show_bonds = True
        self.zoom = 1.0
        self.auto_fit = False
        self._fixed_spans = None
        self.atom_palette = "Mineral"
        self.color_overrides = {}
        self._drag = None

    def set_crystal(self, crystal):
        self.crystal = crystal
        self.frame = None
        self._fixed_spans = None
        self.zoom = 1.0
        self.prepare_geometry()

    def set_repeat(self, repeat):
        self.repeat = repeat
        if self.crystal is not None:
            self.prepare_geometry()
            self.fit_view()

    def set_frame(self, frame):
        self.frame = frame
        if self._fixed_spans is None:
            self._fixed_spans = self.projected_spans()
        self.update()

    def projected_spans(self):
        screen_axes = np.stack((self.frame.right, self.frame.up))
        return np.maximum(np.ptp((self.corners - self.center) @ screen_axes.T, axis=0), 1e-6)

    def pixels_per_angstrom(self):
        spans = self.projected_spans() if self.auto_fit else self._fixed_spans
        return min((self.width() - 78) / spans[0], (self.height() - 78) / spans[1]) * self.zoom

    def fit_view(self):
        if self.frame is not None:
            self._fixed_spans = self.projected_spans()
            self.zoom = 1.0
            self.update()

    def set_auto_fit(self, enabled):
        if enabled:
            self.fit_view()
        elif self.frame is not None:
            # Freeze the current fit, including any manual zoom offset.
            self._fixed_spans = self.projected_spans()
        self.auto_fit = enabled
        self.update()

    def set_atom_palette(self, name):
        self.atom_palette = name
        self.color_overrides.clear()
        self.refresh_colors()

    def refresh_colors(self):
        if self.crystal is None:
            return
        table = {"Jmol": jmol_colors, "CPK": cpk_colors}.get(self.atom_palette)
        self.colors = {}
        for z in np.unique(self.crystal.numbers):
            symbol = chemical_symbols[z]
            if symbol in self.color_overrides:
                color = self.color_overrides[symbol]
            elif table is not None and z < len(table):
                color = QColor.fromRgbF(*table[z])
            elif symbol in MINERAL_COLORS:
                color = QColor(MINERAL_COLORS[symbol])
            else:
                color = QColor.fromHsv(int(z * 137.508) % 360, 149, 185)
            self.colors[int(z)] = color
        self.update()

    def prepare_geometry(self):
        crystal = self.crystal
        frac = np.mod(crystal.frac, 1.0)
        frac[np.isclose(frac, 1, atol=1e-7, rtol=0)] = 0
        frac[np.abs(frac) < 1e-7] = 0
        # Include boundary images so the edges of a conventional cell are complete.
        # Thin very large cells deterministically for interactive drawing.
        max_draw = 4500
        step = max(1, int(np.ceil(len(frac) * self.repeat ** 3 / max_draw)))
        self.sampled = step > 1
        frac = frac[::step]
        numbers = np.asarray(crystal.numbers)[::step]
        occupancy = np.asarray(crystal.occupancy)[::step]
        positions, species, occupancies = [], [], []
        for offset in product(range(self.repeat + 1), repeat=3):
            shifted = frac + np.array(offset)
            keep = np.all(shifted <= self.repeat + 1e-8, axis=1)
            if np.any(keep):
                positions.append(shifted[keep])
                species.append(numbers[keep])
                occupancies.append(occupancy[keep])
        self.corners = np.array(list(product((0, self.repeat), repeat=3))) @ crystal.cell
        self.center = np.mean(self.corners, axis=0)
        self.positions = np.concatenate(positions) @ crystal.cell - self.center
        self.species = np.concatenate(species)
        self.occupancies = np.concatenate(occupancies)
        self.refresh_colors()
        self.radii = covalent_radii[self.species]
        self.bonds = np.empty((0, 2), dtype=int)
        if len(self.positions) < 1600 and not self.sampled:
            cutoff = 2.32 * np.max(self.radii)
            pairs = cKDTree(self.positions).query_pairs(cutoff, output_type="ndarray")
            if len(pairs):
                d = np.linalg.norm(self.positions[pairs[:, 0]] - self.positions[pairs[:, 1]], axis=1)
                limits = 1.16 * (self.radii[pairs[:, 0]] + self.radii[pairs[:, 1]])
                self.bonds = pairs[(d < limits) & (d > 0.38)]
        self.update()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag = event.position()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)

    def mouseMoveEvent(self, event):
        if self._drag is None:
            return
        delta = event.position() - self._drag
        self._drag = event.position()
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            self.rollRequested.emit(0.55 * delta.x())
        else:
            self.orbitRequested.emit(0.55 * delta.x(), 0.55 * delta.y())

    def mouseReleaseEvent(self, event):
        self._drag = None
        self.setCursor(Qt.CursorShape.OpenHandCursor)

    def mouseDoubleClickEvent(self, event):
        self.fit_view()

    def wheelEvent(self, event):
        limit = 1.0 if self.auto_fit else 5.0
        self.zoom = float(np.clip(self.zoom * 1.15 ** (event.angleDelta().y() / 120), 0.25, limit))
        self.update()
        event.accept()

    def paintEvent(self, event):
        p = QPainter(self)
        w, h = self.width(), self.height()
        p.fillRect(self.rect(), Qt.GlobalColor.white)
        if self.crystal is None or self.frame is None:
            return
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        # The detector faces the specimen. The right-handed beam points toward the viewer.
        basis = np.array([self.frame.right, -self.frame.up, self.frame.beam])
        projected = self.positions @ basis.T
        corners = (self.corners - self.center) @ basis.T
        scale = self.pixels_per_angstrom()
        center = np.array([w / 2, h / 2 - 1])
        screen = center + projected[:, :2] * scale
        cell = center + corners[:, :2] * scale
        p.setPen(pen("#968470", 1.15))
        for i in range(8):
            for bit in (1, 2, 4):
                if not (i & bit):
                    p.drawLine(point(cell[i]), point(cell[i | bit]))
        if self.show_bonds:
            p.setPen(pen("#928775", min(5, max(1, 0.07 * scale))))
            for a, b in self.bonds:
                p.drawLine(point(screen[a]), point(screen[b]))
        draw_order = np.argsort(projected[:, 2])
        for i in draw_order:
            x, y = screen[i]
            radius = float(np.clip(self.radii[i] * 0.30 * scale, 1.5, 30))
            if x + radius < 0 or y + radius < 0 or x - radius > w or y - radius > h:
                continue
            base = self.colors[int(self.species[i])]
            if len(draw_order) < 1000:
                shading = QRadialGradient(x - 0.30 * radius, y - 0.34 * radius, 1.45 * radius)
                shading.setColorAt(0, base.lighter(151))
                shading.setColorAt(0.56, base)
                shading.setColorAt(1, base.darker(180))
                p.setBrush(shading)
            else:
                p.setBrush(base)
            p.setOpacity(max(0.24, min(1.0, float(self.occupancies[i]))))
            p.setPen(pen(base.darker(162), 0.75))
            p.drawEllipse(QPointF(x, y), radius, radius)
        p.setOpacity(1)
        scale_bar(p, w, h, scale, "Å", "#594f47")
        origin = np.array([38.0, h - 34.0])
        for axis, name, color in zip(self.crystal.cell, "abc", ("#94452f", "#7e651d", "#34577d")):
            vector = axis @ basis.T / np.linalg.norm(axis)
            end = origin + vector[:2] * 19
            p.setPen(pen(color, 1.8))
            p.drawLine(point(origin), point(end))
            if np.linalg.norm(vector[:2]) < 0.2:
                p.setBrush(QColor(color))
                p.drawEllipse(point(end), 2.4, 2.4)
            p.setFont(font(11, True))
            p.drawText(point(end + [3, -2]), name)
        p.setFont(font(11, True))
        x = 13
        for i, (number, color) in enumerate(self.colors.items()):
            if x > w - 55:
                p.setPen(QColor("#4e463e"))
                p.drawText(x, 18, f"+{len(self.colors) - i}")
                break
            p.setPen(pen(color.darker(150), 0.6))
            p.setBrush(color)
            p.drawEllipse(QPointF(x + 3, 14), 3.5, 3.5)
            text = chemical_symbols[number]
            p.setPen(QColor("#4e463e"))
            p.drawText(x + 11, 18, text)
            x += 24 + p.fontMetrics().horizontalAdvance(text)
        if self.sampled:
            p.setFont(font(11))
            p.drawText(12, 34, "Showing a subset of sites")


class PatternView(QWidget):
    reflectionSelected = Signal(int)
    zoomChanged = Signal(float)

    def __init__(self):
        super().__init__()
        from .detector_view import DetectorOptions, DetectorRaster
        self.setObjectName("patternView")
        self.setAccessibleName("Diffraction pattern")
        self.setAccessibleDescription("Scroll to zoom around the pointer. Double-click to fit the pattern.")
        self.setMinimumSize(220, 180)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.pattern = None
        self.g_max = 1.25
        self.zoom = 1.0
        # Detector coordinates at the viewport center, in units of g_max,
        # with horizontal and vertical components following the screen axes.
        self._focus = np.zeros(2)
        self.show_labels = False
        self.show_grid = True
        self.mode = "spots"
        self.grids = {"spots": True, "detector": False}
        self.detector_options = DetectorOptions()
        self.raster = DetectorRaster()
        self.selected = -1
        self.message = "Calculating…"
        self._screen = np.empty((0, 2))
        self._indices = np.array([], dtype=int)

    def set_pattern(self, pattern, g_max):
        self.pattern = pattern
        self.raster.exposure_key = None
        self._indices = np.array([], dtype=int)
        self.g_max = g_max
        self.message = ""
        self.selected = -1
        self.update()

    def set_message(self, text):
        self.pattern = None
        self.message = text
        self._indices = np.array([], dtype=int)
        self.update()

    def set_mode(self, mode):
        self.grids[self.mode] = self.show_grid
        self.mode = mode
        self.show_grid = self.grids[mode]
        self._indices = np.array([], dtype=int)
        self.update()

    def set_detector_options(self, options):
        self.detector_options = options
        self.update()

    def set_selected(self, index):
        self.selected = index
        self.update()

    def _plot_height(self):
        return self.height() - 45 if self.mode == "detector" else self.height()

    def _projection(self):
        height = self._plot_height()
        radius = max(20, min(self.width() - 42, height - 32) / 2) * self.zoom
        midpoint = np.array([self.width() / 2, height / 2 - 1])
        return midpoint - radius * self._focus, radius

    def _set_zoom(self, value, anchor=None):
        value = float(np.clip(value, 0.25, 16.0))
        if value == self.zoom:
            return
        if anchor is not None:
            origin, radius = self._projection()
            fixed = (np.array([anchor.x(), anchor.y()]) - origin) / radius
            self._focus = fixed + (self._focus - fixed) * (self.zoom / value)
        # Retain the full frame at fit and leave part of the detector in view
        # when scrolling over the margin surrounding the square detector.
        limit = max(0.0, 1.0 - 1.0 / value)
        self._focus = np.clip(self._focus, -limit, limit)
        self.zoom = value
        self._indices = np.array([], dtype=int)
        self.zoomChanged.emit(value)
        self.update()

    def zoom_in(self):
        self._set_zoom(self.zoom * 1.25)

    def zoom_out(self):
        self._set_zoom(self.zoom / 1.25)

    def reset_zoom(self):
        self._set_zoom(1.0)

    def wheelEvent(self, event):
        if event.position().y() >= self._plot_height():
            event.ignore()
            return
        pixels = event.pixelDelta().y()
        steps = pixels / 40.0 if pixels else event.angleDelta().y() / 120.0
        if steps:
            self._set_zoom(self.zoom * 1.25 ** np.clip(steps, -24, 24), event.position())
        event.accept()

    def mouseDoubleClickEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton and event.position().y() < self._plot_height():
            self.reset_zoom()
            event.accept()
        else:
            super().mouseDoubleClickEvent(event)

    def nearest(self, pos):
        if not len(self._indices) or self.pattern is None or pos.y() >= self._plot_height():
            return -1
        delta = self._screen[self._indices] - [pos.x(), pos.y()]
        distances = np.einsum("ij,ij->i", delta, delta)
        matches = self._indices[distances < 81]
        if not len(matches):
            return -1
        closest = np.min(distances)
        # HOLZ reflections can project onto the same detector position.
        matches = self._indices[distances <= closest + 0.8]
        return int(matches[np.argmax(self.pattern.intensity[matches])])

    def mouseMoveEvent(self, event):
        index = self.nearest(event.position())
        details = []
        if (self.mode == "detector" and event.position().y() < self._plot_height()
                and self.pattern is not None and self.raster.exposure_key is not None):
            center, radius = self._projection()
            q = (np.array([event.position().x(), event.position().y()]) - center) / radius
            size = self.raster.exposure.shape[0]
            x, y = np.rint((q + 1) * (size - 1) / 2).astype(int)
            if 0 <= x < size and 0 <= y < size:
                details.append(f"Detector pixel    I / I₀ = {self.raster.exposure[y, x]:.4g}")
        if index >= 0:
            h, k, l = self.pattern.hkl[index]
            details.append(f"({h} {k} {l})    I / I₀ = {self.pattern.intensity[index]:.4g}")
        if not details:
            QToolTip.hideText()
            self.unsetCursor()
            return
        self.setCursor(Qt.CursorShape.PointingHandCursor if index >= 0 else Qt.CursorShape.CrossCursor)
        QToolTip.showText(event.globalPosition().toPoint(), "\n".join(details), self)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            i = self.nearest(event.position())
            if i >= 0:
                self.reflectionSelected.emit(i)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w = self.width()
        detector = self.mode == "detector"
        plot_height = self._plot_height()
        plot_rect = QRectF(0, 0, w, plot_height)
        center, radius = self._projection()
        scale = radius / self.g_max
        p.fillRect(self.rect(), QColor("#202735"))
        p.setClipRect(plot_rect)
        ink = "#e2daca"
        grid_color = QColor("#3c4657")
        if detector and self.pattern is not None:
            self.raster.prepare(self.pattern, self.g_max, self.detector_options)
            background = QColor.fromRgba(self.raster.table[0])
            p.fillRect(QRectF(0, 0, w, plot_height), background)
            # The first and last pixel centres coincide with ±g_max.
            half_pixel = radius / (self.raster.image.width() - 1)
            side = 2 * (radius + half_pixel)
            target = QRectF(center[0] - radius - half_pixel, center[1] - radius - half_pixel, side, side)
            p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
            p.drawImage(target, self.raster.image)
            ink = "#281f1b" if background.lightnessF() > 0.5 else "#f3ecdf"
            grid_color = QColor(ink)
            grid_color.setAlpha(75)
        if self.show_grid:
            p.setPen(QPen(grid_color, 1))
            p.drawLine(QPointF(0, center[1]), QPointF(w, center[1]))
            p.drawLine(QPointF(center[0], 0), QPointF(center[0], plot_height))
            step = scale_length(self.g_max * 0.65 / max(1, self.zoom))
            for q in np.arange(step, self.g_max * 1.02, step):
                p.drawEllipse(point(center), q * scale, q * scale)
        if self.pattern is None:
            scale_bar(p, w, plot_height, scale, "Å⁻¹", ink)
            p.setPen(QColor("#f3ecdf"))
            p.setFont(font(13))
            p.drawText(plot_rect, Qt.AlignmentFlag.AlignCenter, self.message)
            return
        pattern = self.pattern
        screen = center + pattern.xy * np.array([scale, -scale])
        self._screen = screen
        intensities = np.maximum(pattern.intensity, 0.0)
        off_axis = np.linalg.norm(pattern.xy, axis=1) > 1e-8
        peak = np.max(intensities[off_axis], initial=0)
        # A fourth-root display brings out weaker diffracted beams. Table values are linear.
        reference = max(float(peak), 1e-14)
        strength = np.power(np.clip(intensities / reference, 0, 1), 0.25)
        visible = ((screen[:, 0] >= -5) & (screen[:, 0] <= w + 5)
                   & (screen[:, 1] >= -5) & (screen[:, 1] <= plot_height + 5))
        if detector:
            self._indices = np.flatnonzero(visible & (intensities > 0))
        else:
            self._indices = np.flatnonzero(visible & (strength >= 0.023) & (intensities > 1e-16))
        order = self._indices[np.argsort(intensities[self._indices], kind="stable")]
        if not detector:
            self._draw_spots(p, screen, order, strength)
        if self.show_labels:
            p.setFont(font(10))
            boxes = []
            for i in reversed(order):
                if len(boxes) >= 24:
                    break
                x, y = screen[i]
                hkl = " ".join(str(v) for v in pattern.hkl[i])
                box = QRectF(x + 5, y - 17, p.fontMetrics().horizontalAdvance(hkl) + 4, 14)
                if not plot_rect.contains(box) or any(box.intersects(b) for b in boxes):
                    continue
                p.setPen(QColor(ink))
                p.drawText(box, Qt.AlignmentFlag.AlignCenter, hkl)
                boxes.append(box)
        if 0 <= self.selected < len(screen):
            p.setPen(pen("#e49a50", 1.5))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawEllipse(point(screen[self.selected]), 8, 8)
        scale_bar(p, w, plot_height, scale, "Å⁻¹", ink)
        if detector:
            p.setClipping(False)
            self._draw_colorbar(p, plot_height)

    def _draw_colorbar(self, painter, top):
        painter.save()
        w = self.width()
        painter.setPen(QColor("#eee5d7"))
        painter.setFont(font(10))
        mode = "log" if self.detector_options.logarithmic else "linear"
        painter.drawText(QRectF(9, top + 1, w - 18, 13), Qt.AlignmentFlag.AlignCenter,
                         f"I / I₀ per detector pixel  ·  {mode}")
        inset = 14
        painter.drawImage(QRectF(inset, top + 16, w - inset * 2, 8), self.raster.colorbar)
        low, high = self.raster.lower, self.raster.upper
        midpoint = np.sqrt(low * high) if self.detector_options.logarithmic else (low + high) / 2
        rect = QRectF(inset, top + 26, w - 2 * inset, 16)
        for value, align in ((low, Qt.AlignmentFlag.AlignLeft),
                             (midpoint, Qt.AlignmentFlag.AlignHCenter),
                             (high, Qt.AlignmentFlag.AlignRight)):
            painter.drawText(rect, align | Qt.AlignmentFlag.AlignTop, f"{value:.2g}")
        painter.restore()

    @staticmethod
    def _draw_spots(painter, screen, order, strength):
        from PySide6.QtGui import QPolygonF

        # Draw points in groups of similar brightness so a large basis does not
        # require thousands of brush changes and individual Qt paint calls.
        bins = np.rint(strength[order] * 255).astype(np.int16)
        for group in np.split(np.arange(len(order)), np.flatnonzero(np.diff(bins)) + 1):
            if not len(group):
                continue
            s = bins[group[0]] / 255
            color = QColor.fromRgbF(0.72 + 0.27*s, 0.69 + 0.29*s, 0.67 + 0.27*s,
                                    0.16 + 0.84*s)
            pen = QPen(color, 2 * (1.0 + 3.4*s))
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            painter.drawPoints(QPolygonF([point(xy) for xy in screen[order[group]]]))
