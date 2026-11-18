"""Color themes for the Qt controls surrounding the scientific views."""

from dataclasses import asdict, dataclass
from string import Template

from PySide6.QtGui import QColor, QPalette


@dataclass(frozen=True)
class _Theme:
    name: str
    workspace: str
    panel: str
    field: str
    alternate: str
    button: str
    hover: str
    border: str
    text: str
    muted: str
    disabled: str
    heading: str
    banner: str
    trim: str
    banner_text: str
    accent: str
    accent_hover: str
    on_accent: str
    selection: str
    on_selection: str
    error: str


_COLORS = {
    "vintage": _Theme(
        name="Vintage", workspace="#e9e2d5", panel="#fcf8f0", field="#fffdf8",
        alternate="#f2eadc", button="#faf5e9", hover="#efe5d4", border="#b7a68f",
        text="#302f38", muted="#655d53", disabled="#70685d", heading="#293449",
        banner="#293449", trim="#d7b07a", banner_text="#fff4dd", accent="#90482f",
        accent_hover="#773d29", on_accent="#fff8e9", selection="#e7d3b5",
        on_selection="#302b33", error="#923f29",
    ),
    "white": _Theme(
        name="White", workspace="#ffffff", panel="#ffffff", field="#ffffff",
        alternate="#f4f4f4", button="#fafafa", hover="#eaeaea", border="#909090",
        text="#242424", muted="#545454", disabled="#686868", heading="#171717",
        banner="#ffffff", trim="#595959", banner_text="#202020", accent="#353535",
        accent_hover="#171717", on_accent="#ffffff", selection="#dedede",
        on_selection="#181818", error="#a13632",
    ),
    "studio_light": _Theme(
        name="Studio Light", workspace="#e1e8f1", panel="#f5f8fc", field="#ffffff",
        alternate="#e8eff7", button="#f0f5fa", hover="#dbe8f6", border="#9aacc1",
        text="#1e2e40", muted="#475d73", disabled="#5b6a7b", heading="#163553",
        banner="#183854", trim="#a6cdec", banner_text="#f5faff", accent="#245f94",
        accent_hover="#194c79", on_accent="#ffffff", selection="#c5def5",
        on_selection="#162c43", error="#a13738",
    ),
    "midnight": _Theme(
        name="Midnight", workspace="#141a26", panel="#202838", field="#182131",
        alternate="#2a3446", button="#303c50", hover="#3b4960", border="#607189",
        text="#edf1f7", muted="#bfccda", disabled="#9aa7b8", heading="#eef3fa",
        banner="#121927", trim="#e6b878", banner_text="#fff3de", accent="#e6b878",
        accent_hover="#f6cc93", on_accent="#261c10", selection="#455574",
        on_selection="#f4f7ff", error="#ffad98",
    ),
    "high_contrast": _Theme(
        name="High Contrast", workspace="#000000", panel="#111111", field="#000000",
        alternate="#242424", button="#191919", hover="#333333", border="#c5c5c5",
        text="#ffffff", muted="#ececec", disabled="#adadad", heading="#ffffff",
        banner="#000000", trim="#ffe36b", banner_text="#ffffff", accent="#ffe36b",
        accent_hover="#fff29f", on_accent="#000000", selection="#ffe36b",
        on_selection="#000000", error="#ffb8bf",
    ),
}
THEMES = {theme_id: theme.name for theme_id, theme in _COLORS.items()}
DEFAULT_THEME = "white"


_STYLE = Template("""
QWidget { font-size: 13px; }
QMainWindow, QWidget#workspace, QSplitter, QSplitter::handle { background: $workspace; }
QFrame#panel, QWidget#sidebar, QScrollArea#controls { background: $panel; }
QFrame#banner { background: $banner; color: $trim; border-bottom: 3px solid $trim; }
QFrame#panel { border: 1px solid $border; border-radius: 5px; }
QScrollArea#controls { border: 0; border-right: 1px solid $border; }
QLabel#brand { color: $banner_text; font-size: 32px; font-weight: 600; }
QLabel#file { color: $banner_text; font-size: 14px; font-weight: 600; }
QLabel#section { font-size: 21px; font-weight: 600; color: $heading; }
QLabel#group { font-size: 14px; font-weight: 650; color: $accent; }
QLabel#subtle, QLabel#count, QLabel#readout { color: $muted; font-size: 12px; }
QLabel#readout { font-size: 11px; }
QFrame#banner QLabel#subtle { color: $trim; }
QLabel#error { color: $error; font-size: 12px; }
QFrame#rule { background: $border; border: 0; }
QFrame#banner QFrame#rule { background: $trim; }
QPushButton, QToolButton {
    background: $button; color: $text; border: 1px solid $border; border-radius: 4px;
    padding: 5px 11px; font-weight: 550; min-height: 16px;
}
QPushButton:hover, QToolButton:hover { background: $hover; border-color: $accent; }
QPushButton:pressed, QToolButton:pressed, QToolButton:checked {
    background: $selection; color: $on_selection;
}
QPushButton#simulateButton { background: $accent; color: $on_accent; border-color: $accent; padding: 8px; }
QPushButton#simulateButton:hover { background: $accent_hover; }
QPushButton#openButton { background: $banner_text; border-color: $trim; color: $banner; }
QPushButton#openButton:hover { background: $trim; }
QToolButton#themeButton { background: $banner; color: $banner_text; border-color: $trim; }
QToolButton#themeButton:hover, QToolButton#themeButton:pressed {
    background: $trim; color: $banner;
}
QPushButton:focus, QToolButton:focus { outline: 2px solid $accent; }
QToolButton#themeButton:focus, QPushButton#openButton:focus { outline: 2px solid $trim; }
QToolButton::menu-indicator { width: 0; }
QComboBox, QLineEdit, QAbstractSpinBox {
    background: $field; color: $text; border: 1px solid $border; border-radius: 4px;
    padding: 4px 6px; min-height: 19px; selection-background-color: $accent; selection-color: $on_accent;
}
QComboBox:focus, QLineEdit:focus, QAbstractSpinBox:focus { border-color: $accent; }
QLineEdit[invalid="true"] { border-color: $error; }
QCheckBox { spacing: 8px; }
QPushButton:disabled, QToolButton:disabled, QComboBox:disabled, QAbstractSpinBox:disabled,
QLineEdit:disabled { color: $disabled; background: $alternate; border-color: $border; }
QCheckBox:disabled, QLabel:disabled, QMenu::item:disabled { color: $disabled; }
QSplitter::handle:hover { background: $border; }
QHeaderView::section {
    color: $text; background: $alternate; border: 0; border-bottom: 1px solid $border;
    padding: 4px 10px; font-size: 12px; font-weight: 600;
}
QTableView {
    background: $field; alternate-background-color: $alternate; color: $text; border: 0;
    selection-background-color: $selection; selection-color: $on_selection; outline: 0;
}
QStatusBar { background: $workspace; color: $muted; border-top: 1px solid $border; }
QStatusBar::item { border: 0; }
QMenuBar, QMenu { background: $panel; color: $text; }
QMenu { border: 1px solid $border; }
QMenuBar::item:selected, QMenu::item:selected { background: $selection; color: $on_selection; }
QComboBox QAbstractItemView {
    background: $field; color: $text; selection-background-color: $accent; selection-color: $on_accent;
}
QToolTip { background: $banner; color: $banner_text; border: 1px solid $trim; padding: 5px; }
""")


def apply_theme(app, theme_id):
    """Apply the selected palette to existing and future windows. Return its ID."""
    if theme_id not in _COLORS:
        theme_id = DEFAULT_THEME
    theme = _COLORS[theme_id]
    # Fusion also uses these roles to paint native indicators, arrows and dialogs.
    role = QPalette.ColorRole
    palette = QPalette()
    for color, roles in (
        (theme.panel, (role.Window,)),
        (theme.field, (role.Base, role.Light)),
        (theme.alternate, (role.AlternateBase,)),
        (theme.text, (role.Text, role.WindowText, role.ButtonText)),
        (theme.button, (role.Button,)),
        (theme.hover, (role.Midlight,)),
        (theme.border, (role.Mid, role.Dark)),
        (theme.banner, (role.ToolTipBase, role.Shadow)),
        (theme.banner_text, (role.ToolTipText, role.BrightText)),
        (theme.muted, (role.PlaceholderText,)),
        (theme.accent, (role.Highlight, role.Accent, role.Link)),
        (theme.accent_hover, (role.LinkVisited,)),
        (theme.on_accent, (role.HighlightedText,)),
    ):
        for palette_role in roles:
            palette.setColor(palette_role, QColor(color))
    for palette_role in (role.Text, role.WindowText, role.ButtonText, role.PlaceholderText):
        palette.setColor(QPalette.ColorGroup.Disabled, palette_role, QColor(theme.disabled))
    app.setPalette(palette)
    app.setStyleSheet(_STYLE.substitute(asdict(theme)))
    return theme_id
