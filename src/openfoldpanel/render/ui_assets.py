"""Cached accessors for packaged HTML UI assets."""

from __future__ import annotations

from functools import lru_cache
from importlib import resources

from openfoldpanel.render.font_assets import embedded_times_new_roman_css

UI_PACKAGE = "openfoldpanel.UI"
UI_STYLE_FILES = (
    "styles/tokens.css",
    "styles/base.css",
    "styles/layout.css",
    "styles/components.css",
    "styles/figure.css",
    "styles/atmosphere.css",
)
UI_SCRIPT_FILE = "scripts/report.js"
UI_TEMPLATE_FILE = "report.template.html"


def load_ui_template() -> str:
    """Return the packaged report HTML template."""

    return _load_ui_resource(UI_TEMPLATE_FILE)


@lru_cache(maxsize=1)
def load_ui_styles() -> str:
    """Return all packaged CSS concatenated into one inline stylesheet."""

    return "\n\n".join([embedded_times_new_roman_css(), *(_load_ui_resource(path) for path in UI_STYLE_FILES)])


@lru_cache(maxsize=1)
def load_ui_script() -> str:
    """Return the packaged report JavaScript."""

    return _load_ui_resource(UI_SCRIPT_FILE)


@lru_cache(maxsize=None)
def _load_ui_resource(relative_path: str) -> str:
    return resources.files(UI_PACKAGE).joinpath(relative_path).read_text(encoding="utf-8")
