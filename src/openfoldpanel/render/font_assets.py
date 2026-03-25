"""Helpers for embedding packaged font assets into HTML and SVG outputs."""

from __future__ import annotations

import base64
from functools import lru_cache
from importlib import resources

from openfoldpanel.constants import EMBEDDED_TIMES_NEW_ROMAN_ALIAS


UI_PACKAGE = "openfoldpanel.UI"
TIMES_NEW_ROMAN_RESOURCE = "fonts/Times New Roman.ttf"


@lru_cache(maxsize=1)
def embedded_times_new_roman_css() -> str:
    """Return a reusable @font-face rule for the packaged Times New Roman asset."""

    encoded_font = base64.b64encode(
        resources.files(UI_PACKAGE).joinpath(TIMES_NEW_ROMAN_RESOURCE).read_bytes()
    ).decode("ascii")
    return "\n".join(
        [
            "@font-face {",
            f'  font-family: "{EMBEDDED_TIMES_NEW_ROMAN_ALIAS}";',
            f'  src: url("data:font/ttf;base64,{encoded_font}") format("truetype");',
            "  font-style: normal;",
            "  font-weight: 400;",
            "  font-display: swap;",
            "}",
        ]
    )
