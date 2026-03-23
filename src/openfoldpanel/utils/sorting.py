"""Sorting helpers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

_NATURAL_RE = re.compile(r"(\d+)")


def natural_sort_key(value: str | Path) -> list[Any]:
    """Return a mixed string/integer key for natural ordering."""

    text = str(value)
    parts = _NATURAL_RE.split(text)
    key: list[Any] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part.lower())
    return key
