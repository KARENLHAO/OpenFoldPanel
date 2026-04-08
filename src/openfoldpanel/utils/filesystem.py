"""Filesystem and temporary-directory helpers."""

from __future__ import annotations

import shutil
from pathlib import Path


def ensure_directory(path: Path) -> Path:
    """Create a directory if it does not already exist."""

    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_rmtree(path: Path) -> None:
    """Delete a temporary directory only if it exists."""

    if path.exists():
        shutil.rmtree(path, ignore_errors=True)
