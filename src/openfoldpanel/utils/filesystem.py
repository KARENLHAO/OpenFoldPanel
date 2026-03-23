"""Filesystem and temporary-directory helpers."""

from __future__ import annotations

import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator


def ensure_directory(path: Path) -> Path:
    """Create a directory if it does not already exist."""

    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_rmtree(path: Path) -> None:
    """Delete a temporary directory only if it exists."""

    if path.exists():
        shutil.rmtree(path, ignore_errors=True)


@contextmanager
def temporary_workspace(prefix: str = "openfoldpanel_") -> Iterator[Path]:
    """Yield a temporary directory path and clean it up afterwards."""

    temp_dir = Path(tempfile.mkdtemp(prefix=prefix))
    try:
        yield temp_dir
    finally:
        safe_rmtree(temp_dir)
