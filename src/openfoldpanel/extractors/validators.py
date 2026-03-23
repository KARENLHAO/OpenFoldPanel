"""Validation helpers for supported inputs."""

from __future__ import annotations

from pathlib import Path

from openfoldpanel.constants import SUPPORTED_ARCHIVE_SUFFIXES, SUPPORTED_STRUCTURE_SUFFIXES


def normalize_lower_name(path: Path) -> str:
    """Return a lowercased filename for suffix matching."""

    return path.name.lower()


def is_supported_structure_file(path: Path) -> bool:
    """Check whether the file looks like a supported structure file."""

    suffixes = path.suffixes
    if not suffixes:
        return False
    lower = normalize_lower_name(path)
    return any(lower.endswith(ext) for ext in SUPPORTED_STRUCTURE_SUFFIXES)


def is_supported_archive(path: Path) -> bool:
    """Check whether the file looks like a supported archive."""

    lower = normalize_lower_name(path)
    return any(lower.endswith(ext) for ext in SUPPORTED_ARCHIVE_SUFFIXES)


def validate_input_path(path: Path) -> None:
    """Raise a readable error for unsupported inputs."""

    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {path}")
    if path.is_dir():
        raise ValueError(f"Input path must be a structure file or archive, not a directory: {path}")
    if not is_supported_structure_file(path) and not is_supported_archive(path):
        raise ValueError(f"Unsupported input type: {path}")
