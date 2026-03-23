"""Secure archive extraction with zip-slip and tar traversal protection."""

from __future__ import annotations

import tarfile
import zipfile
from pathlib import Path

from openfoldpanel.extractors.validators import is_supported_archive


def _resolve_destination(base_dir: Path, member_name: str) -> Path:
    target = (base_dir / member_name).resolve()
    base = base_dir.resolve()
    if not str(target).startswith(str(base)):
        raise ValueError(f"Refusing to extract path outside destination: {member_name}")
    return target


def extract_archive(archive_path: Path, destination: Path) -> list[Path]:
    """Extract a supported archive into ``destination`` safely."""

    if not is_supported_archive(archive_path):
        raise ValueError(f"Unsupported archive type: {archive_path}")

    destination.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    archive_name = archive_path.name.lower()

    if archive_name.endswith(".zip"):
        with zipfile.ZipFile(archive_path) as zf:
            for info in zf.infolist():
                target = _resolve_destination(destination, info.filename)
                if info.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(info, "r") as src, target.open("wb") as dst:
                    dst.write(src.read())
                extracted.append(target)
        return extracted

    mode = "r"
    if archive_name.endswith((".tar.gz", ".tgz")):
        mode = "r:gz"
    elif archive_name.endswith((".tar.bz2", ".tbz2")):
        mode = "r:bz2"
    elif archive_name.endswith((".tar.xz", ".txz")):
        mode = "r:xz"

    with tarfile.open(archive_path, mode) as tf:
        for member in tf.getmembers():
            target = _resolve_destination(destination, member.name)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            if not member.isfile():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            extracted_file = tf.extractfile(member)
            if extracted_file is None:
                continue
            with extracted_file, target.open("wb") as dst:
                dst.write(extracted_file.read())
            extracted.append(target)
    return extracted
