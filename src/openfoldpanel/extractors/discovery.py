"""Job discovery for structure files and archives."""

from __future__ import annotations

import logging
from pathlib import Path

from openfoldpanel.extractors.validators import is_supported_structure_file
from openfoldpanel.models import JobDefinition
from openfoldpanel.utils.sorting import natural_sort_key


def _collect_structure_files(root: Path) -> tuple[list[Path], list[str]]:
    structure_files: list[Path] = []
    ignored_files: list[str] = []
    for path in sorted(root.rglob("*"), key=natural_sort_key):
        if not path.is_file():
            continue
        if is_supported_structure_file(path):
            structure_files.append(path)
        else:
            ignored_files.append(str(path.relative_to(root)))
    return structure_files, ignored_files


def discover_jobs_from_structure(path: Path) -> list[JobDefinition]:
    """Treat a single structure file as a single job."""

    return [
        JobDefinition(
            name=path.stem,
            root_dir=path.parent,
            structure_files=[path],
            ignored_files=[],
        )
    ]


def discover_jobs_from_extracted_root(root: Path, logger: logging.Logger) -> list[JobDefinition]:
    """Discover jobs following the archive rules described in the README."""

    first_level_dirs = sorted(
        [child for child in root.iterdir() if child.is_dir()],
        key=lambda item: natural_sort_key(item.name),
    )
    root_files = sorted(
        [child for child in root.iterdir() if child.is_file()],
        key=lambda item: natural_sort_key(item.name),
    )

    jobs: list[JobDefinition] = []

    if first_level_dirs:
        for directory in first_level_dirs:
            structure_files, ignored_files = _collect_structure_files(directory)
            if not structure_files:
                if ignored_files:
                    logger.info("Ignoring archive directory without structures: %s", directory.name)
                continue
            jobs.append(
                JobDefinition(
                    name=directory.name,
                    root_dir=directory,
                    structure_files=sorted(structure_files, key=natural_sort_key),
                    ignored_files=ignored_files,
                )
            )
        if jobs:
            for file_path in root_files:
                logger.info("Ignoring root-level file because archive is directory-batched: %s", file_path.name)
            return jobs

    structure_files = [path for path in root_files if is_supported_structure_file(path)]
    ignored_files = [path.name for path in root_files if not is_supported_structure_file(path)]
    for ignored in ignored_files:
        logger.info("Ignoring non-structure file in archive root: %s", ignored)
    if structure_files:
        jobs.append(
            JobDefinition(
                name=root.name,
                root_dir=root,
                structure_files=sorted(structure_files, key=natural_sort_key),
                ignored_files=ignored_files,
            )
        )
    return jobs
