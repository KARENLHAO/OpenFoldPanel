"""Optional homolog search wrappers for BLAST+ or MMseqs2."""

from __future__ import annotations

import logging
from pathlib import Path

from openfoldpanel.utils.subprocess import ExternalToolError, MissingExecutableError, run_command, which


def search_homologs(
    query_fasta: Path,
    database: Path,
    *,
    max_hits: int,
    workdir: Path,
    logger: logging.Logger,
) -> tuple[list[tuple[str, str]], list[str]]:
    """Search a local sequence database for homologs."""

    if which("blastp"):
        return _search_with_blastp(query_fasta, database, max_hits=max_hits, logger=logger)
    if which("mmseqs"):
        return _search_with_mmseqs(query_fasta, database, max_hits=max_hits, workdir=workdir, logger=logger)
    return [], ["Neither blastp nor mmseqs was found; MSA search was skipped."]


def _search_with_blastp(
    query_fasta: Path,
    database: Path,
    *,
    max_hits: int,
    logger: logging.Logger,
) -> tuple[list[tuple[str, str]], list[str]]:
    try:
        result = run_command(
            [
                "blastp",
                "-query",
                str(query_fasta),
                "-db",
                str(database),
                "-max_target_seqs",
                str(max_hits),
                "-outfmt",
                "6 sseqid sseq",
            ]
        )
    except (MissingExecutableError, ExternalToolError) as exc:
        logger.warning("blastp search failed: %s", exc)
        return [], ["blastp search failed; MSA track was skipped."]

    rows: list[tuple[str, str]] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        rows.append((parts[0], parts[1].replace("-", "")))
    return rows, []


def _search_with_mmseqs(
    query_fasta: Path,
    database: Path,
    *,
    max_hits: int,
    workdir: Path,
    logger: logging.Logger,
) -> tuple[list[tuple[str, str]], list[str]]:
    output_tsv = workdir / "mmseqs_hits.tsv"
    tmp_dir = workdir / "mmseqs_tmp"
    try:
        run_command(
            [
                "mmseqs",
                "easy-search",
                str(query_fasta),
                str(database),
                str(output_tsv),
                str(tmp_dir),
                "--format-output",
                "target,alntseq",
                "--max-seqs",
                str(max_hits),
            ]
        )
    except (MissingExecutableError, ExternalToolError) as exc:
        logger.warning("mmseqs search failed: %s", exc)
        return [], ["mmseqs search failed; MSA track was skipped."]

    rows: list[tuple[str, str]] = []
    if output_tsv.exists():
        for line in output_tsv.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rows.append((parts[0], parts[1].replace("-", "")))
    return rows, []
