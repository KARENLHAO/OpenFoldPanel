"""Optional multiple-sequence alignment via Clustal Omega."""

from __future__ import annotations

import logging
from pathlib import Path

from openfoldpanel.models import MSARow
from openfoldpanel.utils.subprocess import ExternalToolError, MissingExecutableError, run_command, which


def write_fasta(rows: list[tuple[str, str]], path: Path) -> None:
    """Write a FASTA file from identifier/sequence rows."""

    with path.open("w", encoding="utf-8") as handle:
        for identifier, sequence in rows:
            handle.write(f">{identifier}\n{sequence}\n")


def align_sequences(rows: list[tuple[str, str]], output_path: Path, logger: logging.Logger) -> tuple[list[MSARow], list[str]]:
    """Align sequences with clustalo when available."""

    if which("clustalo") is None:
        return [], ["clustalo was not found; MSA alignment was skipped."]

    input_fasta = output_path.parent / "msa_input.fasta"
    write_fasta(rows, input_fasta)
    try:
        run_command(
            [
                "clustalo",
                "-i",
                str(input_fasta),
                "-o",
                str(output_path),
                "--force",
                "--outfmt=fasta",
            ]
        )
    except (MissingExecutableError, ExternalToolError) as exc:
        logger.warning("clustalo alignment failed: %s", exc)
        return [], ["clustalo alignment failed; MSA track was skipped."]
    return read_fasta_alignment(output_path), []


def read_fasta_alignment(path: Path) -> list[MSARow]:
    """Read a FASTA alignment file into MSA rows."""

    rows: list[MSARow] = []
    identifier: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if identifier is not None:
                rows.append(MSARow(identifier=identifier, sequence="".join(chunks), is_query=identifier == "query"))
            identifier = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if identifier is not None:
        rows.append(MSARow(identifier=identifier, sequence="".join(chunks), is_query=identifier == "query"))
    return rows
