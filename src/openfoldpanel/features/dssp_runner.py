"""Run DSSP externally and parse per-residue annotations."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from openfoldpanel.utils.subprocess import ExternalToolError, MissingExecutableError, run_command, which


@dataclass(slots=True)
class DSSPResidueFeature:
    chain: str
    seq_id: int
    insertion_code: str
    amino_acid: str
    dssp_code: str
    asa: float | None


def find_dssp_executable() -> str | None:
    """Return a working DSSP executable name if available."""

    for executable in ("mkdssp", "dssp"):
        if which(executable):
            return executable
    return None


def run_dssp(
    structure_path: Path,
    logger: logging.Logger,
    *,
    display_name: str | None = None,
) -> tuple[dict[tuple[str, int, str], DSSPResidueFeature], list[str]]:
    """Execute DSSP and parse the result; degrade gracefully on errors."""

    executable = find_dssp_executable()
    target_name = display_name or structure_path.name
    if executable is None:
        return {}, ["DSSP executable not found; using geometric fallback for secondary structure and accessibility."]

    commands_to_try = [
        [executable, "--output-format", "dssp", str(structure_path)],
        [executable, str(structure_path)],
        [executable, "-i", str(structure_path)],
    ]
    output = ""
    last_error: Exception | None = None
    for command in commands_to_try:
        try:
            result = run_command(command, check=True)
            output = result.stdout
            if output.strip():
                break
        except (MissingExecutableError, ExternalToolError) as exc:
            last_error = exc
            logger.debug("DSSP command variant failed for %s: %s", target_name, " ".join(command))
            continue

    if not output.strip():
        if last_error is not None:
            logger.warning("DSSP failed for %s: %s", target_name, last_error)
        else:
            logger.warning("DSSP produced no output for %s", target_name)
        return {}, [f"DSSP failed for {target_name}; using geometric fallback for secondary structure and accessibility."]

    try:
        parsed = parse_dssp_output(_trim_to_classic_dssp(output))
    except Exception as exc:
        logger.warning("Unable to parse DSSP output for %s: %s", target_name, exc)
        return {}, [f"DSSP output could not be parsed for {target_name}; using geometric fallback."]
    return parsed, []


def _trim_to_classic_dssp(output: str) -> str:
    """Trim mixed-output DSSP stdout down to the classic DSSP text section."""

    lines = output.splitlines()
    for index, line in enumerate(lines):
        if line.startswith("==== Secondary Structure Definition by the program DSSP"):
            return "\n".join(lines[index:])
    return output


def parse_dssp_output(output: str) -> dict[tuple[str, int, str], DSSPResidueFeature]:
    """Parse classic DSSP fixed-width text output."""

    in_body = False
    features: dict[tuple[str, int, str], DSSPResidueFeature] = {}
    for line in output.splitlines():
        if not in_body:
            if line.strip().startswith("#"):
                in_body = True
            continue
        if len(line) < 40:
            continue
        seq_text = line[5:10].strip()
        if not seq_text or "!" in seq_text:
            continue
        insertion_code = line[10].strip()
        chain = line[11].strip() or "_"
        amino_acid = line[13].strip() or "X"
        dssp_code = (line[16].strip() or "C").upper()
        asa_text = line[34:38].strip()
        seq_id = int(seq_text)
        asa = float(asa_text) if asa_text else None
        features[(chain, seq_id, insertion_code)] = DSSPResidueFeature(
            chain=chain,
            seq_id=seq_id,
            insertion_code=insertion_code,
            amino_acid=amino_acid,
            dssp_code=dssp_code,
            asa=asa,
        )
    return features
