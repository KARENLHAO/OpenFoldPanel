"""Shared CSV formatting and writing helpers."""

from __future__ import annotations

import csv
from pathlib import Path

CSV_HEADER_TOKEN_OVERRIDES = {
    "uid": "UID",
    "imgt": "IMGT",
    "plddt": "pLDDT",
    "tm": "TM",
}
CSV_HEADER_NAME_OVERRIDES = {
    "aa": "Resname",
    "partner_aa": "Partner Resname",
    "one_letter": "Resname",
    "structure_name": "Structure",
    "structure_count": "Structure Count",
    "union_count": "Combine Count",
    "intersection_count": "Consensus Count",
    "combine_residue": "Combine Residue",
    "consensus_residue": "Consensus Residue",
}


def write_display_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    """Write a CSV whose header labels are humanized for downstream readers."""

    display_fieldnames = [csv_header_label(name) for name in fieldnames]
    header_map = dict(zip(fieldnames, display_fieldnames, strict=True))
    display_rows = [{header_map[key]: value for key, value in row.items()} for row in rows]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=display_fieldnames)
        writer.writeheader()
        writer.writerows(display_rows)


def write_raw_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    """Write a CSV without modifying the provided field names."""

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def csv_header_label(name: str) -> str:
    """Convert an internal field name into the exported CSV column label."""

    if name in CSV_HEADER_NAME_OVERRIDES:
        return CSV_HEADER_NAME_OVERRIDES[name]
    words = []
    for token in name.split("_"):
        lowered = token.lower()
        if lowered in CSV_HEADER_TOKEN_OVERRIDES:
            words.append(CSV_HEADER_TOKEN_OVERRIDES[lowered])
            continue
        words.append(token.capitalize())
    return " ".join(words)
