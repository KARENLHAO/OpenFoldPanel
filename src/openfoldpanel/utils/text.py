"""Helpers for generating readable display text."""

from __future__ import annotations

import re
from pathlib import Path


def humanize_identifier(value: str) -> str:
    """Convert underscore-heavy identifiers into readable title text."""

    collapsed = re.sub(r"[_\-]+", " ", value).strip()
    collapsed = re.sub(r"\s+", " ", collapsed)
    if not collapsed:
        return "Untitled"

    words: list[str] = []
    for token in collapsed.split(" "):
        if token.isupper() or any(char.isupper() for char in token[1:]) or any(char.isdigit() for char in token):
            words.append(token)
        else:
            words.append(token.capitalize())
    return " ".join(words)


def humanize_label(value: str) -> str:
    """Normalize labels for visible UI text without forcing title case."""

    collapsed = re.sub(r"[_\-]+", " ", value).strip()
    return re.sub(r"\s+", " ", collapsed) or "Untitled"


def humanize_chain_label(chain_id: str) -> str:
    """Return a readable chain label."""

    if chain_id == "_":
        return "未命名链"
    return f"链 {chain_id}"


def humanize_model_name(structure_name: str, chain_id: str) -> str:
    """Return a readable model label for reports and summaries."""

    return f"{humanize_identifier(structure_name)} / {humanize_chain_label(chain_id)}"


def safe_chain_slug(chain_id: str) -> str:
    """Convert a chain identifier into a filesystem-safe fragment."""

    if not chain_id or chain_id == "_":
        return "blank"
    sanitized = re.sub(r"[^A-Za-z0-9]+", "-", chain_id).strip("-")
    return sanitized or "blank"


def summarize_msa_database_path(value: str | Path | None) -> str:
    """Convert an MSA database path into a compact human-facing label."""

    if value is None:
        return "未设置"

    path = Path(value)
    leaf = path.name or str(path).rstrip("/").split("/")[-1]
    if not leaf:
        return "未设置"

    lowered_leaf = leaf.casefold()
    lowered_stem = Path(leaf).stem.casefold()
    aliases = {
        "pdbaa": "PDBAA",
        "swissprot": "SWISSPROT",
        "pdbaa50": "PDBAA50",
        "pdbaa70": "PDBAA70",
        "pdbaa90": "PDBAA90",
        "pdbaa95": "PDBAA95",
    }

    if lowered_leaf in aliases:
        return aliases[lowered_leaf]
    if lowered_stem in aliases:
        return aliases[lowered_stem]
    if lowered_leaf == "uniprot_sprot.fasta" or lowered_stem == "uniprot_sprot":
        return "SWISSPROT"
    return leaf
