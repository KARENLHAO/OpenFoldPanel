"""Helpers for reading per-residue pLDDT-like confidence values."""

from __future__ import annotations

from openfoldpanel.models import ResidueRecord


def residue_plddt(residue: ResidueRecord) -> float | None:
    """Estimate per-residue pLDDT from heavy-atom B-factors when possible."""

    values = [
        atom.bfactor
        for atom in residue.atoms
        if atom.bfactor is not None and atom.element.upper() != "H"
    ]
    if not values:
        return None
    average = sum(values) / len(values)
    if 0.0 <= average <= 100.0:
        return round(average, 2)
    return None
