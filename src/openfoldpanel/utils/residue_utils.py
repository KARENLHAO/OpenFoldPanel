"""Residue classification and sequence helpers."""

from __future__ import annotations

from openfoldpanel.constants import COMMON_IONS, COMMON_SUGARS, NUCLEIC_ACID_RESNAMES, PORPHYRIN_LIKE

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "ASX": "B",
    "GLX": "Z",
    "SEC": "U",
    "PYL": "O",
    "MSE": "M",
}

MAX_ASA_BY_RESIDUE = {
    "A": 121.0,
    "R": 265.0,
    "N": 187.0,
    "D": 187.0,
    "C": 148.0,
    "Q": 214.0,
    "E": 214.0,
    "G": 97.0,
    "H": 216.0,
    "I": 195.0,
    "L": 191.0,
    "K": 230.0,
    "M": 203.0,
    "F": 228.0,
    "P": 154.0,
    "S": 143.0,
    "T": 163.0,
    "W": 264.0,
    "Y": 255.0,
    "V": 165.0,
    "X": 200.0,
}


def three_to_one(resname: str) -> str:
    """Translate a residue name to a one-letter code."""

    return AA3_TO_1.get(resname.upper(), "X")


def is_protein_residue(resname: str) -> bool:
    """Return True if the residue is protein-like."""

    return resname.upper() in AA3_TO_1


def residue_entity_type(resname: str, *, is_hetatm: bool = False) -> str:
    """Classify a residue into a coarse entity type."""

    name = resname.upper()
    if is_protein_residue(name) and not is_hetatm:
        return "protein"
    if name in NUCLEIC_ACID_RESNAMES:
        return "nucleic_acid"
    if name in COMMON_SUGARS:
        return "sugar"
    if name in COMMON_IONS:
        return "ion"
    if name in PORPHYRIN_LIKE:
        return "porphyrin_like"
    return "other_ligand"


def compatible_similarity_group(chars: set[str]) -> bool:
    """Return True when residues belong to the same conservative group."""

    from openfoldpanel.constants import MSA_SIMILARITY_GROUPS

    if len(chars) <= 1:
        return True
    for group in MSA_SIMILARITY_GROUPS:
        if chars.issubset(group):
            return True
    return False
