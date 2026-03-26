"""Normalize supported structure inputs into temporary PDB files."""

from __future__ import annotations

from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.PDB.Structure import Structure

PDB_RESSEQ_MIN = -999
PDB_RESSEQ_MAX = 9999


def normalize_structure_to_pdb(source_path: Path, output_dir: Path) -> Path:
    """Rewrite a supported structure file as a normalized temporary PDB."""

    structure = _load_structure(source_path)
    _validate_pdb_compatibility(structure, source_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{source_path.stem}.pdb"

    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(output_path))
    return output_path


def _load_structure(source_path: Path) -> Structure:
    suffix = source_path.suffix.lower()
    if suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    elif suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format for normalization: {source_path.name}")
    return parser.get_structure(source_path.stem, str(source_path))


def _validate_pdb_compatibility(structure: Structure, source_path: Path) -> None:
    for model in structure:
        for chain in model:
            chain_id = str(chain.id)
            if len(chain_id) > 1:
                raise ValueError(
                    f"{source_path.name} uses chain ID {chain_id!r}, which cannot be represented safely in PDB."
                )
            for residue in chain:
                _, seq_id, insertion_code = residue.id
                if seq_id < PDB_RESSEQ_MIN or seq_id > PDB_RESSEQ_MAX:
                    raise ValueError(
                        f"{source_path.name} uses residue number {seq_id}, which is outside the PDB range."
                    )
                normalized_insertion = str(insertion_code or "").strip()
                if len(normalized_insertion) > 1:
                    raise ValueError(
                        f"{source_path.name} uses insertion code {normalized_insertion!r}, which is outside the PDB limit."
                    )
