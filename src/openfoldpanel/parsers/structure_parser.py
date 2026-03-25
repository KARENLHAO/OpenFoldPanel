"""Parsers for PDB and mmCIF structures using Gemmi when available."""

from __future__ import annotations

import logging
import shlex
from collections import defaultdict
from pathlib import Path

from openfoldpanel.models import AtomRecord, ChainRecord, ParsedStructure, ResidueId, ResidueRecord
from openfoldpanel.utils.residue_utils import residue_entity_type, three_to_one

try:
    import gemmi  # type: ignore
except ImportError:  # pragma: no cover - exercised in environments without gemmi
    gemmi = None


def parse_structure(path: Path, logger: logging.Logger) -> ParsedStructure:
    """Parse a PDB or mmCIF file into a normalized in-memory representation."""

    lower = path.name.lower()
    if gemmi is not None:
        try:
            return _parse_with_gemmi(path)
        except Exception as exc:  # pragma: no cover - best-effort fallback
            logger.warning("Gemmi failed for %s, falling back to pure Python parser: %s", path.name, exc)

    if lower.endswith(".pdb"):
        return _parse_pdb_text(path)
    if lower.endswith((".cif", ".mmcif")):
        return _parse_mmcif_text(path)
    raise ValueError(f"Unsupported structure format: {path}")


def _parse_with_gemmi(path: Path) -> ParsedStructure:
    structure = gemmi.read_structure(str(path))
    model = structure[0]
    chains: dict[str, ChainRecord] = {}
    for chain in model:
        residues: list[ResidueRecord] = []
        entity_type = "other_ligand"
        for residue in chain:
            if not residue:
                continue
            atoms: list[AtomRecord] = []
            resname = residue.name.strip().upper()
            het_flag = residue.het_flag.strip() not in {"", "A"}
            for atom in residue:
                element = atom.element.name.strip() or atom.name[:1].strip()
                atoms.append(
                    AtomRecord(
                        atom_name=atom.name.strip(),
                        element=element.upper(),
                        x=float(atom.pos.x),
                        y=float(atom.pos.y),
                        z=float(atom.pos.z),
                        bfactor=float(atom.b_iso) if atom.b_iso is not None else None,
                        occupancy=float(atom.occ) if atom.occ is not None else None,
                        is_hetatm=het_flag,
                    )
                )
            if not atoms:
                continue
            residue_id = ResidueId(
                chain_id=chain.name or "_",
                seq_id=int(residue.seqid.num),
                insertion_code=residue.seqid.icode.strip(),
            )
            residue_type = residue_entity_type(resname, is_hetatm=het_flag)
            if residue_type == "protein":
                entity_type = "protein"
            elif entity_type != "protein":
                entity_type = residue_type
            residues.append(
                ResidueRecord(
                    residue_id=residue_id,
                    resname=resname,
                    atoms=atoms,
                    residue_type=residue_type,
                    one_letter=three_to_one(resname),
                    auth_seq_id=f"{residue.seqid.num}{residue.seqid.icode.strip()}".strip(),
                )
            )
        if residues:
            chains[chain.name or "_"] = ChainRecord(chain_id=chain.name or "_", residues=residues, entity_type=entity_type)

    structure_format = "mmcif" if path.suffix.lower() in {".cif", ".mmcif"} else "pdb"
    return ParsedStructure(name=path.stem, source_path=path, chains=chains, format=structure_format)


def _parse_pdb_text(path: Path) -> ParsedStructure:
    chain_residues: dict[str, dict[tuple[int, str, str, bool], list[AtomRecord]]] = defaultdict(lambda: defaultdict(list))

    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        is_hetatm = line.startswith("HETATM")
        atom_name = line[12:16].strip()
        resname = line[17:20].strip().upper()
        chain_id = (line[21].strip() or "_")
        seq_id = int(line[22:26].strip())
        insertion_code = line[26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occupancy = float(line[54:60].strip()) if line[54:60].strip() else None
        bfactor = float(line[60:66].strip()) if line[60:66].strip() else None
        element = (line[76:78].strip() or atom_name[:1]).upper()
        atom = AtomRecord(
            atom_name=atom_name,
            element=element,
            x=x,
            y=y,
            z=z,
            occupancy=occupancy,
            bfactor=bfactor,
            is_hetatm=is_hetatm,
        )
        residue_key = (seq_id, insertion_code, resname, is_hetatm)
        chain_residues[chain_id][residue_key].append(atom)

    chains: dict[str, ChainRecord] = {}
    for chain_id, residue_map in chain_residues.items():
        residues: list[ResidueRecord] = []
        entity_type = "other_ligand"
        ordered_keys = sorted(residue_map, key=lambda key: (key[0], key[1], key[2], key[3]))
        for seq_id, insertion_code, resname, is_hetatm in ordered_keys:
            residue_type = residue_entity_type(resname, is_hetatm=is_hetatm)
            if residue_type == "protein":
                entity_type = "protein"
            elif entity_type != "protein":
                entity_type = residue_type
            residues.append(
                ResidueRecord(
                    residue_id=ResidueId(chain_id=chain_id, seq_id=seq_id, insertion_code=insertion_code),
                    resname=resname,
                    atoms=residue_map[(seq_id, insertion_code, resname, is_hetatm)],
                    residue_type=residue_type,
                    one_letter=three_to_one(resname),
                    auth_seq_id=f"{seq_id}{insertion_code}".strip(),
                )
            )
        chains[chain_id] = ChainRecord(chain_id=chain_id, residues=residues, entity_type=entity_type)
    return ParsedStructure(name=path.stem, source_path=path, chains=chains, format="pdb")


def _parse_mmcif_text(path: Path) -> ParsedStructure:
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_fields: list[str] = []
    atom_rows: list[list[str]] = []
    collecting_atom_loop = False

    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line or line.startswith("#"):
            index += 1
            continue
        if line == "loop_":
            collecting_atom_loop = False
            atom_fields = []
            atom_rows = []
            index += 1
            while index < len(lines) and lines[index].strip().startswith("_"):
                field_name = lines[index].strip()
                atom_fields.append(field_name)
                if field_name.startswith("_atom_site."):
                    collecting_atom_loop = True
                index += 1
            if collecting_atom_loop and atom_fields and all(field.startswith("_atom_site.") for field in atom_fields):
                while index < len(lines):
                    candidate = lines[index].strip()
                    if not candidate or candidate == "loop_" or candidate.startswith("_"):
                        break
                    atom_rows.append(shlex.split(lines[index], posix=True))
                    index += 1
                break
            continue
        index += 1

    if not atom_fields or not atom_rows:
        raise ValueError(f"Could not locate _atom_site loop in mmCIF file: {path}")

    field_positions = {field: idx for idx, field in enumerate(atom_fields)}

    def value(row: list[str], *names: str, default: str = "?") -> str:
        for name in names:
            idx = field_positions.get(name)
            if idx is not None and idx < len(row):
                val = row[idx]
                return "" if val in {".", "?"} else val
        return default

    chain_residues: dict[str, dict[tuple[int, str, str, bool], list[AtomRecord]]] = defaultdict(lambda: defaultdict(list))
    for row in atom_rows:
        group = value(row, "_atom_site.group_PDB", default="ATOM").upper()
        is_hetatm = group == "HETATM"
        atom_name = value(row, "_atom_site.auth_atom_id", "_atom_site.label_atom_id")
        resname = value(row, "_atom_site.auth_comp_id", "_atom_site.label_comp_id").upper()
        chain_id = value(row, "_atom_site.auth_asym_id", "_atom_site.label_asym_id", default="_") or "_"
        seq_raw = value(row, "_atom_site.auth_seq_id", "_atom_site.label_seq_id", default="0")
        if not seq_raw:
            continue
        seq_id = int(float(seq_raw))
        insertion_code = value(row, "_atom_site.pdbx_PDB_ins_code", default="")
        x = float(value(row, "_atom_site.Cartn_x", default="0.0"))
        y = float(value(row, "_atom_site.Cartn_y", default="0.0"))
        z = float(value(row, "_atom_site.Cartn_z", default="0.0"))
        occupancy_text = value(row, "_atom_site.occupancy", default="")
        bfactor_text = value(row, "_atom_site.B_iso_or_equiv", default="")
        element = value(row, "_atom_site.type_symbol", default=atom_name[:1]).upper()
        atom = AtomRecord(
            atom_name=atom_name,
            element=element,
            x=x,
            y=y,
            z=z,
            occupancy=float(occupancy_text) if occupancy_text else None,
            bfactor=float(bfactor_text) if bfactor_text else None,
            is_hetatm=is_hetatm,
        )
        residue_key = (seq_id, insertion_code, resname, is_hetatm)
        chain_residues[chain_id][residue_key].append(atom)

    chains: dict[str, ChainRecord] = {}
    for chain_id, residue_map in chain_residues.items():
        residues: list[ResidueRecord] = []
        entity_type = "other_ligand"
        for seq_id, insertion_code, resname, is_hetatm in sorted(residue_map, key=lambda key: (key[0], key[1], key[2], key[3])):
            residue_type = residue_entity_type(resname, is_hetatm=is_hetatm)
            if residue_type == "protein":
                entity_type = "protein"
            elif entity_type != "protein":
                entity_type = residue_type
            residues.append(
                ResidueRecord(
                    residue_id=ResidueId(chain_id=chain_id, seq_id=seq_id, insertion_code=insertion_code),
                    resname=resname,
                    atoms=residue_map[(seq_id, insertion_code, resname, is_hetatm)],
                    residue_type=residue_type,
                    one_letter=three_to_one(resname),
                    auth_seq_id=f"{seq_id}{insertion_code}".strip(),
                )
            )
        chains[chain_id] = ChainRecord(chain_id=chain_id, residues=residues, entity_type=entity_type)
    return ParsedStructure(name=path.stem, source_path=path, chains=chains, format="mmcif")


def select_reference_chain(structure: ParsedStructure, requested_chain: str) -> str:
    """Select the reference chain according to the project rules."""

    protein_chains = _protein_chain_ids(structure)
    if not protein_chains:
        raise ValueError(f"No protein chain found in structure: {structure.source_path}")

    normalized_request = _normalize_reference_request(requested_chain)
    if normalized_request != "AUTO":
        if normalized_request == "ALL":
            return protein_chains[0]
        chain = structure.chains.get(requested_chain)
        if chain is None:
            raise ValueError(f"Requested chain {requested_chain} not found in {structure.source_path.name}")
        if chain.entity_type != "protein":
            raise ValueError(f"Requested chain {requested_chain} is not a protein chain in {structure.source_path.name}")
        return requested_chain

    if "A" in protein_chains:
        return "A"
    return protein_chains[0]


def collect_reference_chains(structure: ParsedStructure, requested_chain: str) -> list[str]:
    """Collect all reference chains to render for a job."""

    protein_chains = _protein_chain_ids(structure)
    if not protein_chains:
        raise ValueError(f"No protein chain found in structure: {structure.source_path}")

    normalized_request = _normalize_reference_request(requested_chain)
    if normalized_request == "ALL":
        return protein_chains
    return [select_reference_chain(structure, requested_chain)]


def get_chain_or_best_match(structure: ParsedStructure, preferred_chain: str, reference_sequence: str) -> ChainRecord | None:
    """Find the preferred chain or the best protein-chain fallback by sequence similarity."""

    if preferred_chain in structure.chains and structure.chains[preferred_chain].entity_type == "protein":
        return structure.chains[preferred_chain]

    best_chain: ChainRecord | None = None
    best_score = float("-inf")
    for chain in structure.chains.values():
        if chain.entity_type != "protein":
            continue
        score = _sequence_identity_score(reference_sequence, chain.sequence)
        if score > best_score:
            best_score = score
            best_chain = chain
    return best_chain


def _sequence_identity_score(seq_a: str, seq_b: str) -> float:
    if not seq_a or not seq_b:
        return 0.0
    matches = sum(1 for a, b in zip(seq_a, seq_b) if a == b)
    return matches / max(len(seq_a), len(seq_b))


def _normalize_reference_request(requested_chain: str) -> str:
    if requested_chain.upper() in {"AUTO", "ALL"}:
        return requested_chain.upper()
    return requested_chain


def _protein_chain_ids(structure: ParsedStructure) -> list[str]:
    return [chain_id for chain_id, chain in structure.chains.items() if chain.entity_type == "protein"]
