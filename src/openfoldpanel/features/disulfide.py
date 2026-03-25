"""Disulfide-bond detection helpers."""

from __future__ import annotations

import math
from dataclasses import dataclass

from openfoldpanel.models import AtomRecord, DisulfideBond, ParsedStructure, ResidueRecord


DISULFIDE_DISTANCE_CUTOFF = 2.2


@dataclass(frozen=True, slots=True)
class _CandidateDisulfide:
    distance: float
    residue_key_a: tuple[str, int, str]
    residue_key_b: tuple[str, int, str]
    residue_index_a: int
    residue_index_b: int | None
    chain_b: str
    bridge_scope: str


def infer_disulfides(
    structure: ParsedStructure,
    residue_by_axis_index: dict[int, ResidueRecord],
    *,
    current_chain_id: str,
) -> list[DisulfideBond]:
    """Infer local intramolecular and interchain disulfides for one rendered chain."""

    local_cysteines: list[tuple[int, ResidueRecord, AtomRecord]] = [
        (index, residue, sg_atom)
        for index, residue in residue_by_axis_index.items()
        if residue.resname == "CYS"
        and residue.residue_id.chain_id == current_chain_id
        and (sg_atom := _sg_atom(residue)) is not None
    ]
    if not local_cysteines:
        return []

    candidates: list[_CandidateDisulfide] = []
    for left_pos, (left_idx, left_residue, left_sg) in enumerate(local_cysteines):
        for right_idx, right_residue, right_sg in local_cysteines[left_pos + 1 :]:
            distance = _sg_distance(left_sg, right_sg)
            if distance > DISULFIDE_DISTANCE_CUTOFF:
                continue
            candidates.append(
                _CandidateDisulfide(
                    distance=distance,
                    residue_key_a=_residue_key(left_residue),
                    residue_key_b=_residue_key(right_residue),
                    residue_index_a=left_idx,
                    residue_index_b=right_idx,
                    chain_b=current_chain_id,
                    bridge_scope="intramolecular",
                )
            )

        for partner_chain_id, partner_residue, partner_sg in _partner_chain_cysteines(structure, current_chain_id):
            distance = _sg_distance(left_sg, partner_sg)
            if distance > DISULFIDE_DISTANCE_CUTOFF:
                continue
            candidates.append(
                _CandidateDisulfide(
                    distance=distance,
                    residue_key_a=_residue_key(left_residue),
                    residue_key_b=_residue_key(partner_residue),
                    residue_index_a=left_idx,
                    residue_index_b=None,
                    chain_b=partner_chain_id,
                    bridge_scope="intermolecular",
                )
            )

    bonds: list[DisulfideBond] = []
    used_residues: set[tuple[str, int, str]] = set()
    for candidate in sorted(candidates, key=lambda item: item.distance):
        if candidate.residue_key_a in used_residues or candidate.residue_key_b in used_residues:
            continue
        used_residues.add(candidate.residue_key_a)
        used_residues.add(candidate.residue_key_b)
        bonds.append(
            DisulfideBond(
                residue_index_a=candidate.residue_index_a,
                residue_index_b=candidate.residue_index_b,
                chain_a=current_chain_id,
                chain_b=candidate.chain_b,
                bridge_scope=candidate.bridge_scope,
            )
        )
    return bonds


def _partner_chain_cysteines(
    structure: ParsedStructure,
    current_chain_id: str,
) -> list[tuple[str, ResidueRecord, AtomRecord]]:
    cysteines: list[tuple[str, ResidueRecord, AtomRecord]] = []
    for chain_id, chain in structure.chains.items():
        if chain_id == current_chain_id or chain.entity_type != "protein":
            continue
        for residue in chain.residues:
            sg_atom = _sg_atom(residue)
            if residue.resname == "CYS" and sg_atom is not None:
                cysteines.append((chain_id, residue, sg_atom))
    return cysteines


def _sg_atom(residue: ResidueRecord) -> AtomRecord | None:
    return next((atom for atom in residue.atoms if atom.atom_name == "SG"), None)


def _sg_distance(left_atom: AtomRecord, right_atom: AtomRecord) -> float:
    return math.dist((left_atom.x, left_atom.y, left_atom.z), (right_atom.x, right_atom.y, right_atom.z))


def _residue_key(residue: ResidueRecord) -> tuple[str, int, str]:
    residue_id = residue.residue_id
    return residue_id.chain_id, residue_id.seq_id, residue_id.insertion_code
