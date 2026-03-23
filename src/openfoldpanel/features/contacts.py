"""Residue-centric contact analysis for the reference chain."""

from __future__ import annotations

import math
from collections import defaultdict

from openfoldpanel.constants import DEFAULT_CONTACT_CUTOFF, DEFAULT_STRONG_CONTACT_CUTOFF
from openfoldpanel.models import ChainRecord, ContactEntry, ContactHit, ParsedStructure, ResidueRecord, SequenceAxisPosition
from openfoldpanel.utils.residue_utils import residue_entity_type


def compute_contacts(
    structure: ParsedStructure,
    reference_chain: ChainRecord,
    residue_by_axis_index: dict[int, ResidueRecord],
    axis: list[SequenceAxisPosition],
    *,
    cutoff: float = DEFAULT_CONTACT_CUTOFF,
    strong_cutoff: float = DEFAULT_STRONG_CONTACT_CUTOFF,
) -> list[ContactEntry]:
    """Compute residue-centric cross-chain and ligand contacts."""

    partner_residues: list[tuple[str, ResidueRecord]] = []
    for chain_id, chain in structure.chains.items():
        if chain_id == reference_chain.chain_id:
            continue
        for residue in chain.residues:
            partner_residues.append((chain_id, residue))

    contacts_by_axis: dict[int, list[ContactHit]] = defaultdict(list)

    for axis_index, residue in residue_by_axis_index.items():
        reference_atoms = [atom for atom in residue.atoms if atom.element.upper() != "H"]
        if not reference_atoms:
            continue
        hits: list[ContactHit] = []
        for partner_chain, partner_residue in partner_residues:
            partner_atoms = [atom for atom in partner_residue.atoms if atom.element.upper() != "H"]
            if not partner_atoms:
                continue
            min_distance = _minimum_atom_distance(reference_atoms, partner_atoms, cutoff)
            if min_distance is None:
                continue
            partner_type = _partner_type_for_residue(partner_residue)
            symbol = _contact_symbol(residue, partner_residue, partner_chain, partner_type)
            strength = "strong" if min_distance < strong_cutoff else "weak"
            hits.append(
                ContactHit(
                    partner_type=partner_type,
                    partner_chain=partner_chain,
                    partner_resname=partner_residue.resname,
                    partner_resid=partner_residue.residue_id.label,
                    min_distance=round(min_distance, 3),
                    symbol=symbol,
                    strength_category=strength,
                )
            )
        if hits:
            hits.sort(key=lambda hit: hit.min_distance)
            contacts_by_axis[axis_index] = hits

    track: list[ContactEntry] = []
    for position in axis:
        hits = contacts_by_axis.get(position.residue_index, [])
        if not hits:
            track.append(
                ContactEntry(
                    residue_index=position.residue_index,
                    partner_type=None,
                    partner_chain=None,
                    partner_resname=None,
                    partner_resid=None,
                    min_distance=None,
                    symbol=None,
                    strength_category=None,
                    is_multi_contact=False,
                    all_contacts=[],
                )
            )
            continue
        primary = hits[0]
        has_protein = any(hit.partner_type == "protein_chain" for hit in hits)
        has_non_protein = any(hit.partner_type != "protein_chain" for hit in hits)
        track.append(
            ContactEntry(
                residue_index=position.residue_index,
                partner_type=primary.partner_type,
                partner_chain=primary.partner_chain,
                partner_resname=primary.partner_resname,
                partner_resid=primary.partner_resid,
                min_distance=primary.min_distance,
                symbol=primary.symbol,
                strength_category=primary.strength_category,
                is_multi_contact=has_protein and has_non_protein,
                all_contacts=hits,
            )
        )
    return track


def _minimum_atom_distance(
    atoms_a: list,
    atoms_b: list,
    cutoff: float,
) -> float | None:
    min_distance = None
    cutoff_sq = cutoff * cutoff
    for atom_a in atoms_a:
        for atom_b in atoms_b:
            dx = atom_a.x - atom_b.x
            dy = atom_a.y - atom_b.y
            dz = atom_a.z - atom_b.z
            dist_sq = dx * dx + dy * dy + dz * dz
            if dist_sq <= cutoff_sq:
                distance = math.sqrt(dist_sq)
                if min_distance is None or distance < min_distance:
                    min_distance = distance
    return min_distance


def _partner_type_for_residue(residue: ResidueRecord) -> str:
    coarse = residue_entity_type(residue.resname, is_hetatm=all(atom.is_hetatm for atom in residue.atoms))
    if coarse == "protein":
        return "protein_chain"
    return coarse


def _contact_symbol(
    reference_residue: ResidueRecord,
    partner_residue: ResidueRecord,
    partner_chain: str,
    partner_type: str,
) -> str:
    if partner_type == "protein_chain":
        same_position = (
            reference_residue.resname == partner_residue.resname
            and reference_residue.residue_id.seq_id == partner_residue.residue_id.seq_id
        )
        return "#" if same_position else partner_chain
    if partner_type == "nucleic_acid":
        return "*"
    if partner_type == "porphyrin_like":
        return ":"
    if partner_type == "sugar":
        return '"'
    if partner_type == "ion":
        return "+"
    return "^"
