"""Disulfide-bond detection helpers."""

from __future__ import annotations

import math

from openfoldpanel.models import DisulfideBond, ResidueRecord


def infer_disulfides(residue_by_axis_index: dict[int, ResidueRecord]) -> list[DisulfideBond]:
    """Infer chain-local disulfide bonds from SG-SG distances within 2.2 Angstrom."""

    cysteines: list[tuple[int, ResidueRecord]] = [
        (index, residue)
        for index, residue in residue_by_axis_index.items()
        if residue.resname == "CYS"
    ]
    bonds: list[DisulfideBond] = []
    for left_idx, left_residue in cysteines:
        left_sg = next((atom for atom in left_residue.atoms if atom.atom_name == "SG"), None)
        if left_sg is None:
            continue
        for right_idx, right_residue in cysteines:
            if right_idx <= left_idx:
                continue
            if left_residue.residue_id.chain_id != right_residue.residue_id.chain_id:
                continue
            right_sg = next((atom for atom in right_residue.atoms if atom.atom_name == "SG"), None)
            if right_sg is None:
                continue
            distance = math.dist(
                (left_sg.x, left_sg.y, left_sg.z),
                (right_sg.x, right_sg.y, right_sg.z),
            )
            if distance <= 2.2:
                bonds.append(
                    DisulfideBond(
                        residue_index_a=left_idx,
                        residue_index_b=right_idx,
                        chain_a=left_residue.residue_id.chain_id,
                        chain_b=right_residue.residue_id.chain_id,
                    )
                )
    return bonds
