from __future__ import annotations

from openfoldpanel.features.disulfide import infer_disulfides
from openfoldpanel.models import AtomRecord, DisulfideBond, ResidueId, ResidueRecord


def _residue(
    chain_id: str,
    seq_id: int,
    resname: str,
    *,
    sg_x: float | None,
) -> ResidueRecord:
    atoms = []
    if sg_x is not None:
        atoms.append(AtomRecord(atom_name="SG", element="S", x=sg_x, y=0.0, z=0.0))
    return ResidueRecord(
        residue_id=ResidueId(chain_id=chain_id, seq_id=seq_id),
        resname=resname,
        atoms=atoms,
        residue_type="protein",
        one_letter="C" if resname == "CYS" else "A",
        auth_seq_id=str(seq_id),
    )


def test_infer_disulfides_detects_chain_local_sg_pairs_once():
    residue_by_axis_index = {
        0: _residue("A", 1, "CYS", sg_x=0.0),
        1: _residue("A", 2, "CYS", sg_x=2.1),
        2: _residue("A", 3, "CYS", sg_x=5.5),
    }

    bonds = infer_disulfides(residue_by_axis_index)

    assert bonds == [DisulfideBond(residue_index_a=0, residue_index_b=1, chain_a="A", chain_b="A")]


def test_infer_disulfides_skips_cross_chain_pairs_missing_sg_and_non_cysteines():
    residue_by_axis_index = {
        0: _residue("A", 1, "CYS", sg_x=0.0),
        1: _residue("B", 2, "CYS", sg_x=2.0),
        2: _residue("A", 3, "CYS", sg_x=None),
        3: _residue("A", 4, "SER", sg_x=1.8),
    }

    bonds = infer_disulfides(residue_by_axis_index)

    assert bonds == []
