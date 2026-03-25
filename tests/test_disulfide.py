from __future__ import annotations

from pathlib import Path

from openfoldpanel.features.disulfide import infer_disulfides
from openfoldpanel.models import AtomRecord, ChainRecord, DisulfideBond, ParsedStructure, ResidueId, ResidueRecord


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


def _structure(**chains: list[ResidueRecord]) -> ParsedStructure:
    return ParsedStructure(
        name="demo",
        source_path=Path("demo.pdb"),
        chains={
            chain_id: ChainRecord(chain_id=chain_id, residues=residues, entity_type="protein")
            for chain_id, residues in chains.items()
        },
        format="pdb",
    )


def test_infer_disulfides_detects_chain_local_sg_pairs_once():
    chain_a = [
        _residue("A", 1, "CYS", sg_x=0.0),
        _residue("A", 2, "CYS", sg_x=2.1),
        _residue("A", 3, "CYS", sg_x=5.5),
    ]
    residue_by_axis_index = {index: residue for index, residue in enumerate(chain_a)}

    bonds = infer_disulfides(_structure(A=chain_a), residue_by_axis_index, current_chain_id="A")

    assert bonds == [
        DisulfideBond(
            residue_index_a=0,
            residue_index_b=1,
            chain_a="A",
            chain_b="A",
            bridge_scope="intramolecular",
        )
    ]


def test_infer_disulfides_detects_interchain_pairs_for_current_chain():
    chain_a = [_residue("A", 1, "CYS", sg_x=0.0)]
    chain_b = [_residue("B", 2, "CYS", sg_x=2.0)]
    residue_by_axis_index = {0: chain_a[0]}

    bonds = infer_disulfides(_structure(A=chain_a, B=chain_b), residue_by_axis_index, current_chain_id="A")

    assert bonds == [
        DisulfideBond(
            residue_index_a=0,
            residue_index_b=None,
            chain_a="A",
            chain_b="B",
            bridge_scope="intermolecular",
        )
    ]


def test_infer_disulfides_skips_pairs_missing_sg_and_non_cysteines():
    chain_a = [
        _residue("A", 1, "CYS", sg_x=0.0),
        _residue("A", 3, "CYS", sg_x=None),
        _residue("A", 4, "SER", sg_x=1.8),
    ]
    chain_b = [_residue("B", 2, "CYS", sg_x=None)]
    residue_by_axis_index = {index: residue for index, residue in enumerate(chain_a)}

    bonds = infer_disulfides(_structure(A=chain_a, B=chain_b), residue_by_axis_index, current_chain_id="A")

    assert bonds == []
