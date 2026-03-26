from __future__ import annotations

import pytest

from openfoldpanel.parsers.structure_preprocessor import normalize_structure_to_pdb
from tests.conftest import write_test_mmcif, write_test_pdb


def test_normalize_structure_to_pdb_rewrites_mmcif_input(tmp_path):
    source = write_test_mmcif(tmp_path / "ranked_2.cif")
    output_dir = tmp_path / "normalized"

    normalized = normalize_structure_to_pdb(source, output_dir)

    assert normalized.name == "ranked_2.pdb"
    assert normalized.exists()
    assert normalized.read_text().startswith("ATOM")


def test_normalize_structure_to_pdb_rewrites_existing_pdb_input(tmp_path):
    source = write_test_pdb(tmp_path / "ranked_1.pdb")
    output_dir = tmp_path / "normalized"

    normalized = normalize_structure_to_pdb(source, output_dir)

    assert normalized.name == "ranked_1.pdb"
    assert normalized.exists()
    assert normalized.read_text().startswith("ATOM")


def test_normalize_structure_to_pdb_rejects_multi_character_chain_ids(tmp_path):
    source = write_test_mmcif(tmp_path / "bad_chain.cif", chain_id="AA")

    with pytest.raises(ValueError) as exc_info:
        normalize_structure_to_pdb(source, tmp_path / "normalized")

    assert "chain ID" in str(exc_info.value)
    assert "AA" in str(exc_info.value)


def test_normalize_structure_to_pdb_rejects_out_of_range_residue_numbers(tmp_path):
    source = write_test_mmcif(tmp_path / "bad_resseq.cif", seq_id_start=10000)

    with pytest.raises(ValueError) as exc_info:
        normalize_structure_to_pdb(source, tmp_path / "normalized")

    assert "residue number" in str(exc_info.value)
    assert "10000" in str(exc_info.value)
