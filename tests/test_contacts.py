from __future__ import annotations

from openfoldpanel.features.contacts import compute_contacts
from openfoldpanel.parsers.sequence_mapper import build_sequence_axis
from openfoldpanel.parsers.structure_parser import parse_structure
from tests.conftest import write_test_pdb


def test_contacts_capture_primary_and_multi_contact(tmp_path):
    path = write_test_pdb(tmp_path / "contacts.pdb")
    structure = parse_structure(path, __import__("logging").getLogger("test"))
    reference_chain = structure.chains["A"]
    axis = build_sequence_axis(reference_chain)
    residue_by_axis_index = {index: residue for index, residue in enumerate(reference_chain.residues)}

    contacts = compute_contacts(structure, reference_chain, residue_by_axis_index, axis, cutoff=3.7, strong_cutoff=3.2)

    assert contacts[0].symbol == "+"
    assert contacts[0].strength_category == "strong"
    assert contacts[0].is_multi_contact is True
    assert any(hit.partner_type == "protein_chain" for hit in contacts[0].all_contacts)
    assert any(hit.partner_type == "ion" for hit in contacts[0].all_contacts)
