from __future__ import annotations

from openfoldpanel.features.dssp_runner import DSSPResidueFeature
from openfoldpanel.features.secondary_structure import build_secondary_structure_track
from openfoldpanel.models import AtomRecord, ResidueId, ResidueRecord, SequenceAxisPosition


def _atom(atom_name: str, x: float, y: float, z: float) -> AtomRecord:
    return AtomRecord(atom_name=atom_name, element=atom_name[:1], x=x, y=y, z=z)


def _residue(seq_id: int, atom_positions: dict[str, tuple[float, float, float]]) -> ResidueRecord:
    return ResidueRecord(
        residue_id=ResidueId(chain_id="A", seq_id=seq_id),
        resname="ALA",
        atoms=[_atom(atom_name, *coords) for atom_name, coords in atom_positions.items()],
        residue_type="protein",
        one_letter="A",
        auth_seq_id=str(seq_id),
    )


def _axis(length: int) -> list[SequenceAxisPosition]:
    return [SequenceAxisPosition(index, "A", index + 1, "", "ALA", "A", str(index + 1)) for index in range(length)]


def _dssp_feature(residue: ResidueRecord, dssp_code: str = "T") -> DSSPResidueFeature:
    return DSSPResidueFeature(
        chain=residue.residue_id.chain_id,
        seq_id=residue.residue_id.seq_id,
        insertion_code=residue.residue_id.insertion_code,
        amino_acid=residue.one_letter,
        dssp_code=dssp_code,
        asa=42.0,
    )


def _build_track(
    residues: list[ResidueRecord],
    *,
    dssp_code: str = "T",
) -> list[str]:
    residue_by_axis_index = {index: residue for index, residue in enumerate(residues)}
    dssp_by_residue = {
        (
            residue.residue_id.chain_id,
            residue.residue_id.seq_id,
            residue.residue_id.insertion_code,
        ): _dssp_feature(residue, dssp_code)
        for residue in residues
    }
    track = build_secondary_structure_track(_axis(len(residues)), residue_by_axis_index, dssp_by_residue)
    return [entry.category for entry in track]


def _build_track_without_dssp(residues: list[ResidueRecord]) -> list[str]:
    residue_by_axis_index = {index: residue for index, residue in enumerate(residues)}
    track = build_secondary_structure_track(_axis(len(residues)), residue_by_axis_index, {})
    return [entry.category for entry in track]


def test_build_secondary_structure_track_maps_dssp_helix_codes_to_three_subtypes():
    residues = [
        _residue(1, {"N": (0.0, 0.0, 0.0), "CA": (1.0, 0.0, 0.0), "C": (2.0, 0.0, 0.0), "O": (3.0, 0.0, 0.0)}),
        _residue(2, {"N": (4.0, 0.0, 0.0), "CA": (5.0, 0.0, 0.0), "C": (6.0, 0.0, 0.0), "O": (7.0, 0.0, 0.0)}),
        _residue(3, {"N": (8.0, 0.0, 0.0), "CA": (9.0, 0.0, 0.0), "C": (10.0, 0.0, 0.0), "O": (11.0, 0.0, 0.0)}),
    ]
    residue_by_axis_index = {index: residue for index, residue in enumerate(residues)}
    dssp_by_residue = {
        ("A", 1, ""): _dssp_feature(residues[0], "G"),
        ("A", 2, ""): _dssp_feature(residues[1], "H"),
        ("A", 3, ""): _dssp_feature(residues[2], "I"),
    }

    track = build_secondary_structure_track(_axis(len(residues)), residue_by_axis_index, dssp_by_residue)

    assert [entry.category for entry in track] == ["three_ten_helix", "alpha_helix", "pi_helix"]


def test_build_secondary_structure_track_marks_three_ten_helix_from_backbone_geometry_without_dssp():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (1.0, 1.0, 0.0), "CA": (1.9, 0.7, 0.0), "C": (2.9, 1.0, 0.0), "O": (3.5, 1.0, 0.0)}),
        _residue(3, {"N": (3.0, 1.0, 0.0), "CA": (3.8, 0.2, 0.0), "C": (4.8, 1.0, 0.0), "O": (5.4, 1.0, 0.0)}),
        _residue(4, {"N": (0.0, 2.9, 0.0), "CA": (4.8, -1.0, 0.0), "C": (5.8, 0.0, 0.0), "O": (6.4, 0.0, 0.0)}),
    ]

    assert _build_track_without_dssp(residues) == [
        "three_ten_helix",
        "three_ten_helix",
        "three_ten_helix",
        "three_ten_helix",
    ]


def test_build_secondary_structure_track_marks_pi_helix_from_backbone_geometry_without_dssp():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (2.0, 1.0, 0.0), "CA": (2.4, -0.5, 0.0), "C": (3.0, 1.0, 0.0), "O": (3.6, 1.0, 0.0)}),
        _residue(3, {"N": (4.0, 1.0, 0.0), "CA": (4.4, -0.4, 0.0), "C": (5.0, 1.0, 0.0), "O": (5.6, 1.0, 0.0)}),
        _residue(4, {"N": (6.0, 1.0, 0.0), "CA": (6.4, -0.3, 0.0), "C": (7.0, 1.0, 0.0), "O": (7.6, 1.0, 0.0)}),
        _residue(5, {"N": (8.0, 1.0, 0.0), "CA": (8.4, -0.2, 0.0), "C": (9.0, 1.0, 0.0), "O": (9.6, 1.0, 0.0)}),
        _residue(6, {"N": (0.0, 3.0, 0.0), "CA": (10.6, -1.0, 0.0), "C": (11.6, 0.0, 0.0), "O": (12.2, 0.0, 0.0)}),
    ]

    assert _build_track_without_dssp(residues) == [
        "pi_helix",
        "pi_helix",
        "pi_helix",
        "pi_helix",
        "pi_helix",
        "pi_helix",
    ]


def test_build_secondary_structure_track_keeps_explicit_dssp_alpha_helix_when_geometry_could_suggest_other_subtypes():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (2.0, 1.0, 0.0), "CA": (2.4, -0.5, 0.0), "C": (3.0, 1.0, 0.0), "O": (3.6, 1.0, 0.0)}),
        _residue(3, {"N": (4.0, 1.0, 0.0), "CA": (4.4, -0.4, 0.0), "C": (5.0, 1.0, 0.0), "O": (5.6, 1.0, 0.0)}),
        _residue(4, {"N": (6.0, 1.0, 0.0), "CA": (6.4, -0.3, 0.0), "C": (7.0, 1.0, 0.0), "O": (7.6, 1.0, 0.0)}),
        _residue(5, {"N": (0.0, 3.0, 0.0), "CA": (8.4, -0.2, 0.0), "C": (9.0, 1.0, 0.0), "O": (9.6, 1.0, 0.0)}),
        _residue(6, {"N": (0.0, 3.2, 0.0), "CA": (10.6, -1.0, 0.0), "C": (11.6, 0.0, 0.0), "O": (12.2, 0.0, 0.0)}),
    ]

    assert _build_track(residues, dssp_code="H") == ["alpha_helix"] * 6


def test_build_secondary_structure_track_prefers_shorter_hydrogen_bond_distance_before_subtype_priority():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (2.0, 1.0, 0.0), "CA": (2.4, -0.5, 0.0), "C": (3.0, 1.0, 0.0), "O": (3.6, 1.0, 0.0)}),
        _residue(3, {"N": (4.0, 1.0, 0.0), "CA": (4.4, -0.4, 0.0), "C": (5.0, 1.0, 0.0), "O": (5.6, 1.0, 0.0)}),
        _residue(4, {"N": (0.0, 3.4, 0.0), "CA": (6.4, -0.3, 0.0), "C": (7.0, 1.0, 0.0), "O": (7.6, 1.0, 0.0)}),
        _residue(5, {"N": (0.0, 2.8, 0.0), "CA": (8.4, -0.2, 0.0), "C": (9.0, 1.0, 0.0), "O": (9.6, 1.0, 0.0)}),
    ]

    assert _build_track_without_dssp(residues) == [
        "alpha_helix",
        "alpha_helix",
        "alpha_helix",
        "alpha_helix",
        "alpha_helix",
    ]


def test_build_secondary_structure_track_marks_beta_turn_from_strict_backbone_geometry():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (1.0, 1.0, 0.0), "CA": (2.0, 1.0, 0.0), "C": (3.0, 1.0, 0.0), "O": (3.5, 1.0, 0.0)}),
        _residue(3, {"N": (3.0, 1.0, 0.0), "CA": (4.0, 1.0, 0.0), "C": (5.0, 1.0, 0.0), "O": (5.5, 1.0, 0.0)}),
        _residue(4, {"N": (0.0, 3.0, 0.0), "CA": (-3.0, 3.0, 0.0), "C": (-2.0, 3.0, 0.0), "O": (-1.5, 3.0, 0.0)}),
    ]

    assert _build_track(residues) == ["beta_turn", "beta_turn", "beta_turn", "beta_turn"]


def test_build_secondary_structure_track_marks_alpha_turn_from_strict_backbone_geometry():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (1.0, 1.0, 0.0), "CA": (2.0, 1.0, 0.0), "C": (3.0, 1.0, 0.0), "O": (3.5, 1.0, 0.0)}),
        _residue(3, {"N": (3.0, 1.0, 0.0), "CA": (4.0, 1.0, 0.0), "C": (5.0, 1.0, 0.0), "O": (5.5, 1.0, 0.0)}),
        _residue(4, {"N": (8.0, 8.0, 0.0), "CA": (8.0, 8.0, 0.0), "C": (9.0, 8.0, 0.0), "O": (9.5, 8.0, 0.0)}),
        _residue(5, {"N": (0.0, 3.0, 0.0), "CA": (-3.0, 3.0, 0.0), "C": (-2.0, 3.0, 0.0), "O": (-1.5, 3.0, 0.0)}),
    ]

    assert _build_track(residues) == ["alpha_turn", "alpha_turn", "alpha_turn", "alpha_turn", "alpha_turn"]


def test_build_secondary_structure_track_drops_turn_to_coil_when_hydrogen_bond_geometry_fails():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (1.0, 1.0, 0.0), "CA": (2.0, 1.0, 0.0), "C": (3.0, 1.0, 0.0), "O": (3.5, 1.0, 0.0)}),
        _residue(3, {"N": (3.0, 1.0, 0.0), "CA": (4.0, 1.0, 0.0), "C": (5.0, 1.0, 0.0), "O": (5.5, 1.0, 0.0)}),
        _residue(4, {"N": (0.0, 3.0, 0.0), "CA": (0.0, 2.0, 0.0), "C": (-2.0, 3.0, 0.0), "O": (-1.5, 3.0, 0.0)}),
    ]

    assert _build_track(residues) == ["coil", "coil", "coil", "coil"]


def test_build_secondary_structure_track_drops_turn_to_coil_when_backbone_atoms_are_missing():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (1.0, 1.0, 0.0), "CA": (2.0, 1.0, 0.0), "C": (3.0, 1.0, 0.0), "O": (3.5, 1.0, 0.0)}),
        _residue(3, {"N": (3.0, 1.0, 0.0), "CA": (4.0, 1.0, 0.0), "C": (5.0, 1.0, 0.0), "O": (5.5, 1.0, 0.0)}),
        _residue(4, {"CA": (-3.0, 3.0, 0.0), "C": (-2.0, 3.0, 0.0), "O": (-1.5, 3.0, 0.0)}),
    ]

    assert _build_track(residues) == ["coil", "coil", "coil", "coil"]


def test_build_secondary_structure_track_prefers_alpha_turn_when_alpha_and_beta_overlap():
    residues = [
        _residue(1, {"N": (-1.0, 0.0, 0.0), "CA": (0.0, -1.0, 0.0), "C": (1.0, 0.0, 0.0), "O": (0.0, 0.0, 0.0)}),
        _residue(2, {"N": (1.0, 1.0, 0.0), "CA": (2.0, 1.0, 0.0), "C": (3.0, 1.0, 0.0), "O": (3.5, 1.0, 0.0)}),
        _residue(3, {"N": (3.0, 1.0, 0.0), "CA": (4.0, 1.0, 0.0), "C": (5.0, 1.0, 0.0), "O": (5.5, 1.0, 0.0)}),
        _residue(4, {"N": (0.0, 3.0, 0.0), "CA": (-3.0, 3.0, 0.0), "C": (-2.0, 3.0, 0.0), "O": (-1.5, 3.0, 0.0)}),
        _residue(5, {"N": (0.0, 3.0, 0.0), "CA": (-2.0, 4.0, 0.0), "C": (-1.0, 4.0, 0.0), "O": (-0.5, 4.0, 0.0)}),
    ]

    assert _build_track(residues) == ["alpha_turn", "alpha_turn", "alpha_turn", "alpha_turn", "alpha_turn"]
