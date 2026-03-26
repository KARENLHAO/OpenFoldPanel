"""Secondary-structure normalization helpers."""

from __future__ import annotations

import math

from openfoldpanel.features.dssp_runner import DSSPResidueFeature
from openfoldpanel.models import AtomRecord, ResidueRecord, SecondaryStructureEntry, SequenceAxisPosition

ALPHA_TURN_CA_DISTANCE_MAX = 6.5
BETA_TURN_CA_DISTANCE_MAX = 7.0
BACKBONE_HBOND_DISTANCE_MAX = 3.5
BACKBONE_HBOND_ANGLE_MIN = 90.0
HELIX_SUBTYPE_PRIORITY = {
    "pi_helix": 0,
    "alpha_helix": 1,
    "three_ten_helix": 2,
}
HELIX_WINDOW_DEFINITIONS = (
    ("three_ten_helix", 4, 3),
    ("alpha_helix", 5, 4),
    ("pi_helix", 6, 5),
)


def normalize_dssp_code(code: str | None) -> str:
    """Map raw DSSP codes to the renderer's coarse categories."""

    if code is None:
        return "missing"
    upper = code.upper()
    if upper == "H":
        return "alpha_helix"
    if upper == "G":
        return "three_ten_helix"
    if upper == "I":
        return "pi_helix"
    if upper in {"E", "B"}:
        return "strand"
    return "coil"


def build_secondary_structure_track(
    axis: list[SequenceAxisPosition],
    residue_by_axis_index: dict[int, ResidueRecord],
    dssp_by_residue: dict[tuple[str, int, str], DSSPResidueFeature],
) -> list[SecondaryStructureEntry]:
    """Build a secondary-structure track aligned to the reference axis."""

    track: list[SecondaryStructureEntry] = []
    explicit_dssp_positions: set[int] = set()
    explicit_dssp_helix_positions: set[int] = set()
    for position in axis:
        residue = residue_by_axis_index.get(position.residue_index)
        feature = None
        if residue is not None:
            feature = dssp_by_residue.get(
                (
                    residue.residue_id.chain_id,
                    residue.residue_id.seq_id,
                    residue.residue_id.insertion_code,
                )
            )
        if feature:
            dssp_code = feature.dssp_code
            category = normalize_dssp_code(dssp_code)
            explicit_dssp_positions.add(position.residue_index)
            if dssp_code and dssp_code.upper() in {"G", "H", "I"}:
                explicit_dssp_helix_positions.add(position.residue_index)
        else:
            dssp_code = _approximate_secondary_code(position.residue_index, residue_by_axis_index)
            category = _approximate_secondary_category(dssp_code)
        track.append(
            SecondaryStructureEntry(
                residue_index=position.residue_index,
                dssp_code=dssp_code,
                category=category,
            )
        )
    _assign_helix_subtypes(track, residue_by_axis_index, explicit_dssp_positions, explicit_dssp_helix_positions)
    _assign_turn_subtypes(track, residue_by_axis_index)
    _finalize_helix_candidates(track)
    return track


def _approximate_secondary_category(code: str | None) -> str:
    if code == "H":
        return "helix_candidate"
    return normalize_dssp_code(code)


def _assign_helix_subtypes(
    track: list[SecondaryStructureEntry],
    residue_by_axis_index: dict[int, ResidueRecord],
    explicit_dssp_positions: set[int],
    explicit_dssp_helix_positions: set[int],
) -> None:
    """Refine helix-like residues into 3_10, alpha, or pi helix subtypes."""

    candidates: list[tuple[float, int, int, int, str]] = []
    for start in range(len(track)):
        for subtype, length, offset in HELIX_WINDOW_DEFINITIONS:
            distance = _qualifies_helix_window(
                track,
                residue_by_axis_index,
                explicit_dssp_positions,
                explicit_dssp_helix_positions,
                start=start,
                length=length,
                hbond_offset=offset,
            )
            if distance is None:
                continue
            candidates.append((distance, HELIX_SUBTYPE_PRIORITY[subtype], start, length, subtype))

    for residue_index in range(len(track)):
        entry = track[residue_index]
        if entry.category not in {"coil", "helix_candidate"}:
            continue
        if entry.residue_index in explicit_dssp_positions:
            continue
        covering = [
            candidate
            for candidate in candidates
            if candidate[2] <= residue_index < candidate[2] + candidate[3]
        ]
        if not covering:
            continue
        _, _, _, _, subtype = min(covering, key=lambda candidate: (candidate[0], candidate[1], candidate[2]))
        entry.category = subtype


def _assign_turn_subtypes(
    track: list[SecondaryStructureEntry],
    residue_by_axis_index: dict[int, ResidueRecord],
) -> None:
    """Promote eligible coil-like windows into alpha-turn or beta-turn segments."""

    for start in range(len(track) - 4):
        if _qualifies_turn_window(track, residue_by_axis_index, start=start, length=5, max_ca_distance=ALPHA_TURN_CA_DISTANCE_MAX):
            for index in range(start, start + 5):
                track[index].category = "alpha_turn"

    for start in range(len(track) - 3):
        if any(track[index].category == "alpha_turn" for index in range(start, start + 4)):
            continue
        if _qualifies_turn_window(track, residue_by_axis_index, start=start, length=4, max_ca_distance=BETA_TURN_CA_DISTANCE_MAX):
            for index in range(start, start + 4):
                track[index].category = "beta_turn"


def _qualifies_turn_window(
    track: list[SecondaryStructureEntry],
    residue_by_axis_index: dict[int, ResidueRecord],
    *,
    start: int,
    length: int,
    max_ca_distance: float,
) -> bool:
    window = track[start : start + length]
    if len(window) != length:
        return False
    if any(entry.category not in {"coil", "helix_candidate"} for entry in window):
        return False
    residues: list[ResidueRecord] = []
    for entry in window:
        residue = residue_by_axis_index.get(entry.residue_index)
        if residue is None:
            return False
        residues.append(residue)

    start_residue = residues[0]
    end_residue = residues[-1]
    start_ca = _atom_by_name(start_residue, "CA")
    start_c = _atom_by_name(start_residue, "C")
    start_o = _atom_by_name(start_residue, "O")
    end_n = _atom_by_name(end_residue, "N")
    end_ca = _atom_by_name(end_residue, "CA")
    if any(atom is None for atom in (start_ca, start_c, start_o, end_n, end_ca)):
        return False

    assert start_ca is not None and start_c is not None and start_o is not None and end_n is not None and end_ca is not None
    if _distance(start_ca, end_ca) >= max_ca_distance:
        return False
    if _distance(start_o, end_n) > BACKBONE_HBOND_DISTANCE_MAX:
        return False
    if _angle(end_n, start_o, start_c) < BACKBONE_HBOND_ANGLE_MIN:
        return False
    if _angle(start_o, end_n, end_ca) < BACKBONE_HBOND_ANGLE_MIN:
        return False
    return True


def _qualifies_helix_window(
    track: list[SecondaryStructureEntry],
    residue_by_axis_index: dict[int, ResidueRecord],
    explicit_dssp_positions: set[int],
    explicit_dssp_helix_positions: set[int],
    *,
    start: int,
    length: int,
    hbond_offset: int,
) -> float | None:
    window = track[start : start + length]
    if len(window) != length:
        return None
    if any(
        entry.residue_index in explicit_dssp_positions and entry.residue_index not in explicit_dssp_helix_positions
        for entry in window
    ):
        return None
    if any(entry.residue_index in explicit_dssp_helix_positions for entry in window):
        explicit_subtypes = {entry.category for entry in window if entry.residue_index in explicit_dssp_helix_positions}
        if len(explicit_subtypes) != 1:
            return None
        if subtype_for_offset(hbond_offset) not in explicit_subtypes:
            return None
    if any(entry.category in {"strand", "alpha_turn", "beta_turn", "missing"} for entry in window):
        return None

    residues: list[ResidueRecord] = []
    for entry in window:
        residue = residue_by_axis_index.get(entry.residue_index)
        if residue is None:
            return None
        if _atom_by_name(residue, "CA") is None:
            return None
        residues.append(residue)

    start_residue = residues[0]
    end_residue = residues[hbond_offset]
    start_o = _atom_by_name(start_residue, "O")
    end_n = _atom_by_name(end_residue, "N")
    if start_o is None or end_n is None:
        return None

    distance = _distance(start_o, end_n)
    if distance > BACKBONE_HBOND_DISTANCE_MAX:
        return None
    return distance


def subtype_for_offset(hbond_offset: int) -> str:
    return {
        3: "three_ten_helix",
        4: "alpha_helix",
        5: "pi_helix",
    }[hbond_offset]


def _finalize_helix_candidates(track: list[SecondaryStructureEntry]) -> None:
    for entry in track:
        if entry.category == "helix_candidate":
            entry.category = "alpha_helix"


def _atom_by_name(residue: ResidueRecord, atom_name: str) -> AtomRecord | None:
    return next((atom for atom in residue.atoms if atom.atom_name == atom_name), None)


def _distance(left: AtomRecord, right: AtomRecord) -> float:
    return math.sqrt((left.x - right.x) ** 2 + (left.y - right.y) ** 2 + (left.z - right.z) ** 2)


def _angle(left: AtomRecord, center: AtomRecord, right: AtomRecord) -> float:
    left_vector = (left.x - center.x, left.y - center.y, left.z - center.z)
    right_vector = (right.x - center.x, right.y - center.y, right.z - center.z)
    left_norm = math.sqrt(sum(component * component for component in left_vector))
    right_norm = math.sqrt(sum(component * component for component in right_vector))
    if left_norm == 0.0 or right_norm == 0.0:
        return 0.0
    cosine = sum(left_component * right_component for left_component, right_component in zip(left_vector, right_vector, strict=True))
    cosine /= left_norm * right_norm
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def _approximate_secondary_code(index: int, residue_by_axis_index: dict[int, ResidueRecord]) -> str | None:
    """Approximate a coarse secondary-structure code from CA geometry."""

    current = residue_by_axis_index.get(index)
    plus_two = residue_by_axis_index.get(index + 2)
    plus_three = residue_by_axis_index.get(index + 3)
    if current is None:
        return None
    current_ca = next((atom for atom in current.atoms if atom.atom_name == "CA"), None)
    if current_ca is None:
        return None
    if plus_three is not None:
        plus_three_ca = next((atom for atom in plus_three.atoms if atom.atom_name == "CA"), None)
        if plus_three_ca is not None:
            distance = ((current_ca.x - plus_three_ca.x) ** 2 + (current_ca.y - plus_three_ca.y) ** 2 + (current_ca.z - plus_three_ca.z) ** 2) ** 0.5
            if 4.8 <= distance <= 6.2:
                return "H"
    if plus_two is not None:
        plus_two_ca = next((atom for atom in plus_two.atoms if atom.atom_name == "CA"), None)
        if plus_two_ca is not None:
            distance = ((current_ca.x - plus_two_ca.x) ** 2 + (current_ca.y - plus_two_ca.y) ** 2 + (current_ca.z - plus_two_ca.z) ** 2) ** 0.5
            if distance >= 6.2:
                return "E"
    return "C"
