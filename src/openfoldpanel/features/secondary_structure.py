"""Secondary-structure normalization helpers."""

from __future__ import annotations

from openfoldpanel.features.dssp_runner import DSSPResidueFeature
from openfoldpanel.models import ResidueRecord, SecondaryStructureEntry, SequenceAxisPosition


def normalize_dssp_code(code: str | None) -> str:
    """Map raw DSSP codes to the renderer's coarse categories."""

    if code is None:
        return "missing"
    upper = code.upper()
    if upper in {"H", "G", "I"}:
        return "helix"
    if upper in {"E", "B"}:
        return "strand"
    if upper == "T":
        return "turn"
    return "coil"


def build_secondary_structure_track(
    axis: list[SequenceAxisPosition],
    residue_by_axis_index: dict[int, ResidueRecord],
    dssp_by_residue: dict[tuple[str, int, str], DSSPResidueFeature],
) -> list[SecondaryStructureEntry]:
    """Build a secondary-structure track aligned to the reference axis."""

    track: list[SecondaryStructureEntry] = []
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
        dssp_code = feature.dssp_code if feature else _approximate_secondary_code(position.residue_index, residue_by_axis_index)
        track.append(
            SecondaryStructureEntry(
                residue_index=position.residue_index,
                dssp_code=dssp_code,
                category=normalize_dssp_code(dssp_code),
            )
        )
    return track


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
