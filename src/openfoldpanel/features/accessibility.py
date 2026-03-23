"""Accessibility calculations derived from DSSP output."""

from __future__ import annotations

from openfoldpanel.constants import ACCESSIBILITY_THRESHOLDS
from openfoldpanel.features.dssp_runner import DSSPResidueFeature
from openfoldpanel.models import AccessibilityEntry, ResidueRecord, SequenceAxisPosition
from openfoldpanel.utils.residue_utils import MAX_ASA_BY_RESIDUE


def classify_relative_accessibility(value: float | None) -> str | None:
    """Bucket a relative accessibility value into a named category."""

    if value is None:
        return None
    for category, (lower, upper) in ACCESSIBILITY_THRESHOLDS.items():
        if lower is None and value < upper:
            return category
        if upper is None and value > lower:
            return category
        if lower is not None and upper is not None and lower <= value <= upper:
            return category
    return None


def build_accessibility_track(
    axis: list[SequenceAxisPosition],
    residue_by_axis_index: dict[int, ResidueRecord],
    dssp_by_residue: dict[tuple[str, int, str], DSSPResidueFeature],
) -> list[AccessibilityEntry]:
    """Build relative and absolute accessibility values along the reference axis."""

    approximate_values = _approximate_relative_accessibility(residue_by_axis_index)
    track: list[AccessibilityEntry] = []
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
        absolute = feature.asa if feature else None
        if absolute is None:
            relative = approximate_values.get(position.residue_index)
        else:
            max_asa = MAX_ASA_BY_RESIDUE.get(position.one_letter, MAX_ASA_BY_RESIDUE["X"])
            relative = round(absolute / max_asa, 4) if max_asa else None
        track.append(
            AccessibilityEntry(
                residue_index=position.residue_index,
                absolute=absolute,
                relative=relative,
                category=classify_relative_accessibility(relative),
            )
        )
    return track


def _approximate_relative_accessibility(
    residue_by_axis_index: dict[int, ResidueRecord],
) -> dict[int, float]:
    """Approximate relative accessibility from local heavy-atom crowding."""

    all_atoms: list[tuple[int, float, float, float]] = []
    for index, residue in residue_by_axis_index.items():
        for atom in residue.atoms:
            if atom.element.upper() == "H":
                continue
            all_atoms.append((index, atom.x, atom.y, atom.z))

    approximations: dict[int, float] = {}
    for index, residue in residue_by_axis_index.items():
        count = 0
        for atom in residue.atoms:
            if atom.element.upper() == "H":
                continue
            for other_index, x, y, z in all_atoms:
                if other_index == index:
                    continue
                dx = atom.x - x
                dy = atom.y - y
                dz = atom.z - z
                if dx * dx + dy * dy + dz * dz <= 100.0:
                    count += 1
        relative = max(0.0, min(1.2, 1.2 - (count / 60.0)))
        approximations[index] = round(relative, 4)
    return approximations
