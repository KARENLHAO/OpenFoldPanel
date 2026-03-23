"""Kyte-Doolittle hydropathy calculations."""

from __future__ import annotations

from openfoldpanel.constants import HYDROPATHY_THRESHOLDS
from openfoldpanel.models import HydropathyEntry, SequenceAxisPosition

KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
    "X": 0.0,
}


def classify_hydropathy(value: float | None) -> str | None:
    """Bucket a hydropathy value into a named category."""

    if value is None:
        return None
    for category, (lower, upper) in HYDROPATHY_THRESHOLDS.items():
        if lower is None and value < upper:
            return category
        if upper is None and value > lower:
            return category
        if lower is not None and upper is not None and lower <= value <= upper:
            return category
    return None


def compute_hydropathy(axis: list[SequenceAxisPosition], window: int) -> list[HydropathyEntry]:
    """Compute a smoothed Kyte-Doolittle profile along the sequence axis."""

    if window <= 0:
        raise ValueError("hydropathy window must be positive")
    half_window = window // 2
    values = [KYTE_DOOLITTLE.get(position.one_letter, KYTE_DOOLITTLE["X"]) for position in axis]

    track: list[HydropathyEntry] = []
    for index, position in enumerate(axis):
        start = max(0, index - half_window)
        end = min(len(values), index + half_window + 1)
        averaged = round(sum(values[start:end]) / (end - start), 4) if start < end else None
        track.append(
            HydropathyEntry(
                residue_index=index,
                residue=position.one_letter,
                value=averaged,
                category=classify_hydropathy(averaged),
            )
        )
    return track
