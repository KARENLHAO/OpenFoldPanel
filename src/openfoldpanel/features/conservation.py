"""Compute simple ESPript-style conservation annotations."""

from __future__ import annotations

from openfoldpanel.models import ConservationEntry, MSARow
from openfoldpanel.utils.residue_utils import compatible_similarity_group


def compute_conservation(rows: list[MSARow]) -> list[ConservationEntry]:
    """Compute identity/similarity fractions column-wise for an alignment."""

    if not rows:
        return []
    length = max(len(row.sequence) for row in rows)
    conservation: list[ConservationEntry] = []
    for index in range(length):
        residues = [row.sequence[index] if index < len(row.sequence) else "-" for row in rows]
        non_gap = [residue for residue in residues if residue != "-"]
        if not non_gap:
            conservation.append(
                ConservationEntry(
                    residue_index=index,
                    identity_fraction=0.0,
                    similarity_fraction=0.0,
                    style="default",
                )
            )
            continue
        unique = set(non_gap)
        max_count = max(non_gap.count(residue) for residue in unique)
        identity_fraction = max_count / len(non_gap)
        similarity_fraction = 1.0 if compatible_similarity_group(unique) else identity_fraction
        style = classify_conservation(identity_fraction, similarity_fraction)
        conservation.append(
            ConservationEntry(
                residue_index=index,
                identity_fraction=round(identity_fraction, 4),
                similarity_fraction=round(similarity_fraction, 4),
                style=style,
            )
        )
    return conservation


def classify_conservation(identity_fraction: float, similarity_fraction: float) -> str:
    """Convert conservation fractions into renderer styles."""

    if identity_fraction >= 1.0:
        return "identity"
    if identity_fraction >= 0.7 or similarity_fraction >= 0.7:
        return "similar"
    return "default"
