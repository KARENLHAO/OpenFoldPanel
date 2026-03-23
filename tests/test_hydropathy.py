from __future__ import annotations

from openfoldpanel.features.hydropathy import compute_hydropathy
from openfoldpanel.models import SequenceAxisPosition


def test_hydropathy_window_and_categories():
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", resname, aa, str(index + 1))
        for index, (resname, aa) in enumerate([("ILE", "I"), ("VAL", "V"), ("GLU", "E"), ("LYS", "K"), ("LEU", "L")])
    ]
    track = compute_hydropathy(axis, window=3)

    assert len(track) == 5
    assert track[0].category == "hydrophobic"
    assert track[2].category == "intermediate"
    assert track[3].category == "intermediate"
    assert track[0].value == round((4.5 + 4.2) / 2.0, 4)
