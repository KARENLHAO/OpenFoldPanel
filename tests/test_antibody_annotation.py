from __future__ import annotations

import pytest

import openfoldpanel.features.antibody_annotation as antibody_module
from openfoldpanel.models import SequenceAxisPosition


def _axis(sequence: str) -> list[SequenceAxisPosition]:
    return [
        SequenceAxisPosition(index, "A", index + 1, "", "ALA", residue, str(index + 1))
        for index, residue in enumerate(sequence)
    ]


def _ordered_region(length: int) -> dict[str, str]:
    return {str(index): "A" for index in range(length)}


def test_annotate_antibody_chain_extracts_all_supported_schemes(monkeypatch):
    class FakeChain:
        def __init__(self, sequence, scheme):
            self.chain_type = "H"
            self.regions = {
                "FR1": _ordered_region(2),
                "CDR1": _ordered_region(3),
                "FR2": _ordered_region(2),
                "CDR2": _ordered_region(2),
                "FR3": _ordered_region(1),
                "CDR3": _ordered_region(2),
                "FR4": _ordered_region(1),
            }

    class FakeChainParseError(Exception):
        pass

    monkeypatch.setattr(antibody_module, "_load_abnumber_types", lambda: (FakeChain, FakeChainParseError))

    annotations, warnings = antibody_module.annotate_antibody_chain(
        "ABCDEFGHIJKLM",
        _axis("ABCDEFGHIJKLM"),
        chain_id="A",
    )

    assert warnings == []
    assert set(annotations) == {"kabat", "imgt", "chothia"}
    assert annotations["kabat"].chain_type == "heavy"
    assert [(region.name, region.start, region.end) for region in annotations["kabat"].regions] == [
        ("CDR1", 2, 5),
        ("CDR2", 7, 9),
        ("CDR3", 10, 12),
    ]
    assert annotations["imgt"].regions[0].display_label == "CDR1 - IMGT"


def test_annotate_antibody_chain_returns_empty_without_warning_for_non_antibody_sequences(monkeypatch):
    class FakeChainParseError(Exception):
        pass

    class FakeChain:
        def __init__(self, sequence, scheme):
            raise FakeChainParseError(f'Variable chain sequence not recognized: "{sequence}"')

    monkeypatch.setattr(antibody_module, "_load_abnumber_types", lambda: (FakeChain, FakeChainParseError))

    annotations, warnings = antibody_module.annotate_antibody_chain(
        "PEPTIDE",
        _axis("PEPTIDE"),
        chain_id="A",
    )

    assert annotations == {}
    assert warnings == []


def test_annotate_antibody_chain_warns_for_each_scheme_when_abnumber_is_unavailable(monkeypatch):
    monkeypatch.setattr(
        antibody_module,
        "_load_abnumber_types",
        lambda: (_ for _ in ()).throw(ModuleNotFoundError("No module named 'abnumber'")),
    )

    annotations, warnings = antibody_module.annotate_antibody_chain(
        "PEPTIDE",
        _axis("PEPTIDE"),
        chain_id="B",
    )

    assert annotations == {}
    assert len(warnings) == 3
    assert any("reference chain B" in warning for warning in warnings)
    assert any("Kabat" in warning for warning in warnings)
    assert any("IMGT" in warning for warning in warnings)
    assert any("Chothia" in warning for warning in warnings)


def test_annotate_antibody_chain_keeps_successful_schemes_when_one_fails(monkeypatch):
    class FakeChainParseError(Exception):
        pass

    class FakeChain:
        def __init__(self, sequence, scheme):
            if scheme == "imgt":
                raise FakeChainParseError("Failed to number IMGT")
            self.chain_type = "L"
            self.regions = {
                "FR1": _ordered_region(1),
                "CDR1": _ordered_region(1),
                "FR2": _ordered_region(1),
                "CDR2": _ordered_region(1),
                "FR3": _ordered_region(1),
                "CDR3": _ordered_region(1),
                "FR4": _ordered_region(1),
            }

    monkeypatch.setattr(antibody_module, "_load_abnumber_types", lambda: (FakeChain, FakeChainParseError))

    annotations, warnings = antibody_module.annotate_antibody_chain(
        "ABCDEFG",
        _axis("ABCDEFG"),
        chain_id="A",
    )

    assert set(annotations) == {"kabat", "chothia"}
    assert len(warnings) == 1
    assert "IMGT" in warnings[0]


def test_normalize_antibody_scheme_rejects_invalid_values():
    with pytest.raises(ValueError):
        antibody_module.normalize_antibody_scheme("aho")
