from __future__ import annotations

import logging
from pathlib import Path

import openfoldpanel.pipeline as pipeline_module
from openfoldpanel.models import (
    AccessibilityEntry,
    ContactEntry,
    JobPanelData,
    MSAData,
    MSARow,
    ModelTracks,
    PipelineConfig,
    SecondaryStructureEntry,
    SequenceAxisPosition,
    dataclass_to_dict,
)
from openfoldpanel.pipeline import (
    _build_msa_data,
    _build_leading_display_overrides,
    _leading_gap_display_override,
    _project_alignment_to_query_axis,
    _select_display_msa_rows,
)
from openfoldpanel.render.layout import build_render_config
from openfoldpanel.render.svg_renderer import render_panel_svg


def test_project_alignment_preserves_raw_projected_leading_gap():
    rows = [
        MSARow(identifier="Query Sequence", sequence="--SWNP", is_query=True),
        MSARow(identifier="sp|P01674|KV3AM_MOUSE", sequence="MQ-WNP"),
    ]

    projected = _project_alignment_to_query_axis(rows)

    assert projected[0].sequence == "SWNP"
    assert projected[1].sequence == "-WNP"


def test_leading_display_override_uses_left_neighbor_even_when_later_gaps_exist():
    override = _leading_gap_display_override("-WN-P", "MQ-WN-P", [2, 3, 4, 5, 6])

    assert override == "Q"


def test_leading_display_override_uses_inserted_residue_between_first_and_second_query_columns():
    override = _leading_gap_display_override("-WNP", "-PWNP", [0, 2, 3, 4])

    assert override == "P"


def test_leading_display_override_is_none_when_no_left_residue_exists():
    override = _leading_gap_display_override("-WNP", "-WNP", [0, 1, 2, 3])

    assert override is None


def test_build_leading_display_overrides_tracks_rows_without_mutating_sequences():
    aligned_rows = [
        MSARow(identifier="Query Sequence", sequence="--SWNP", is_query=True),
        MSARow(identifier="sp|P01674|KV3AM_MOUSE", sequence="MQ-WN-P"),
    ]
    projected_rows = _project_alignment_to_query_axis(aligned_rows)

    overrides = _build_leading_display_overrides(aligned_rows, projected_rows)

    assert projected_rows[1].sequence == "-WN-"
    assert overrides == [None, "Q"]


def test_select_display_rows_keeps_first_homologs_in_original_order():
    rows = [
        MSARow(identifier="Query Sequence", sequence="ABCDEFGH", is_query=True),
        MSARow(identifier="h1", sequence="ABCD-FGH"),
        MSARow(identifier="h2", sequence="ABCDEFGH"),
        MSARow(identifier="h3", sequence="ABCXEFGH"),
        MSARow(identifier="h4", sequence="A-CDEFGH"),
    ]

    selected, filtered_count = _select_display_msa_rows(rows, max_homologs_displayed=2)

    assert [row.identifier for row in selected] == ["Query Sequence", "h1", "h2"]
    assert filtered_count == 0


def test_svg_uses_display_only_leading_override_without_mutating_msa_data():
    config = build_render_config(columns=4, font_size=12, max_homologs_displayed=1)
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", "ALA", residue, str(index + 1))
        for index, residue in enumerate("AAAA")
    ]
    model = ModelTracks(
        name="ranked_0_A",
        source_path="model.pdb",
        chain="A",
        secondary_structure=[SecondaryStructureEntry(index, "C", "coil") for index in range(4)],
        plddt=[90.0] * 4,
        accessibility=[AccessibilityEntry(index, None, 0.5, "accessible") for index in range(4)],
        contacts=[ContactEntry(index, None, None, None, None, None, None, None) for index in range(4)],
        display_name="ranked_0_A",
    )
    panel = JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=[model],
        msa=MSAData(
            enabled=True,
            query="AAAA",
            rows=[
                MSARow(identifier="Query Sequence", sequence="AAAA", is_query=True),
                MSARow(identifier="sp|Q15116|PDCD1_HUMAN", sequence="-AAA"),
            ],
            leading_display_overrides=[None, "P"],
        ),
        hydropathy=[],
        render_config=config,
    )

    svg, _ = render_panel_svg(panel)

    assert panel.msa.rows[1].sequence == "-AAA"
    assert ">P<" in svg


def test_dataclass_serialization_omits_display_only_leading_override_sidecar():
    payload = dataclass_to_dict(
        MSAData(
            enabled=True,
            query="AAAA",
            rows=[MSARow(identifier="Query Sequence", sequence="AAAA", is_query=True)],
            leading_display_overrides=[None],
        )
    )

    assert "leading_display_overrides" not in payload


def test_build_msa_data_uses_display_limit_as_search_limit(monkeypatch, tmp_path):
    captured: dict[str, object] = {}

    def fake_search_homologs(query_fasta, database, *, max_homologs_displayed, evalue, workdir, logger):
        captured["limit"] = max_homologs_displayed
        captured["evalue"] = evalue
        return [], []

    monkeypatch.setattr(pipeline_module, "search_homologs", fake_search_homologs)

    config = PipelineConfig(
        input_path=Path("demo.pdb"),
        outdir=tmp_path / "out",
        chain="ALL",
        columns=80,
        max_homologs_displayed=5,
        evalue="1e-6",
        font_size=12,
        hyd_window=3,
        msa_db=Path("/tmp/demo_db"),
        disable_msa=False,
        keep_temp=False,
        contact_cutoff=3.7,
        strong_contact_cutoff=3.2,
        verbose=False,
    )
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", "ALA", residue, str(index + 1))
        for index, residue in enumerate("ACDEFG")
    ]

    msa = _build_msa_data(
        axis=axis,
        config=config,
        workdir=tmp_path / "chain_a",
        logger=logging.getLogger("test"),
        warnings=[],
        job_name="demo",
        reference_chain_id="A",
    )

    assert captured["limit"] == 5
    assert captured["evalue"] == "1e-6"
    assert msa.enabled is False
    assert [row.identifier for row in msa.rows] == ["Query Sequence"]
