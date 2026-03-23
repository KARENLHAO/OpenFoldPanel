from __future__ import annotations

from openfoldpanel.features.hydropathy import compute_hydropathy
from openfoldpanel.models import (
    AccessibilityEntry,
    ContactEntry,
    JobPanelData,
    MSAData,
    MSARow,
    ModelTracks,
    SecondaryStructureEntry,
    SequenceAxisPosition,
)
from openfoldpanel.render.layout import build_panel_layout, build_render_config
from openfoldpanel.render.svg_renderer import render_panel_svg


def _build_tick_panel(seq_ids: list[int]) -> JobPanelData:
    config = build_render_config(columns=len(seq_ids), font_size=12)
    axis = [
        SequenceAxisPosition(index, "A", seq_id, "", "ALA", "A", str(seq_id))
        for index, seq_id in enumerate(seq_ids)
    ]
    hydropathy = compute_hydropathy(axis, window=3)
    model = ModelTracks(
        name="ranked_0_A",
        source_path="model.pdb",
        chain="A",
        secondary_structure=[SecondaryStructureEntry(index, "C", "coil") for index in range(len(seq_ids))],
        plddt=[90.0] * len(seq_ids),
        accessibility=[AccessibilityEntry(index, None, 0.5, "accessible") for index in range(len(seq_ids))],
        contacts=[ContactEntry(index, None, None, None, None, None, None, None) for index in range(len(seq_ids))],
        display_name="Tick Demo / 链 A",
    )
    return JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=[model],
        msa=MSAData(enabled=False, query="A" * len(seq_ids), rows=[MSARow(identifier="query_sequence", sequence="A" * len(seq_ids), is_query=True)]),
        hydropathy=hydropathy,
        render_config=config,
    )


def test_layout_wraps_blocks_and_preserves_label_width():
    base_config = build_render_config(columns=5, font_size=12)
    assert base_config.cell_width > 12 * 0.92
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", "ALA", "A", str(index + 1))
        for index in range(12)
    ]
    hydropathy = compute_hydropathy(axis, window=3)
    model = ModelTracks(
        name="ranked_0_A",
        source_path="model.pdb",
        chain="A",
        secondary_structure=[SecondaryStructureEntry(index, "C", "coil") for index in range(12)],
        plddt=[90.0] * 12,
        accessibility=[AccessibilityEntry(index, None, 0.5, "accessible") for index in range(12)],
        contacts=[ContactEntry(index, None, None, None, None, None, None, None) for index in range(12)],
        display_name="Very Long Demonstration Model - Chain A",
    )
    panel = JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=[model],
        msa=MSAData(enabled=False, query="A" * 12, rows=[MSARow(identifier="query_sequence", sequence="A" * 12, is_query=True)]),
        hydropathy=hydropathy,
        render_config=base_config,
    )

    layout = build_panel_layout(panel)
    assert len(layout.blocks) == 3
    assert layout.annotation_height > 0
    assert layout.render_config.label_width >= base_config.label_width
    assert layout.blocks[0].end - layout.blocks[0].start == 5
    assert any(row.kind == "secondary" and row.label == "ranked_0_A" for row in layout.rows)
    assert any(row.kind == "confidence" and row.label == "置信度" for row in layout.rows)
    assert any(row.kind == "msa_query" and row.label == "查询序列" for row in layout.rows)
    secondary_index = next(index for index, row in enumerate(layout.rows) if row.kind == "secondary")
    assert layout.rows[secondary_index + 1].kind == "confidence"

    first_block = layout.blocks[0]
    grid_end = (
        first_block.x
        + layout.render_config.margin
        + layout.render_config.label_width
        + (first_block.end - first_block.start) * layout.render_config.cell_width
    )
    assert layout.width - grid_end >= layout.render_config.font_size


def test_turn_segments_render_as_curves_without_text_overlap():
    config = build_render_config(columns=8, font_size=12)
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", "ALA", "A", str(index + 1))
        for index in range(8)
    ]
    hydropathy = compute_hydropathy(axis, window=3)
    secondary = [
        SecondaryStructureEntry(0, "C", "coil"),
        SecondaryStructureEntry(1, "T", "turn"),
        SecondaryStructureEntry(2, "T", "turn"),
        SecondaryStructureEntry(3, "T", "turn"),
        SecondaryStructureEntry(4, "C", "coil"),
        SecondaryStructureEntry(5, "C", "coil"),
        SecondaryStructureEntry(6, "C", "coil"),
        SecondaryStructureEntry(7, "C", "coil"),
    ]
    model = ModelTracks(
        name="ranked_0_A",
        source_path="model.pdb",
        chain="A",
        secondary_structure=secondary,
        plddt=[90.0] * 8,
        accessibility=[AccessibilityEntry(index, None, 0.5, "accessible") for index in range(8)],
        contacts=[ContactEntry(index, None, None, None, None, None, None, None) for index in range(8)],
        display_name="Turn Demo / 链 A",
    )
    panel = JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=[model],
        msa=MSAData(enabled=False, query="A" * 8, rows=[MSARow(identifier="query_sequence", sequence="A" * 8, is_query=True)]),
        hydropathy=hydropathy,
        render_config=config,
    )

    svg, _ = render_panel_svg(panel)
    assert 'class="turn-track"' in svg
    assert ">TT<" not in svg


def test_layout_limits_displayed_homolog_rows_without_extra_secondary_annotation_band():
    config = build_render_config(columns=8, font_size=12, msa_display_rows=1)
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", "ALA", "A", str(index + 1))
        for index in range(8)
    ]
    hydropathy = compute_hydropathy(axis, window=3)
    secondary = [
        SecondaryStructureEntry(0, "E", "strand"),
        SecondaryStructureEntry(1, "E", "strand"),
        SecondaryStructureEntry(2, "C", "coil"),
        SecondaryStructureEntry(3, "E", "strand"),
        SecondaryStructureEntry(4, "E", "strand"),
        SecondaryStructureEntry(5, "C", "coil"),
        SecondaryStructureEntry(6, "E", "strand"),
        SecondaryStructureEntry(7, "E", "strand"),
    ]
    model = ModelTracks(
        name="ranked_0_A",
        source_path="model.pdb",
        chain="A",
        secondary_structure=secondary,
        plddt=[90.0] * 8,
        accessibility=[AccessibilityEntry(index, None, 0.5, "accessible") for index in range(8)],
        contacts=[ContactEntry(index, None, None, None, None, None, None, None) for index in range(8)],
        display_name="ranked_0_A",
    )
    panel = JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=[model],
        msa=MSAData(
            enabled=True,
            query="ACDEFGHI",
            rows=[
                MSARow(identifier="Query Sequence", sequence="ACDEFGHI", is_query=True),
                MSARow(identifier="Hit 1", sequence="ACDEFGHI"),
                MSARow(identifier="Hit 2", sequence="ACDEFGHI"),
            ],
        ),
        hydropathy=hydropathy,
        render_config=config,
    )

    layout = build_panel_layout(panel)
    homolog_rows = [row for row in layout.rows if row.kind == "msa_homolog"]
    assert len(homolog_rows) == 1

    svg, _ = render_panel_svg(panel)
    assert ">β1<" in svg
    assert ">β2<" in svg
    assert ">β3<" in svg


def test_structure_annotations_render_only_once_for_top_secondary_track():
    config = build_render_config(columns=8, font_size=12)
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", "ALA", "A", str(index + 1))
        for index in range(8)
    ]
    hydropathy = compute_hydropathy(axis, window=3)
    secondary = [
        SecondaryStructureEntry(0, "H", "helix"),
        SecondaryStructureEntry(1, "H", "helix"),
        SecondaryStructureEntry(2, "E", "strand"),
        SecondaryStructureEntry(3, "E", "strand"),
        SecondaryStructureEntry(4, "C", "coil"),
        SecondaryStructureEntry(5, "E", "strand"),
        SecondaryStructureEntry(6, "E", "strand"),
        SecondaryStructureEntry(7, "C", "coil"),
    ]
    models = [
        ModelTracks(
            name=f"ranked_{index}_A",
            source_path=f"model_{index}.pdb",
            chain="A",
            secondary_structure=secondary,
            plddt=[90.0] * 8,
            accessibility=[AccessibilityEntry(position, None, 0.5, "accessible") for position in range(8)],
            contacts=[ContactEntry(position, None, None, None, None, None, None, None) for position in range(8)],
            display_name=f"Top Annotation Demo {index}",
        )
        for index in range(3)
    ]
    panel = JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=models,
        msa=MSAData(enabled=False, query="AAAAAAAA", rows=[MSARow(identifier="query_sequence", sequence="AAAAAAAA", is_query=True)]),
        hydropathy=hydropathy,
        render_config=config,
    )

    svg, _ = render_panel_svg(panel)

    assert svg.count(">α1<") == 1
    assert svg.count(">β1<") == 1
    assert svg.count(">β2<") == 1
    assert "η1" not in svg


def test_confidence_track_renders_one_cell_per_residue():
    config = build_render_config(columns=4, font_size=12)
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", "ALA", "A", str(index + 1))
        for index in range(4)
    ]
    hydropathy = compute_hydropathy(axis, window=3)
    model = ModelTracks(
        name="ranked_0_A",
        source_path="model.pdb",
        chain="A",
        secondary_structure=[SecondaryStructureEntry(index, "C", "coil") for index in range(4)],
        plddt=[95.0, 82.0, 61.0, 24.0],
        accessibility=[AccessibilityEntry(index, None, 0.5, "accessible") for index in range(4)],
        contacts=[ContactEntry(index, None, None, None, None, None, None, None) for index in range(4)],
        display_name="Confidence Demo / 链 A",
    )
    panel = JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=[model],
        msa=MSAData(enabled=False, query="AAAA", rows=[MSARow(identifier="query_sequence", sequence="AAAA", is_query=True)]),
        hydropathy=hydropathy,
        render_config=config,
    )

    svg, _ = render_panel_svg(panel)

    assert svg.count('class="confidence-cell"') == 4
    assert "#174d86" in svg
    assert "#5387aa" in svg
    assert "#c59133" in svg
    assert "#b56a39" in svg


def test_ticks_prefer_round_tens_when_terminal_label_would_overlap():
    panel = _build_tick_panel(list(range(81, 112)))

    svg, _ = render_panel_svg(panel)

    assert ">90<" in svg
    assert ">100<" in svg
    assert ">110<" in svg
    assert ">111<" not in svg


def test_ticks_keep_terminal_label_when_it_does_not_overlap():
    panel = _build_tick_panel(list(range(81, 114)))

    svg, _ = render_panel_svg(panel)

    assert ">110<" in svg
    assert ">113<" in svg
