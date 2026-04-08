"""Layout planning for FoldScript-style flat panels."""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import replace

from openfoldpanel.constants import (
    COLORS,
    DEFAULT_BLOCK_GAP,
    DEFAULT_CELL_WIDTH_RATIO,
    DEFAULT_FONT_FAMILY,
    DEFAULT_HEADING_FONT_FAMILY,
    DEFAULT_LABEL_COLUMNS,
    DEFAULT_MARGIN,
    DEFAULT_TICK_HEIGHT,
)
from openfoldpanel.models import JobPanelData, RenderConfig


@dataclass(slots=True)
class LayoutRow:
    kind: str
    label: str
    model_index: int | None = None
    msa_row_index: int | None = None


@dataclass(slots=True)
class LayoutBlock:
    block_index: int
    start: int
    end: int
    x: float
    y: float
    width: float
    height: float


@dataclass(slots=True)
class PanelLayout:
    width: float
    height: float
    rows: list[LayoutRow]
    blocks: list[LayoutBlock]
    row_positions: list[float]
    row_heights: list[float]
    annotation_y: float
    annotation_height: float
    tick_y: float
    tick_height: float
    render_config: RenderConfig


def build_render_config(columns: int, font_size: int, max_homologs_displayed: int = 5) -> RenderConfig:
    """Build the shared render configuration from CLI parameters."""

    cell_width = round(font_size * max(DEFAULT_CELL_WIDTH_RATIO, 1.08), 2)
    row_height = round(font_size * 1.42, 2)
    label_width = round(cell_width * DEFAULT_LABEL_COLUMNS, 2)
    return RenderConfig(
        columns=columns,
        max_homologs_displayed=max_homologs_displayed,
        font_size=font_size,
        cell_width=cell_width,
        row_height=row_height,
        label_width=label_width,
        margin=DEFAULT_MARGIN,
        colors=dict(COLORS),
        font_family=DEFAULT_FONT_FAMILY,
        heading_font_family=DEFAULT_HEADING_FONT_FAMILY,
    )


def build_rows(panel_data: JobPanelData) -> list[LayoutRow]:
    """Create ordered rows for the FoldScript-style panel."""

    rows: list[LayoutRow] = []
    for index, model in enumerate(panel_data.models):
        rows.append(LayoutRow(kind="secondary", label=model.name, model_index=index))

    if panel_data.antibody_numberings:
        rows.append(LayoutRow(kind="antibody_numbering", label="Antibody Numbering"))

    for index, _model in enumerate(panel_data.models):
        rows.append(LayoutRow(kind="confidence", label="Confidence", model_index=index))

    query_index = next((index for index, row in enumerate(panel_data.msa.rows) if row.is_query), 0)
    rows.append(LayoutRow(kind="msa_query", label="Query Sequence", msa_row_index=query_index))

    for homolog_index, row_index in enumerate(_displayed_homolog_indices(panel_data), start=1):
        homolog_row = panel_data.msa.rows[row_index]
        label = homolog_row.identifier.strip() or f"Homolog {homolog_index}"
        rows.append(LayoutRow(kind="msa_homolog", label=label, msa_row_index=row_index))

    rows.append(LayoutRow(kind="accessibility", label="Accessibility"))
    rows.append(LayoutRow(kind="hydropathy", label="Hydropathy"))

    for index, model in enumerate(panel_data.models):
        rows.append(LayoutRow(kind="contacts", label=model.name, model_index=index))
    return rows


def build_panel_layout(panel_data: JobPanelData) -> PanelLayout:
    """Build block positions and panel dimensions for FoldScript-like rendering."""

    config = panel_data.render_config
    rows = build_rows(panel_data)
    label_width = max(config.label_width, _estimate_label_width(rows, config))
    if label_width != config.label_width:
        config = replace(config, label_width=label_width)

    row_heights = [_row_height(row.kind, config) for row in rows]
    annotation_height = _annotation_height(config)
    tick_height = DEFAULT_TICK_HEIGHT
    annotation_y = config.margin

    row_positions: list[float] = []
    cursor = annotation_y + annotation_height
    previous_group: str | None = None
    tick_y = cursor
    for row, height in zip(rows, row_heights, strict=True):
        group = _row_group(row.kind)
        if previous_group is None:
            cursor += _section_gap(config, 0.22, 3.0)
        elif previous_group != group:
            if previous_group == "secondary" and group == "msa":
                cursor += _section_gap(config, 0.30, 4.0)
                tick_y = cursor
                cursor += tick_height + _section_gap(config, 0.22, 3.0)
            else:
                cursor += _section_gap(config, 0.42, 6.0)
        else:
            cursor += _row_spacing(row.kind, config)
        row_positions.append(cursor)
        cursor += height
        previous_group = group

    sequence_length = len(panel_data.sequence_axis)
    columns = max(1, config.columns)
    num_blocks = max(1, (sequence_length + columns - 1) // columns)
    right_gutter = _right_gutter(config)
    block_width = config.margin * 2 + config.label_width + columns * config.cell_width + right_gutter
    block_height = cursor + config.margin

    blocks: list[LayoutBlock] = []
    y = config.margin
    for block_index in range(num_blocks):
        start = block_index * columns
        end = min(sequence_length, start + columns)
        width = config.margin * 2 + config.label_width + (end - start) * config.cell_width + right_gutter
        blocks.append(
            LayoutBlock(
                block_index=block_index,
                start=start,
                end=end,
                x=config.margin,
                y=y,
                width=width,
                height=block_height,
            )
        )
        y += block_height + DEFAULT_BLOCK_GAP

    total_width = max(block_width, max(block.width for block in blocks))
    total_height = blocks[-1].y + blocks[-1].height + config.margin if blocks else block_height + config.margin * 2
    return PanelLayout(
        width=total_width,
        height=total_height,
        rows=rows,
        blocks=blocks,
        row_positions=row_positions,
        row_heights=row_heights,
        annotation_y=annotation_y,
        annotation_height=annotation_height,
        tick_y=tick_y,
        tick_height=tick_height,
        render_config=config,
    )


def _displayed_homolog_indices(panel_data: JobPanelData) -> list[int]:
    max_rows = max(0, panel_data.render_config.max_homologs_displayed)
    if max_rows == 0:
        return []
    homolog_rows = [index for index, row in enumerate(panel_data.msa.rows) if not row.is_query]
    return homolog_rows[:max_rows]


def _estimate_label_width(rows: list[LayoutRow], config: RenderConfig) -> float:
    longest_label = max((len(row.label) for row in rows), default=DEFAULT_LABEL_COLUMNS)
    return round(max(config.label_width, longest_label * config.font_size * 0.72), 2)


def _row_group(kind: str) -> str:
    if kind in {"secondary", "confidence", "antibody_numbering"}:
        return "secondary"
    if kind.startswith("msa_"):
        return "msa"
    if kind in {"accessibility", "hydropathy"}:
        return "tracks"
    return "contacts"


def _row_height(kind: str, config: RenderConfig) -> float:
    if kind == "secondary":
        return round(config.font_size * 1.52, 2)
    if kind == "antibody_numbering":
        return round(max(config.font_size * 1.62, 18.0), 2)
    if kind == "confidence":
        return round(max(config.font_size * 0.82, 8.4), 2)
    if kind.startswith("msa_"):
        return round(config.font_size * 1.54, 2)
    if kind in {"accessibility", "hydropathy"}:
        return round(max(config.font_size * 0.78, 8.0), 2)
    return round(config.font_size * 1.46, 2)


def _row_spacing(kind: str, config: RenderConfig) -> float:
    if kind in {"confidence", "accessibility", "hydropathy"}:
        return 2.0
    return _section_gap(config, 0.08, 1.5)


def _annotation_height(config: RenderConfig) -> float:
    return round(max(config.font_size * 1.28, 15.0), 2)


def _section_gap(config: RenderConfig, ratio: float, minimum: float) -> float:
    return round(max(config.font_size * ratio, minimum), 2)


def _right_gutter(config: RenderConfig) -> float:
    """Reserve trailing space so end-of-block labels do not clip at the SVG edge."""

    return round(max(config.margin, config.font_size * 2.6), 2)
