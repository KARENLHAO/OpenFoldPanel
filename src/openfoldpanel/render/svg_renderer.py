"""Native SVG rendering for FoldScript-style flat panels."""

from __future__ import annotations

import html
from collections.abc import Iterable
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path

from openfoldpanel.constants import PLDDT_THRESHOLDS
from openfoldpanel.models import ContactEntry, JobPanelData, SecondaryStructureEntry
from openfoldpanel.render.font_assets import embedded_times_new_roman_css
from openfoldpanel.render.glyphs import helix_path, strand_points, turn_path
from openfoldpanel.render.layout import LayoutBlock, LayoutRow, PanelLayout, build_panel_layout
from openfoldpanel.utils.residue_utils import compatible_similarity_group


@dataclass(slots=True)
class StructureAnnotation:
    category: str
    start: int
    end: int
    label: str


@dataclass(slots=True)
class TickCandidate:
    axis_index: int
    x: float
    label: str
    priority: int


def render_panel_svg(panel_data: JobPanelData) -> tuple[str, PanelLayout]:
    """Render a panel into SVG markup without writing to disk."""

    layout = build_panel_layout(panel_data)
    svg = _render_svg_string(panel_data, layout)
    return svg, layout


def render_svg(panel_data: JobPanelData, output_path: Path) -> PanelLayout:
    """Render the panel to an SVG file."""

    svg, layout = render_panel_svg(panel_data)
    output_path.write_text(svg, encoding="utf-8")
    return layout


def _render_svg_string(panel_data: JobPanelData, layout: PanelLayout) -> str:
    config = layout.render_config
    pieces = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{layout.width:.2f}" height="{layout.height:.2f}" viewBox="0 0 {layout.width:.2f} {layout.height:.2f}">',
        f'<rect width="100%" height="100%" fill="{config.colors["surface"]}"/>',
        _style_block(config),
    ]
    for block in layout.blocks:
        pieces.append(_render_block(panel_data, layout, block))
    pieces.append("</svg>")
    return "\n".join(pieces)


def _style_block(config) -> str:
    shared_font = config.font_family
    return (
        "<style>"
        f"{embedded_times_new_roman_css()}"
        f'.model-label{{font-family:{shared_font};font-size:{config.font_size + 0.5}px;fill:{config.colors["strand_fill"]};font-style:italic;font-weight:700;}}'
        f'.track-label{{font-family:{config.font_family};font-size:{config.font_size}px;fill:{config.colors["strand_fill"]};font-style:italic;font-weight:700;}}'
        f'.homolog-label{{font-family:{config.font_family};font-size:{config.font_size - 0.3}px;fill:{config.colors["text"]};font-weight:600;}}'
        f'.annotation-label{{font-family:{config.heading_font_family};font-size:{config.font_size + 1.2}px;font-weight:700;}}'
        f'.tick-label{{font-family:{shared_font};font-size:{config.font_size}px;fill:{config.colors["strand_fill"]};font-weight:700;}}'
        f'.sequence-text{{font-family:{shared_font};font-size:{config.font_size + 0.9}px;font-weight:700;dominant-baseline:middle;text-anchor:middle;}}'
        f'.contact-text{{font-family:{shared_font};font-size:{config.font_size + 0.9}px;font-weight:700;dominant-baseline:middle;text-anchor:middle;}}'
        '.confidence-cell{shape-rendering:crispEdges;}'
        f'.turn-track{{fill:none;stroke:{config.colors["turn_text"]};stroke-width:1.9;stroke-linecap:round;stroke-linejoin:round;}}'
        f'.soft-rule{{stroke:{config.colors["grid"]};stroke-width:0.9;}}'
        "</style>"
    )


def _render_block(panel_data: JobPanelData, layout: PanelLayout, block: LayoutBlock) -> str:
    config = layout.render_config
    grid_x = config.margin + config.label_width
    pieces = [f'<g transform="translate({block.x:.2f},{block.y:.2f})">']

    pieces.append(_render_structure_annotation(panel_data, block, layout))
    pieces.append(_render_sequence_ticks(panel_data, block, layout))

    for row, y, height in zip(layout.rows, layout.row_positions, layout.row_heights, strict=True):
        pieces.append(_render_row_label(row, y, height, config))
        pieces.append(_render_row(panel_data, row, block, y, height, layout))

    grid_end = grid_x + (block.end - block.start) * config.cell_width
    pieces.append(f'<line class="soft-rule" x1="{grid_x:.2f}" y1="{layout.tick_y + layout.tick_height:.2f}" x2="{grid_end:.2f}" y2="{layout.tick_y + layout.tick_height:.2f}"/>')
    pieces.append("</g>")
    return "\n".join(piece for piece in pieces if piece)


def _render_row_label(row: LayoutRow, y: float, height: float, config) -> str:
    label_class = "model-label" if row.kind in {"secondary", "contacts"} else "track-label"
    if row.kind == "msa_homolog":
        label_class = "homolog-label"
    label_y = y + height * 0.74
    return f'<text class="{label_class}" x="{config.margin:.2f}" y="{label_y:.2f}">{html.escape(row.label)}</text>'


def _render_structure_annotation(panel_data: JobPanelData, block: LayoutBlock, layout: PanelLayout) -> str:
    if not panel_data.models:
        return ""
    config = layout.render_config
    grid_x = config.margin + config.label_width
    band_y = layout.annotation_y
    pieces: list[str] = []
    for segment in _secondary_annotations(panel_data.models[0].secondary_structure):
        if segment.end <= block.start or segment.start >= block.end:
            continue
        clipped_start = max(segment.start, block.start)
        clipped_end = min(segment.end, block.end)
        x = grid_x + (clipped_start - block.start) * config.cell_width
        width = (clipped_end - clipped_start) * config.cell_width
        label_x = x + width / 2.0
        color = _structure_color(segment.category, config.colors)
        pieces.append(
            f'<text class="annotation-label" x="{label_x:.2f}" y="{band_y + config.font_size + 1.4:.2f}" '
            f'fill="{color}" text-anchor="middle">{html.escape(segment.label)}</text>'
        )
    return "\n".join(pieces)


def _render_sequence_ticks(panel_data: JobPanelData, block: LayoutBlock, layout: PanelLayout) -> str:
    config = layout.render_config
    grid_x = config.margin + config.label_width
    tick_base_y = layout.tick_y + layout.tick_height - 2
    query_row_index = next((index for index, row in enumerate(layout.rows) if row.kind == "msa_query"), None)
    query_top = layout.row_positions[query_row_index] if query_row_index is not None else layout.tick_y + layout.tick_height + 2
    candidates_by_axis: dict[int, TickCandidate] = {}
    for local_index, axis_index in enumerate(range(block.start, block.end)):
        position = panel_data.sequence_axis[axis_index]
        x = grid_x + local_index * config.cell_width
        if local_index == 0:
            _register_tick_candidate(candidates_by_axis, TickCandidate(axis_index=axis_index, x=x, label=position.label, priority=2))
        if position.seq_id % 10 == 0:
            _register_tick_candidate(candidates_by_axis, TickCandidate(axis_index=axis_index, x=x, label=position.label, priority=3))
        if axis_index == block.end - 1:
            _register_tick_candidate(candidates_by_axis, TickCandidate(axis_index=axis_index, x=x, label=position.label, priority=1))

    pieces: list[str] = []
    for candidate in _select_tick_candidates(candidates_by_axis.values(), config.font_size):
        pieces.append(f'<text class="tick-label" x="{candidate.x + 1:.2f}" y="{tick_base_y:.2f}">{html.escape(candidate.label)}</text>')
        pieces.append(
            f'<line class="soft-rule" x1="{candidate.x:.2f}" y1="{layout.tick_y + 1:.2f}" x2="{candidate.x:.2f}" y2="{query_top - 2:.2f}"/>'
        )
    return "\n".join(pieces)


def _register_tick_candidate(candidates_by_axis: dict[int, TickCandidate], candidate: TickCandidate) -> None:
    existing = candidates_by_axis.get(candidate.axis_index)
    if existing is None or candidate.priority > existing.priority:
        candidates_by_axis[candidate.axis_index] = candidate


def _select_tick_candidates(candidates: Iterable[TickCandidate], font_size: int) -> list[TickCandidate]:
    accepted: list[TickCandidate] = []
    ordered_candidates = sorted(candidates, key=lambda item: (-item.priority, item.axis_index))
    for candidate in ordered_candidates:
        if any(_tick_candidates_overlap(candidate, other, font_size) for other in accepted):
            continue
        accepted.append(candidate)
    return sorted(accepted, key=lambda item: item.axis_index)


def _tick_candidates_overlap(left: TickCandidate, right: TickCandidate, font_size: int) -> bool:
    left_start, left_end = _tick_label_bounds(left, font_size)
    right_start, right_end = _tick_label_bounds(right, font_size)
    return left_start < right_end and right_start < left_end


def _tick_label_bounds(candidate: TickCandidate, font_size: int) -> tuple[float, float]:
    text_start = candidate.x + 1.0
    text_width = len(candidate.label) * font_size * 0.64 + 2.0
    return text_start, text_start + text_width


def _render_row(panel_data: JobPanelData, row: LayoutRow, block: LayoutBlock, y: float, height: float, layout: PanelLayout) -> str:
    if row.kind == "secondary":
        return _render_secondary_row(panel_data, row.model_index, block, y, height, layout)
    if row.kind == "confidence":
        return _render_confidence_row(panel_data, row.model_index, block, y, height, layout)
    if row.kind in {"msa_query", "msa_homolog"}:
        return _render_msa_row(panel_data, row, block, y, height, layout)
    if row.kind == "accessibility":
        return _render_accessibility_row(panel_data, block, y, height, layout)
    if row.kind == "hydropathy":
        return _render_hydropathy_row(panel_data, block, y, height, layout)
    if row.kind == "contacts":
        return _render_contacts_row(panel_data, row.model_index, block, y, height, layout)
    return ""


def _segment_ranges(track: list[SecondaryStructureEntry], block: LayoutBlock, category: str) -> list[tuple[int, int]]:
    indices = [entry.residue_index for entry in track[block.start:block.end] if entry.category == category]
    segments: list[tuple[int, int]] = []
    for _, group in groupby(enumerate(indices), lambda pair: pair[1] - pair[0]):
        chunk = [item for _, item in group]
        segments.append((chunk[0], chunk[-1] + 1))
    return segments


def _render_secondary_row(panel_data: JobPanelData, model_index: int | None, block: LayoutBlock, y: float, height: float, layout: PanelLayout) -> str:
    if model_index is None:
        return ""
    config = layout.render_config
    model = panel_data.models[model_index]
    grid_x = config.margin + config.label_width
    row_y = y + height * 0.16
    row_height = height * 0.62
    pieces = []

    helix_segments = _segment_ranges(model.secondary_structure, block, "helix")
    strand_segments = _segment_ranges(model.secondary_structure, block, "strand")
    turn_segments = _segment_ranges(model.secondary_structure, block, "turn")

    for start, end in helix_segments:
        x = grid_x + (start - block.start) * config.cell_width
        width = (end - start) * config.cell_width
        pieces.append(f'<path d="{helix_path(x, row_y, width, row_height)}" fill="none" stroke="{config.colors["helix_fill"]}" stroke-width="1.9"/>')
    for start, end in strand_segments:
        x = grid_x + (start - block.start) * config.cell_width
        width = (end - start) * config.cell_width
        pieces.append(f'<polygon points="{strand_points(x, row_y, width, row_height)}" fill="{config.colors["strand_fill"]}" stroke="{config.colors["strand_stroke"]}" stroke-width="0.8"/>')
    turn_height = row_height * 0.46
    turn_y = y + (height - turn_height) / 2.0
    turn_padding = min(max(config.cell_width * 0.12, 0.6), config.cell_width * 0.24)
    for start, end in turn_segments:
        x = grid_x + (start - block.start) * config.cell_width + turn_padding
        width = max((end - start) * config.cell_width - turn_padding * 2, config.cell_width * 0.56)
        pieces.append(f'<path class="turn-track" d="{turn_path(x, turn_y, width, turn_height)}"/>')
    return "\n".join(pieces)


def _render_confidence_row(panel_data: JobPanelData, model_index: int | None, block: LayoutBlock, y: float, height: float, layout: PanelLayout) -> str:
    if model_index is None:
        return ""
    config = layout.render_config
    grid_x = config.margin + config.label_width
    plddt = panel_data.models[model_index].plddt
    pieces: list[str] = []
    for local_index, axis_index in enumerate(range(block.start, block.end)):
        if axis_index >= len(plddt):
            continue
        x = grid_x + local_index * config.cell_width
        color = _plddt_color(plddt[axis_index], config.colors)
        pieces.append(
            f'<rect class="confidence-cell" x="{x:.2f}" y="{y + 1.0:.2f}" width="{config.cell_width:.2f}" height="{height - 2.0:.2f}" fill="{color}" stroke="none"/>'
        )
    return "\n".join(pieces)


def _render_msa_row(panel_data: JobPanelData, row: LayoutRow, block: LayoutBlock, y: float, height: float, layout: PanelLayout) -> str:
    if row.msa_row_index is None or row.msa_row_index >= len(panel_data.msa.rows):
        return ""
    config = layout.render_config
    msa_row = panel_data.msa.rows[row.msa_row_index]
    query_row = next((entry for entry in panel_data.msa.rows if entry.is_query), panel_data.msa.rows[0])
    grid_x = config.margin + config.label_width
    pieces: list[str] = []
    for local_index, axis_index in enumerate(range(block.start, block.end)):
        if axis_index >= len(msa_row.sequence):
            continue
        residue = msa_row.sequence[axis_index]
        if (
            not msa_row.is_query
            and axis_index == 0
            and row.msa_row_index < len(panel_data.msa.leading_display_overrides)
            and panel_data.msa.leading_display_overrides[row.msa_row_index] is not None
            and residue == "-"
        ):
            residue = panel_data.msa.leading_display_overrides[row.msa_row_index] or residue
        query_residue = query_row.sequence[axis_index] if axis_index < len(query_row.sequence) else residue
        x = grid_x + local_index * config.cell_width
        bg_color, text_color = _msa_style(residue, query_residue, row.kind, config.colors)
        pieces.append(
            f'<rect x="{x:.2f}" y="{y + 0.8:.2f}" width="{config.cell_width:.2f}" height="{height - 1.6:.2f}" fill="{bg_color}" stroke="{config.colors["grid"]}" stroke-width="0.45"/>'
        )
        pieces.append(
            f'<text class="sequence-text" x="{x + config.cell_width / 2.0:.2f}" y="{y + height / 2.0 + 0.35:.2f}" fill="{text_color}">{html.escape(residue)}</text>'
        )
    return "\n".join(pieces)


def _render_accessibility_row(panel_data: JobPanelData, block: LayoutBlock, y: float, height: float, layout: PanelLayout) -> str:
    config = layout.render_config
    grid_x = config.margin + config.label_width
    track = panel_data.models[0].accessibility if panel_data.models else []
    pieces: list[str] = []
    for local_index, axis_index in enumerate(range(block.start, block.end)):
        if axis_index >= len(track):
            continue
        entry = track[axis_index]
        x = grid_x + local_index * config.cell_width
        color = _accessibility_color(entry.category, config.colors)
        pieces.append(
            f'<rect x="{x:.2f}" y="{y + 1.0:.2f}" width="{config.cell_width:.2f}" height="{height - 2.0:.2f}" fill="{color}" stroke="none"/>'
        )
    return "\n".join(pieces)


def _render_hydropathy_row(panel_data: JobPanelData, block: LayoutBlock, y: float, height: float, layout: PanelLayout) -> str:
    config = layout.render_config
    grid_x = config.margin + config.label_width
    pieces: list[str] = []
    for local_index, axis_index in enumerate(range(block.start, block.end)):
        if axis_index >= len(panel_data.hydropathy):
            continue
        entry = panel_data.hydropathy[axis_index]
        x = grid_x + local_index * config.cell_width
        pieces.append(
            f'<rect x="{x:.2f}" y="{y + 1.0:.2f}" width="{config.cell_width:.2f}" height="{height - 2.0:.2f}" fill="{_hydropathy_color(entry.category, config.colors)}" stroke="none"/>'
        )
    return "\n".join(pieces)


def _render_contacts_row(panel_data: JobPanelData, model_index: int | None, block: LayoutBlock, y: float, height: float, layout: PanelLayout) -> str:
    if model_index is None:
        return ""
    config = layout.render_config
    grid_x = config.margin + config.label_width
    model = panel_data.models[model_index]
    track = model.contacts
    disulfide_scopes = _disulfide_scopes_by_residue(model.disulfides)
    pieces: list[str] = []
    for local_index, axis_index in enumerate(range(block.start, block.end)):
        if axis_index >= len(track):
            continue
        entry = track[axis_index]
        x = grid_x + local_index * config.cell_width
        disulfide_scope = disulfide_scopes.get(axis_index)
        if disulfide_scope is not None:
            pieces.append(_render_contact_cell_background(x, y, height, config))
            pieces.append(
                _render_contact_symbol(
                    entry,
                    x,
                    y,
                    height,
                    config,
                    symbol_override="S",
                    color_override=_disulfide_color(disulfide_scope, config.colors),
                )
            )
        elif entry.symbol:
            pieces.append(_render_contact_cell_background(x, y, height, config))
            pieces.append(_render_contact_symbol(entry, x, y, height, config))
        elif entry.is_multi_contact:
            pieces.append(
                f'<rect x="{x + 0.4:.2f}" y="{y + 0.9:.2f}" width="{config.cell_width - 0.8:.2f}" height="{height - 1.8:.2f}" fill="none" stroke="{config.colors["contact_multi_outline"]}" stroke-width="0.9"/>'
            )
    return "\n".join(pieces)


def _render_contact_cell_background(x: float, y: float, height: float, config) -> str:
    return (
        f'<rect x="{x + 0.2:.2f}" y="{y + 0.7:.2f}" width="{config.cell_width - 0.4:.2f}" '
        f'height="{height - 1.4:.2f}" fill="{config.colors["contact_bg"]}" stroke="none"/>'
    )


def _render_contact_symbol(
    entry: ContactEntry,
    x: float,
    y: float,
    height: float,
    config,
    *,
    symbol_override: str | None = None,
    color_override: str | None = None,
) -> str:
    symbol = symbol_override if symbol_override is not None else (entry.symbol or "")
    color = color_override or (
        config.colors["contact_strong"] if entry.strength_category == "strong" else config.colors["contact_weak"]
    )
    outline = ""
    if entry.is_multi_contact:
        outline = (
            f'<rect x="{x + 0.4:.2f}" y="{y + 0.7:.2f}" width="{config.cell_width - 0.8:.2f}" height="{height - 1.4:.2f}" '
            f'fill="none" stroke="{config.colors["contact_multi_outline"]}" stroke-width="0.85"/>'
        )
    text = (
        f'<text class="contact-text" x="{x + config.cell_width / 2.0:.2f}" '
        f'y="{y + height / 2.0 + 0.35:.2f}" fill="{color}">{html.escape(symbol)}</text>'
    )
    return f"{outline}\n{text}" if outline else text


def _disulfide_scopes_by_residue(disulfides) -> dict[int, str]:
    scopes: dict[int, str] = {}
    for bond in disulfides:
        scopes[bond.residue_index_a] = bond.bridge_scope
        if bond.residue_index_b is not None:
            scopes[bond.residue_index_b] = bond.bridge_scope
    return scopes


def _disulfide_color(scope: str, colors: dict[str, str]) -> str:
    if scope == "intermolecular":
        return colors["disulfide_inter_symbol"]
    return colors["disulfide_symbol"]


def _secondary_annotations(track: list[SecondaryStructureEntry]) -> list[StructureAnnotation]:
    annotations: list[StructureAnnotation] = []
    counters = {"strand": 0, "helix": 0}
    index = 0
    while index < len(track):
        category = track[index].category
        if category not in {"strand", "helix"}:
            index += 1
            continue
        end = index + 1
        while end < len(track) and track[end].category == category:
            end += 1
        counters[category] += 1
        symbol = {"strand": "β", "helix": "α"}[category]
        annotations.append(StructureAnnotation(category=category, start=index, end=end, label=f"{symbol}{counters[category]}"))
        index = end
    return annotations


def _structure_color(category: str, colors: dict[str, str]) -> str:
    if category == "strand":
        return colors["strand_fill"]
    if category == "helix":
        return colors["helix_fill"]
    return colors["turn_text"]


def _msa_style(residue: str, query_residue: str, kind: str, colors: dict[str, str]) -> tuple[str, str]:
    if kind == "msa_query":
        return colors["msa_query_bg"], colors["msa_query_text"]
    if residue == query_residue:
        return colors["msa_identity_bg"], colors["msa_identity_text"]
    if compatible_similarity_group({residue, query_residue}):
        return colors["msa_similar_bg"], colors["msa_similar_text"]
    return colors["msa_default_bg"], colors["msa_default_text"]


def _accessibility_color(category: str | None, colors: dict[str, str]) -> str:
    if category == "buried":
        return colors["accessibility_buried"]
    if category == "intermediate":
        return colors["accessibility_intermediate"]
    if category == "accessible":
        return colors["accessibility_accessible"]
    if category == "highly_exposed":
        return colors["accessibility_highly_exposed"]
    return colors["coil"]


def _hydropathy_color(category: str | None, colors: dict[str, str]) -> str:
    if category == "hydrophobic":
        return colors["hydropathy_hydrophobic"]
    if category == "hydrophilic":
        return colors["hydropathy_hydrophilic"]
    return colors["hydropathy_intermediate"]


def _plddt_color(value: float | None, colors: dict[str, str]) -> str:
    if value is None:
        return colors["coil"]
    if value >= PLDDT_THRESHOLDS["very_high"]:
        return colors["plddt_very_high"]
    if value >= PLDDT_THRESHOLDS["confident"]:
        return colors["plddt_confident"]
    if value >= PLDDT_THRESHOLDS["low"]:
        return colors["plddt_low"]
    return colors["plddt_very_low"]
