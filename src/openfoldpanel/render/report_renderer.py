"""HTML shell assembly around static UI assets and rendered SVG panels."""

from __future__ import annotations

import html
import json

from openfoldpanel.constants import SUPPORTED_ANTIBODY_NUMBERING
from openfoldpanel.features.antibody_annotation import antibody_scheme_label
from openfoldpanel.models import JobPanelData, JobReportData, PipelineConfig, RenderConfig
from openfoldpanel.render.svg_renderer import render_panel_svg
from openfoldpanel.render.ui_assets import load_ui_script, load_ui_styles, load_ui_template
from openfoldpanel.utils.text import humanize_chain_label, humanize_job_status, safe_chain_slug, summarize_msa_database_path


def reference_chain_pdf_name(reference_chain: str) -> str:
    """Return the PDF filename for a rendered reference chain."""

    return f"reference-chain-{safe_chain_slug(reference_chain)}.pdf"


def render_reference_chain_report_svg(panel_data: JobPanelData) -> str:
    """Render one FoldScript-style panel SVG for PDF export."""

    panel_svg, _ = render_panel_svg(panel_data, antibody_scheme=panel_data.default_antibody_numbering_scheme)
    return panel_svg


def render_html_report(report_data: JobReportData, config: PipelineConfig) -> str:
    """Assemble the static HTML report from packaged frontend resources."""

    panel_views = [_build_panel_view(panel_data, config) for panel_data in report_data.chain_panels]
    default_view = next(
        (panel for panel in panel_views if panel["referenceChain"] == report_data.default_reference_chain),
        panel_views[0],
    )
    template = load_ui_template()
    replacements = {
        "__OFP_PAGE_TITLE__": html.escape(f"{report_data.job_name} - OpenFoldPanel Report"),
        "__OFP_REPORT_TITLE__": html.escape(report_data.job_name),
        "__OFP_DEFAULT_CHAIN_ID__": html.escape(report_data.default_reference_chain),
        "__OFP_DEFAULT_CHAIN_LABEL__": html.escape(default_view["chainLabel"]),
        "__OFP_DEFAULT_PANEL_WIDTH__": f'{default_view["panelWidth"]:.2f}',
        "__OFP_LEGEND_CONTACT_STRONG_TEXT__": html.escape(_contact_strong_legend_text(config)),
        "__OFP_LEGEND_CONTACT_WEAK_TEXT__": html.escape(_contact_weak_legend_text(config)),
        "__OFP_ANTIBODY_LEGEND_HIDDEN__": "" if default_view["availableAntibodySchemes"] else "hidden",
        "__OFP_THEME_VARS__": _build_theme_vars(report_data.chain_panels[0].render_config),
        "__OFP_INLINE_STYLES__": load_ui_styles(),
        "__OFP_INLINE_SCRIPT__": load_ui_script(),
        "__OFP_REPORT_PAYLOAD__": _serialize_json_for_script(
            {
                "jobName": report_data.job_name,
                "defaultReferenceChain": report_data.default_reference_chain,
                "warnings": list(report_data.warnings),
                "chainPanels": [
                    {
                        "referenceChain": panel["referenceChain"],
                        "chainLabel": panel["chainLabel"],
                        "panelWidth": panel["panelWidth"],
                        "summaryItems": panel["summaryItems"],
                        "warnings": panel["warnings"],
                        "availableAntibodySchemes": panel["availableAntibodySchemes"],
                        "defaultAntibodyScheme": panel["defaultAntibodyScheme"],
                    }
                    for panel in panel_views
                ],
                "defaultAntibodyScheme": panel_views[0]["defaultAntibodyScheme"],
            }
        ),
        "__OFP_CHAIN_TEMPLATES__": "\n".join(
            markup for panel in panel_views for markup in panel["templateMarkups"]
        ),
    }
    for token, value in replacements.items():
        template = template.replace(token, value)
    return template


def _build_panel_view(panel_data: JobPanelData, config: PipelineConfig) -> dict[str, object]:
    available_antibody_schemes = [
        scheme for scheme in SUPPORTED_ANTIBODY_NUMBERING if scheme in panel_data.antibody_numberings
    ]
    template_schemes = available_antibody_schemes or [panel_data.default_antibody_numbering_scheme]
    default_antibody_scheme = (
        panel_data.default_antibody_numbering_scheme
        if panel_data.default_antibody_numbering_scheme in template_schemes
        else template_schemes[0]
    )
    template_markups: list[str] = []
    panel_width = 0.0
    for scheme in template_schemes:
        panel_svg, layout = render_panel_svg(panel_data, antibody_scheme=scheme)
        if panel_width == 0.0:
            panel_width = layout.width
        template_markups.append(
            _render_chain_template(
                reference_chain=panel_data.reference_chain,
                chain_label=humanize_chain_label(panel_data.reference_chain),
                panel_width=layout.width,
                panel_svg=_strip_xml_declaration(panel_svg),
                antibody_scheme=scheme,
            )
        )
    chain_label = humanize_chain_label(panel_data.reference_chain)
    return {
        "referenceChain": panel_data.reference_chain,
        "chainLabel": chain_label,
        "panelWidth": round(panel_width, 2),
        "summaryItems": _build_summary_items(panel_data, config),
        "warnings": list(panel_data.warnings),
        "availableAntibodySchemes": list(available_antibody_schemes),
        "defaultAntibodyScheme": default_antibody_scheme,
        "templateMarkups": template_markups,
    }


def _render_chain_template(
    *,
    reference_chain: str,
    chain_label: str,
    panel_width: float,
    panel_svg: str,
    antibody_scheme: str,
) -> str:
    return (
        f'<template id="ofp-chain-{safe_chain_slug(reference_chain)}-{safe_chain_slug(antibody_scheme)}" '
        f'data-chain-figure="{html.escape(reference_chain)}" '
        f'data-chain-label="{html.escape(chain_label)}" '
        f'data-antibody-scheme="{html.escape(antibody_scheme)}" '
        f'data-antibody-scheme-label="{html.escape(antibody_scheme_label(antibody_scheme))}" '
        f'data-panel-width="{panel_width:.2f}">'
        f"{panel_svg}"
        "</template>"
    )


def _build_summary_items(panel_data: JobPanelData, config: PipelineConfig) -> list[dict[str, str]]:
    return [
        {"label": "Reference Chain", "value": humanize_chain_label(panel_data.reference_chain)},
        {"label": "Residue Span", "value": _sequence_span_label(panel_data)},
        {"label": "Model Count", "value": str(len(panel_data.models))},
        {"label": "Output Status", "value": humanize_job_status(panel_data.status)},
        {
            "label": "Hydropathy Window",
            "value": str(config.hyd_window),
            "tooltip": "Sliding-window size for Kyte-Doolittle hydropathy. CLI flag: --hyd-window",
        },
        {
            "label": "E-value Threshold",
            "value": config.evalue,
            "tooltip": "Significance cutoff used to filter homolog hits. Smaller values are stricter. CLI flag: --evalue",
        },
        {
            "label": "Weak Contact Cutoff",
            "value": f"{config.contact_cutoff:g} A",
            "tooltip": "Distance cutoff used to classify weak contacts, in angstroms. CLI flag: --contact-cutoff",
        },
        {
            "label": "Strong Contact Cutoff",
            "value": f"{config.strong_contact_cutoff:g} A",
            "tooltip": "Distance cutoff used to classify strong contacts, in angstroms. CLI flag: --strong-contact-cutoff",
        },
        {
            "label": "Homolog Display Limit",
            "value": str(config.max_homologs_displayed),
            "tooltip": "Maximum number of homolog sequences to retrieve and display. CLI flag: --max-homologs-displayed",
        },
        {
            "label": "Database",
            "value": summarize_msa_database_path(config.msa_db),
            "tooltip": "Database label derived from the tail of the configured path. CLI flag: --msa-db",
        },
    ]


def _sequence_span_label(panel_data: JobPanelData) -> str:
    if not panel_data.sequence_axis:
        return "No residues available"
    first = panel_data.sequence_axis[0].label
    last = panel_data.sequence_axis[-1].label
    return f"{first} - {last}"


def _contact_strong_legend_text(config: PipelineConfig) -> str:
    return f"Based on the shortest non-hydrogen atom distance: below {config.strong_contact_cutoff:g} A."


def _contact_weak_legend_text(config: PipelineConfig) -> str:
    return (
        f"Based on the shortest non-hydrogen atom distance: "
        f"between {config.strong_contact_cutoff:g} A and {config.contact_cutoff:g} A, inclusive."
    )


def _strip_xml_declaration(svg_markup: str) -> str:
    if svg_markup.startswith("<?xml"):
        return svg_markup.split("\n", 1)[1]
    return svg_markup


def _serialize_json_for_script(payload: dict[str, object]) -> str:
    text = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    return text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def _build_theme_vars(config: RenderConfig) -> str:
    colors = config.colors
    return "\n".join(
        [
            f"  --ofp-font-sans: {config.font_family};",
            f"  --ofp-font-serif: {config.heading_font_family};",
            '  --ofp-font-mono: "Liberation Mono", "Nimbus Mono PS", "Courier New", monospace;',
            f"  --ofp-color-page: {colors['background']};",
            f"  --ofp-color-surface: {colors['surface']};",
            f"  --ofp-color-border: {colors['border']};",
            f"  --ofp-color-border-strong: {colors['accent_border']};",
            f"  --ofp-color-grid: {colors['grid']};",
            f"  --ofp-color-text: {colors['text']};",
            f"  --ofp-color-muted: {colors['muted_text']};",
            f"  --ofp-color-accent: {colors['accent']};",
            f"  --ofp-color-accent-soft: {colors['accent_soft']};",
            f"  --ofp-color-warning: {colors['warning']};",
            f"  --ofp-color-warning-bg: {colors['warning_bg']};",
            f"  --ofp-color-warning-border: {colors['warning_border']};",
            f"  --ofp-color-alpha-helix: {colors['alpha_helix']};",
            f"  --ofp-color-three-ten-helix: {colors['three_ten_helix']};",
            f"  --ofp-color-pi-helix: {colors['pi_helix']};",
            f"  --ofp-color-strand-fill: {colors['strand_fill']};",
            f"  --ofp-color-strand-stroke: {colors['strand_stroke']};",
            f"  --ofp-color-alpha-turn: {colors['alpha_turn_text']};",
            f"  --ofp-color-beta-turn: {colors['beta_turn_text']};",
            f"  --ofp-color-msa-query-bg: {colors['msa_query_bg']};",
            f"  --ofp-color-msa-query-text: {colors['msa_query_text']};",
            f"  --ofp-color-msa-identity-bg: {colors['msa_identity_bg']};",
            f"  --ofp-color-msa-identity-text: {colors['msa_identity_text']};",
            f"  --ofp-color-msa-similar-bg: {colors['msa_similar_bg']};",
            f"  --ofp-color-msa-similar-text: {colors['msa_similar_text']};",
            f"  --ofp-color-msa-default-bg: {colors['msa_default_bg']};",
            f"  --ofp-color-msa-default-text: {colors['msa_default_text']};",
            f"  --ofp-color-accessibility-buried: {colors['accessibility_buried']};",
            f"  --ofp-color-accessibility-intermediate: {colors['accessibility_intermediate']};",
            f"  --ofp-color-accessibility-accessible: {colors['accessibility_accessible']};",
            f"  --ofp-color-accessibility-highly-exposed: {colors['accessibility_highly_exposed']};",
            f"  --ofp-color-hydrophobic: {colors['hydropathy_hydrophobic']};",
            f"  --ofp-color-hydropathy-intermediate: {colors['hydropathy_intermediate']};",
            f"  --ofp-color-hydrophilic: {colors['hydropathy_hydrophilic']};",
            f"  --ofp-color-confidence-very-high: {colors['plddt_very_high']};",
            f"  --ofp-color-confidence-confident: {colors['plddt_confident']};",
            f"  --ofp-color-confidence-low: {colors['plddt_low']};",
            f"  --ofp-color-confidence-very-low: {colors['plddt_very_low']};",
            f"  --ofp-color-contact-bg: {colors['contact_bg']};",
            f"  --ofp-color-contact-strong: {colors['contact_strong']};",
            f"  --ofp-color-contact-weak: {colors['contact_weak']};",
            f"  --ofp-color-contact-multi-outline: {colors['contact_multi_outline']};",
            f"  --ofp-color-disulfide: {colors['disulfide_symbol']};",
            f"  --ofp-color-disulfide-inter: {colors['disulfide_inter_symbol']};",
            f"  --ofp-color-atmosphere: {colors['accent_border']};",
        ]
    )
