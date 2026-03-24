"""HTML shell assembly around static UI assets and rendered SVG panels."""

from __future__ import annotations

import html
import json
from functools import lru_cache
from importlib import resources

from openfoldpanel.models import JobPanelData, JobReportData, PipelineConfig, RenderConfig
from openfoldpanel.render.svg_renderer import render_panel_svg
from openfoldpanel.utils.text import humanize_chain_label, safe_chain_slug, summarize_msa_database_path


UI_PACKAGE = "openfoldpanel.UI"
UI_STYLE_FILES = (
    "styles/tokens.css",
    "styles/base.css",
    "styles/layout.css",
    "styles/components.css",
    "styles/figure.css",
    "styles/atmosphere.css",
)
UI_SCRIPT_FILE = "scripts/report.js"
UI_TEMPLATE_FILE = "report.template.html"


def reference_chain_pdf_name(reference_chain: str) -> str:
    """Return the PDF filename for a rendered reference chain."""

    return f"reference-chain-{safe_chain_slug(reference_chain)}.pdf"


def render_reference_chain_report_svg(job_name: str, panel_data: JobPanelData, default_reference_chain: str) -> str:
    """Render one FoldScript-style panel SVG for PDF export."""

    panel_svg, _ = render_panel_svg(panel_data)
    return panel_svg


def render_html_report(report_data: JobReportData, config: PipelineConfig) -> str:
    """Assemble the static HTML report from packaged frontend resources."""

    panel_views = [_build_panel_view(panel_data, config) for panel_data in report_data.chain_panels]
    default_view = next(
        (panel for panel in panel_views if panel["referenceChain"] == report_data.default_reference_chain),
        panel_views[0],
    )
    template = _load_ui_resource(UI_TEMPLATE_FILE)
    replacements = {
        "__OFP_PAGE_TITLE__": html.escape(f"{report_data.job_name} - OpenFoldPanel 图板"),
        "__OFP_REPORT_TITLE__": html.escape(report_data.job_name),
        "__OFP_DEFAULT_CHAIN_ID__": html.escape(report_data.default_reference_chain),
        "__OFP_DEFAULT_CHAIN_LABEL__": html.escape(default_view["chainLabel"]),
        "__OFP_DEFAULT_PANEL_WIDTH__": f'{default_view["panelWidth"]:.2f}',
        "__OFP_LEGEND_CONTACT_STRONG_TEXT__": html.escape(_contact_strong_legend_text(config)),
        "__OFP_LEGEND_CONTACT_WEAK_TEXT__": html.escape(_contact_weak_legend_text(config)),
        "__OFP_THEME_VARS__": _build_theme_vars(report_data.chain_panels[0].render_config),
        "__OFP_INLINE_STYLES__": _load_ui_styles(),
        "__OFP_INLINE_SCRIPT__": _load_ui_script(),
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
                    }
                    for panel in panel_views
                ],
            }
        ),
        "__OFP_CHAIN_TEMPLATES__": "\n".join(panel["templateMarkup"] for panel in panel_views),
    }
    for token, value in replacements.items():
        template = template.replace(token, value)
    return template


def _build_panel_view(panel_data: JobPanelData, config: PipelineConfig) -> dict[str, object]:
    panel_svg, layout = render_panel_svg(panel_data)
    chain_label = humanize_chain_label(panel_data.reference_chain)
    return {
        "referenceChain": panel_data.reference_chain,
        "chainLabel": chain_label,
        "panelWidth": round(layout.width, 2),
        "summaryItems": _build_summary_items(panel_data, config),
        "warnings": list(panel_data.warnings),
        "templateMarkup": _render_chain_template(
            reference_chain=panel_data.reference_chain,
            chain_label=chain_label,
            panel_width=layout.width,
            panel_svg=_strip_xml_declaration(panel_svg),
        ),
    }


def _render_chain_template(*, reference_chain: str, chain_label: str, panel_width: float, panel_svg: str) -> str:
    return (
        f'<template id="ofp-chain-{safe_chain_slug(reference_chain)}" '
        f'data-chain-figure="{html.escape(reference_chain)}" '
        f'data-chain-label="{html.escape(chain_label)}" '
        f'data-panel-width="{panel_width:.2f}">'
        f"{panel_svg}"
        "</template>"
    )


def _build_summary_items(panel_data: JobPanelData, config: PipelineConfig) -> list[dict[str, str]]:
    return [
        {"label": "参考链", "value": humanize_chain_label(panel_data.reference_chain)},
        {"label": "残基范围", "value": _sequence_span_label(panel_data)},
        {"label": "模型数量", "value": str(len(panel_data.models))},
        {"label": "输出状态", "value": _status_label(panel_data.status)},
        {
            "label": "疏水性窗口",
            "value": str(config.hyd_window),
            "tooltip": "Kyte-Doolittle 疏水性计算的滑动窗口大小。对应参数：--hyd-window",
        },
        {
            "label": "显著性阈值",
            "value": config.evalue,
            "tooltip": "命中显著性筛选阈值，值越小越严格。对应参数：--evalue",
        },
        {
            "label": "弱接触阈值",
            "value": f"{config.contact_cutoff:g} A",
            "tooltip": "判定弱接触的距离阈值，单位为埃。对应参数：--contact-cutoff",
        },
        {
            "label": "强接触阈值",
            "value": f"{config.strong_contact_cutoff:g} A",
            "tooltip": "判定强接触的距离阈值，单位为埃。对应参数：--strong-contact-cutoff",
        },
        {
            "label": "同源显示上限",
            "value": str(config.max_homologs_displayed),
            "tooltip": "同源序列检索和展示的最大条数。对应参数：--max-homologs-displayed",
        },
        {
            "label": "数据库",
            "value": summarize_msa_database_path(config.msa_db).upper(),
            "tooltip": "检索使用的数据库名称，取自参数路径末段。对应参数：--msa-db",
        },
    ]


def _sequence_span_label(panel_data: JobPanelData) -> str:
    if not panel_data.sequence_axis:
        return "无可用残基"
    first = panel_data.sequence_axis[0].label
    last = panel_data.sequence_axis[-1].label
    return f"{first} - {last}"


def _contact_strong_legend_text(config: PipelineConfig) -> str:
    return f"按最短非氢原子距离判定：小于 {config.strong_contact_cutoff:g} A。"


def _contact_weak_legend_text(config: PipelineConfig) -> str:
    return (
        f"按最短非氢原子距离判定："
        f"{config.strong_contact_cutoff:g} A 到 {config.contact_cutoff:g} A 之间（含边界）。"
    )


def _status_label(status: str) -> str:
    labels = {
        "success": "成功",
        "partial_success": "部分成功",
        "failed": "失败",
    }
    return labels.get(status, status)


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
            f"  --ofp-color-helix-fill: {colors['helix_fill']};",
            f"  --ofp-color-strand-fill: {colors['strand_fill']};",
            f"  --ofp-color-strand-stroke: {colors['strand_stroke']};",
            f"  --ofp-color-turn: {colors['turn_text']};",
            f"  --ofp-color-msa-query-bg: {colors['msa_query_bg']};",
            f"  --ofp-color-msa-query-text: {colors['msa_query_text']};",
            f"  --ofp-color-msa-identity-bg: {colors['msa_identity_bg']};",
            f"  --ofp-color-msa-identity-text: {colors['msa_identity_text']};",
            f"  --ofp-color-msa-similar-bg: {colors['msa_similar_bg']};",
            f"  --ofp-color-msa-similar-text: {colors['msa_similar_text']};",
            f"  --ofp-color-msa-default-bg: {colors['msa_default_bg']};",
            f"  --ofp-color-msa-default-text: {colors['msa_default_text']};",
            f"  --ofp-color-accessibility-buried: {colors['accessibility_buried']};",
            f"  --ofp-color-accessibility-accessible: {colors['accessibility_accessible']};",
            f"  --ofp-color-hydrophobic: {colors['hydropathy_hydrophobic']};",
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
            f"  --ofp-color-atmosphere: {colors['accent_border']};",
        ]
    )


@lru_cache(maxsize=1)
def _load_ui_styles() -> str:
    return "\n\n".join(_load_ui_resource(path) for path in UI_STYLE_FILES)


@lru_cache(maxsize=1)
def _load_ui_script() -> str:
    return _load_ui_resource(UI_SCRIPT_FILE)


@lru_cache(maxsize=None)
def _load_ui_resource(relative_path: str) -> str:
    return resources.files(UI_PACKAGE).joinpath(relative_path).read_text(encoding="utf-8")
