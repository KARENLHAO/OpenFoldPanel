"""Report-style SVG and HTML renderers for multi-chain outputs."""

from __future__ import annotations

import html
from typing import Iterable

from openfoldpanel.models import JobPanelData, JobReportData
from openfoldpanel.render.svg_renderer import render_panel_svg
from openfoldpanel.utils.text import humanize_chain_label, safe_chain_slug


REPORT_MARGIN = 34
REPORT_CARD_HEIGHT = 82
REPORT_PANEL_PADDING = 18
REPORT_HEADER_HEIGHT = 120
REPORT_LEGEND_HEIGHT = 128


def reference_chain_pdf_name(reference_chain: str) -> str:
    """Return the PDF filename for a rendered reference chain."""

    return f"reference-chain-{safe_chain_slug(reference_chain)}.pdf"


def render_reference_chain_report_svg(job_name: str, panel_data: JobPanelData, default_reference_chain: str) -> str:
    """Render one FoldScript-style panel SVG for PDF export."""

    panel_svg, _ = render_panel_svg(panel_data)
    return panel_svg


def render_html_report(report_data: JobReportData) -> str:
    """Render a minimal static HTML shell around FoldScript-style panels."""

    primary_panel = report_data.chain_panels[0]
    config = primary_panel.render_config
    rendered_sections = [_render_chain_section(panel) for panel in report_data.chain_panels]
    default_panel_width = next(
        (panel_width for _markup, panel_width, chain_id in rendered_sections if chain_id == report_data.default_reference_chain),
        rendered_sections[0][1],
    )
    chain_options = "\n".join(
        _render_chain_option(panel.reference_chain, report_data.default_reference_chain)
        for panel in report_data.chain_panels
    )
    chain_sections = "\n".join(markup for markup, _panel_width, _chain_id in rendered_sections)
    global_notes = _render_global_notes(report_data.warnings)

    return f"""<!DOCTYPE html>
<html lang="zh-CN">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>{html.escape(report_data.job_name)} - OpenFoldPanel 图板</title>
    <style>
      :root {{
        --page-bg: #ffffff;
        --text: {config.colors["text"]};
        --muted: {config.colors["muted_text"]};
        --rule: {config.colors["border"]};
        --accent: {config.colors["strand_fill"]};
        --font: "WenQuanYi Zen Hei", "Noto Sans CJK SC", "PingFang SC", "Microsoft YaHei", {config.font_family};
        --mono: "Liberation Mono", "Nimbus Mono PS", "Courier New", monospace;
        --page-inline-pad: clamp(12px, 2vw, 22px);
        --page-block-pad-top: clamp(22px, 3vw, 26px);
        --page-block-pad-bottom: 40px;
        --panel-readable-width: 760px;
        --rail-width: 296px;
        --layout-gap: 18px;
      }}

      * {{ box-sizing: border-box; }}
      body {{
        margin: 0;
        min-height: 100vh;
        font-family: var(--font);
        color: var(--text);
        background: var(--page-bg);
      }}

      .page {{
        --page-target-width: calc(var(--active-panel-width, {default_panel_width:.2f}px) + var(--rail-width) + var(--layout-gap) + (var(--page-inline-pad) * 2));
        width: min(100%, var(--page-target-width));
        margin: 0 auto;
        padding: var(--page-block-pad-top) var(--page-inline-pad) var(--page-block-pad-bottom);
        transition: width 180ms ease-out, padding 180ms ease-out;
      }}

      .toolbar {{
        display: flex;
        align-items: end;
        justify-content: space-between;
        gap: 18px;
        flex-wrap: wrap;
        margin-bottom: 18px;
        width: 100%;
      }}

      .toolbar-copy {{
        display: grid;
        gap: 6px;
      }}

      .eyebrow {{
        margin: 0;
        color: var(--muted);
        font-size: 0.76rem;
        letter-spacing: 0.18em;
        text-transform: uppercase;
      }}

      .title {{
        margin: 0;
        font-family: var(--mono);
        font-size: 1.36rem;
        font-style: italic;
        color: var(--accent);
      }}

      .current-chain {{
        margin: 0;
        color: var(--muted);
        font-size: 0.95rem;
      }}

      .current-chain strong {{
        color: var(--text);
      }}

      .chain-control {{
        display: grid;
        gap: 8px;
      }}

      .chain-label {{
        margin: 0;
        color: var(--muted);
        font-size: 0.8rem;
      }}

      .chain-select {{
        min-width: 220px;
        border: 1px solid var(--rule);
        background: #fff;
        color: var(--text);
        font: inherit;
        padding: 10px 12px;
      }}

      .chain-stack {{
        display: grid;
        gap: 12px;
        width: 100%;
      }}

      .chain-report {{
        display: none;
        width: 100%;
      }}

      .chain-report.is-active {{
        display: block;
      }}

      .chain-layout {{
        display: flex;
        align-items: flex-start;
        gap: var(--layout-gap);
        width: min(100%, calc(var(--panel-width) + var(--rail-width) + var(--layout-gap)));
        max-width: 100%;
      }}

      .main-stage {{
        flex: 1 1 auto;
        min-width: 0;
        width: min(var(--panel-width), calc(100% - var(--rail-width) - var(--layout-gap)));
      }}

      .supporting-rail {{
        display: grid;
        gap: 14px;
        flex: 0 0 var(--rail-width);
        width: var(--rail-width);
        align-self: start;
        position: sticky;
        top: 18px;
      }}

      .figure-wrap {{
        width: 100%;
        overflow-x: auto;
        overflow-y: hidden;
        border-top: 1px solid #f1f1f1;
        padding-top: 8px;
      }}

      .figure-sheet {{
        width: min(100%, var(--panel-width));
        min-width: min(var(--panel-width), var(--panel-readable-width));
      }}

      .figure-sheet svg {{
        display: block;
        width: 100%;
        height: auto;
      }}

      .rail-card,
      .global-notes {{
        border: 1px solid var(--rule);
        background: #fff;
      }}

      .global-notes {{
        margin-top: 18px;
      }}

      .rail-card summary,
      .global-notes summary,
      .rail-card-title {{
        display: block;
        margin: 0;
        padding: 14px 16px 12px;
        font-size: 0.98rem;
        font-weight: 700;
        color: var(--text);
      }}

      .rail-card summary,
      .global-notes summary {{
        cursor: pointer;
        list-style: none;
      }}

      .rail-card summary::-webkit-details-marker,
      .global-notes summary::-webkit-details-marker {{
        display: none;
      }}

      .rail-card summary::after,
      .global-notes summary::after {{
        content: "+";
        float: right;
        color: var(--muted);
        font-weight: 400;
      }}

      .rail-body {{
        padding: 0 16px 16px;
      }}

      .rail-card[open] summary::after,
      .global-notes[open] summary::after {{
        content: "−";
      }}

      .global-notes ul {{
        margin: 0;
        padding-left: 18px;
        line-height: 1.68;
      }}

      .summary-grid {{
        display: grid;
        gap: 10px;
      }}

      .summary-item {{
        display: grid;
        gap: 4px;
        padding-bottom: 10px;
        border-bottom: 1px solid #efefef;
      }}

      .summary-item:last-child {{
        padding-bottom: 0;
        border-bottom: 0;
      }}

      .summary-label {{
        color: var(--muted);
        font-size: 0.84rem;
      }}

      .summary-value {{
        color: var(--text);
        font-size: 1rem;
        font-weight: 700;
      }}

      .legend-groups {{
        display: grid;
        gap: 14px;
      }}

      .legend-group {{
        display: grid;
        gap: 8px;
      }}

      .legend-group-title {{
        margin: 0;
        color: var(--muted);
        font-size: 0.82rem;
        font-weight: 700;
      }}

      .legend-items {{
        display: grid;
        gap: 7px;
      }}

      .legend-item {{
        display: flex;
        align-items: center;
        gap: 10px;
        font-size: 0.92rem;
      }}

      .legend-swatch {{
        width: 12px;
        height: 12px;
        flex: 0 0 12px;
        border: 1px solid rgba(0, 0, 0, 0.08);
      }}

      @media (max-width: 1449px) {{
        :root {{
          --page-inline-pad: clamp(12px, 1.8vw, 18px);
          --rail-width: 288px;
          --layout-gap: 16px;
        }}
      }}

      @media (max-width: 1100px) {{
        :root {{
          --layout-gap: 14px;
        }}

        .page {{
          width: 100%;
        }}

        .chain-layout {{
          display: grid;
          width: 100%;
          grid-template-columns: 1fr;
        }}

        .main-stage {{
          width: 100%;
        }}

        .supporting-rail {{
          position: static;
          flex: 0 1 auto;
          width: 100%;
          grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
        }}
      }}

      @media (max-width: 720px) {{
        .page {{
          width: 100%;
          padding-inline: 12px;
        }}

        .toolbar {{
          align-items: start;
          gap: 14px;
        }}

        .chain-control {{
          width: 100%;
        }}

        .chain-select {{
          width: 100%;
        }}
      }}
    </style>
  </head>
  <body>
    <main class="page" data-report-page style="--active-panel-width:{default_panel_width:.2f}px;">
      <section class="toolbar" aria-label="图板控制区">
        <div class="toolbar-copy">
          <p class="eyebrow">OpenFoldPanel / FoldScript 风格图板</p>
          <h1 class="title">{html.escape(report_data.job_name)}</h1>
          <p class="current-chain">当前参考链：<strong data-current-chain-label>{html.escape(humanize_chain_label(report_data.default_reference_chain))}</strong></p>
        </div>
        <label class="chain-control">
          <span class="chain-label">参考链选择</span>
          <select class="chain-select" data-chain-select aria-label="参考链选择器">
            {chain_options}
          </select>
        </label>
      </section>

      <section class="chain-stack" aria-label="参考链图板">
        {chain_sections}
      </section>

      {global_notes}
    </main>
    <script>
      const chainSelect = document.querySelector('[data-chain-select]');
      const sections = Array.from(document.querySelectorAll('[data-chain-section]'));
      const currentChainLabel = document.querySelector('[data-current-chain-label]');
      const reportPage = document.querySelector('[data-report-page]');

      function activateChain(chainId, pushHash = true) {{
        if (chainSelect) {{
          chainSelect.value = chainId;
        }}
        let activeSection = null;
        sections.forEach((section) => {{
          const active = section.dataset.chainSection === chainId;
          section.classList.toggle('is-active', active);
          section.toggleAttribute('hidden', !active);
          if (active) {{
            activeSection = section;
          }}
        }});
        if (activeSection && currentChainLabel) {{
          currentChainLabel.textContent = activeSection.dataset.chainLabel;
        }}
        if (activeSection && reportPage && activeSection.dataset.panelWidth) {{
          reportPage.style.setProperty('--active-panel-width', `${{activeSection.dataset.panelWidth}}px`);
        }}
        if (pushHash) {{
          window.history.replaceState(null, '', `#chain-${{chainId}}`);
        }}
      }}

      if (chainSelect) {{
        chainSelect.addEventListener('change', (event) => activateChain(event.target.value));
      }}

      const requestedChain = window.location.hash.replace('#chain-', '');
      const availableChains = sections.map((section) => section.dataset.chainSection);
      const initialChain = availableChains.includes(requestedChain)
        ? requestedChain
        : '{html.escape(report_data.default_reference_chain)}';
      activateChain(initialChain, false);
    </script>
  </body>
</html>
"""


def _render_chain_option(reference_chain: str, default_reference_chain: str) -> str:
    label = humanize_chain_label(reference_chain)
    selected = " selected" if reference_chain == default_reference_chain else ""
    return f'<option value="{html.escape(reference_chain)}"{selected}>{html.escape(label)}</option>'


def _render_chain_section(panel_data: JobPanelData) -> tuple[str, float, str]:
    panel_svg, layout = render_panel_svg(panel_data)
    report_svg = _strip_xml_declaration(panel_svg)
    summary = _render_chain_summary(panel_data)
    legend = _render_html_legend(panel_data.render_config.colors)
    warnings = _render_chain_warnings(panel_data)
    panel_width = layout.width
    markup = (
        f'<article class="chain-report" id="chain-{html.escape(panel_data.reference_chain)}" '
        f'data-chain-section="{html.escape(panel_data.reference_chain)}" '
        f'data-chain-label="{html.escape(humanize_chain_label(panel_data.reference_chain))}" '
        f'data-panel-width="{panel_width:.2f}" '
        f'style="--panel-width:{panel_width:.2f}px;" '
        'hidden>'
        '<div class="chain-layout">'
        '<section class="main-stage">'
        f'<div class="figure-wrap"><div class="figure-sheet">{report_svg}</div></div>'
        '</section>'
        '<aside class="supporting-rail" aria-label="链摘要与图例">'
        f"{summary}"
        f"{legend}"
        f"{warnings}"
        '</aside>'
        '</div>'
        "</article>"
    )
    return markup, panel_width, panel_data.reference_chain


def _msa_label(enabled: bool) -> str:
    return "完整比对" if enabled else "仅查询序列"


def _render_stage_metrics(panel_data: JobPanelData) -> str:
    metrics = [
        ("参考链", humanize_chain_label(panel_data.reference_chain)),
        ("序列长度", str(len(panel_data.sequence_axis))),
        ("模型数量", str(len(panel_data.models))),
        ("序列比对", _msa_label(panel_data.msa.enabled)),
    ]
    items = "".join(
        (
            '<div class="metric">'
            f'<p class="metric-label">{html.escape(label)}</p>'
            f'<p class="metric-value">{html.escape(value)}</p>'
            "</div>"
        )
        for label, value in metrics
    )
    return f'<section class="stage-meta" aria-label="当前链概览">{items}</section>'


def _render_chain_summary(panel_data: JobPanelData) -> str:
    items = [
        ("参考链", humanize_chain_label(panel_data.reference_chain)),
        ("残基范围", _sequence_span_label(panel_data)),
        ("模型数量", str(len(panel_data.models))),
        ("输出状态", _status_label(panel_data.status)),
    ]
    body = "".join(
        (
            '<div class="summary-item">'
            f'<span class="summary-label">{html.escape(label)}</span>'
            f'<span class="summary-value">{html.escape(value)}</span>'
            "</div>"
        )
        for label, value in items
    )
    return (
        '<details class="rail-card" open>'
        '<summary>链摘要</summary>'
        f'<div class="rail-body"><div class="summary-grid">{body}</div></div>'
        '</details>'
    )


def _render_html_legend(colors: dict[str, str]) -> str:
    groups = [
        (
            "二级结构与置信度",
            [
                ("螺旋", colors["helix_fill"]),
                ("折叠片", colors["strand_fill"]),
                ("转角", colors["turn_text"]),
                ("高置信度", colors["plddt_very_high"]),
            ],
        ),
        (
            "表面与理化性质",
            [
                ("埋藏", colors["accessibility_buried"]),
                ("暴露", colors["accessibility_accessible"]),
                ("疏水", colors["hydropathy_hydrophobic"]),
                ("亲水", colors["hydropathy_hydrophilic"]),
            ],
        ),
        (
            "接触信息",
            [
                ("强接触", colors["contact_strong"]),
                ("弱接触", colors["contact_weak"]),
                ("多重接触", colors["contact_multi_outline"]),
                ("离子或配体", colors["warning"]),
            ],
        ),
    ]
    group_markup = []
    for title, items in groups:
        item_markup = "".join(
            (
                '<div class="legend-item">'
                f'<span class="legend-swatch" style="background:{html.escape(color)}"></span>'
                f'<span>{html.escape(label)}</span>'
                "</div>"
            )
            for label, color in items
        )
        group_markup.append(
            (
                '<section class="legend-group">'
                f'<h3 class="legend-group-title">{html.escape(title)}</h3>'
                f'<div class="legend-items">{item_markup}</div>'
                "</section>"
            )
        )
    return (
        '<section class="rail-card rail-card-static" aria-label="图例">'
        '<h2 class="rail-card-title">图例</h2>'
        f'<div class="rail-body"><div class="legend-groups">{"".join(group_markup)}</div></div>'
        '</section>'
    )


def _render_chain_warnings(panel_data: JobPanelData) -> str:
    if not panel_data.warnings:
        return ""
    items = "".join(f"<li>{html.escape(note)}</li>" for note in panel_data.warnings)
    return (
        '<details class="global-notes">'
        f'<summary>当前链提示（{len(panel_data.warnings)}）</summary>'
        f'<div class="rail-body"><ul>{items}</ul></div>'
        '</details>'
    )


def _render_global_notes(warnings: Iterable[str]) -> str:
    warning_list = list(warnings)
    if not warning_list:
        return ""
    items = "".join(f"<li>{html.escape(note)}</li>" for note in warning_list)
    return (
        '<details class="global-notes" open>'
        '<summary>任务级提示</summary>'
        f'<div class="rail-body"><ul>{items}</ul></div>'
        '</details>'
    )


def _sequence_span_label(panel_data: JobPanelData) -> str:
    if not panel_data.sequence_axis:
        return "无可用残基"
    first = panel_data.sequence_axis[0].label
    last = panel_data.sequence_axis[-1].label
    return f"{first} - {last}"


def _status_label(status: str) -> str:
    labels = {
        "success": "成功",
        "partial_success": "部分成功",
        "failed": "失败",
    }
    return labels.get(status, status)


def _render_metric_cards(panel_data: JobPanelData, report_width: float, y: float, is_default_reference: bool) -> str:
    width = report_width - REPORT_MARGIN * 2
    card_gap = 14
    card_width = (width - card_gap * 3) / 4
    x_positions = [REPORT_MARGIN + index * (card_width + card_gap) for index in range(4)]
    cards = [
        ("参考链", humanize_chain_label(panel_data.reference_chain)),
        ("序列长度", str(len(panel_data.sequence_axis))),
        ("模型数量", str(len(panel_data.models))),
        ("序列比对", _msa_label(panel_data.msa.enabled)),
    ]
    pieces = []
    for index, (label, value) in enumerate(cards):
        x = x_positions[index]
        pieces.append(
            f'<g transform="translate({x:.2f},{y:.2f})">'
            f'<rect class="report-card" width="{card_width:.2f}" height="{REPORT_CARD_HEIGHT:.2f}" rx="0"/>'
            f'<text class="report-card-label" x="14" y="22">{html.escape(label)}</text>'
            f'<text class="report-card-value" x="14" y="56">{html.escape(value)}</text>'
            "</g>"
        )
    badge_label = "默认展示链" if is_default_reference else "附加参考链"
    pieces.append(
        f'<text class="report-side-note" x="{report_width - REPORT_MARGIN:.2f}" y="{y - 10:.2f}" text-anchor="end">{html.escape(badge_label)}</text>'
    )
    return "\n".join(pieces)


def _render_legend(report_width: float, y: float, colors: dict[str, str]) -> str:
    x = REPORT_MARGIN
    width = report_width - REPORT_MARGIN * 2
    columns = [
        (
            "轨道图例",
            [
                ("螺旋", colors["helix_fill"]),
                ("折叠片", colors["strand_fill"]),
                ("转角", colors["turn_text"]),
                ("高置信度", colors["plddt_very_high"]),
            ],
        ),
        (
            "表面与理化性质",
            [
                ("埋藏", colors["accessibility_buried"]),
                ("暴露", colors["accessibility_accessible"]),
                ("疏水", colors["hydropathy_hydrophobic"]),
                ("亲水", colors["hydropathy_hydrophilic"]),
            ],
        ),
        (
            "接触信息",
            [
                ("强接触", colors["contact_strong"]),
                ("弱接触", colors["contact_weak"]),
                ("多重接触", colors["contact_multi_outline"]),
                ("离子或配体", colors["warning"]),
            ],
        ),
    ]
    column_gap = 18
    column_width = (width - column_gap * 2) / 3
    pieces = [f'<g transform="translate({x:.2f},{y:.2f})"><rect class="legend-card" width="{width:.2f}" height="{REPORT_LEGEND_HEIGHT:.2f}" rx="0"/>']
    for index, (title, items) in enumerate(columns):
        column_x = 16 + index * (column_width + column_gap)
        pieces.append(f'<text class="legend-title" x="{column_x:.2f}" y="26">{html.escape(title)}</text>')
        for row_index, (label, color) in enumerate(items):
            item_y = 48 + row_index * 18
            pieces.append(
                f'<rect x="{column_x:.2f}" y="{item_y - 9:.2f}" width="12" height="12" fill="{color}" stroke="{color}"/>'
                f'<text class="legend-label" x="{column_x + 20:.2f}" y="{item_y:.2f}">{html.escape(label)}</text>'
            )
    pieces.append("</g>")
    return "\n".join(pieces)


def _render_panel_card(panel_svg: str, report_width: float, panel_height: float, y: float) -> str:
    outer_width = report_width - REPORT_MARGIN * 2
    panel_x = REPORT_MARGIN + REPORT_PANEL_PADDING
    panel_y = y + REPORT_PANEL_PADDING
    nested_svg = _position_nested_svg(panel_svg, x=panel_x, y=panel_y)
    return (
        f'<g transform="translate({REPORT_MARGIN:.2f},{y:.2f})">'
        f'<rect class="panel-shadow" x="8" y="10" width="{outer_width:.2f}" height="{panel_height + REPORT_PANEL_PADDING * 2:.2f}" rx="0"/>'
        f'<rect class="panel-card" width="{outer_width:.2f}" height="{panel_height + REPORT_PANEL_PADDING * 2:.2f}" rx="0"/>'
        "</g>\n"
        f"{nested_svg}"
    )


def _render_warning_box(
    warnings: Iterable[str],
    *,
    x: float,
    y: float,
    width: float,
    colors: dict[str, str],
    font_family: str,
) -> str:
    warning_list = list(warnings)
    if not warning_list:
        return ""
    height = _warning_box_height(warning_list)
    pieces = [
        f'<g transform="translate({x:.2f},{y:.2f})">',
        f'<rect width="{width:.2f}" height="{height:.2f}" fill="{colors["warning_bg"]}" stroke="{colors["warning_border"]}"/>',
        f"<text x='14' y='24' font-family='{font_family}' font-size='11' fill='{colors['warning']}' letter-spacing='1.3'>链级提示</text>",
    ]
    for index, warning in enumerate(warning_list):
        item_y = 46 + index * 17
        pieces.append(
            (
                f"<text x='14' y='{item_y:.2f}' font-family='{font_family}' font-size='12' "
                f"fill='{colors['warning']}'>{html.escape('- ' + warning)}</text>"
            )
        )
    pieces.append("</g>")
    return "\n".join(pieces)


def _warning_box_height(warnings: Iterable[str]) -> float:
    warning_list = list(warnings)
    if not warning_list:
        return 0
    return 54 + max(0, len(warning_list) - 1) * 17


def _position_nested_svg(svg_markup: str, *, x: float, y: float) -> str:
    nested_svg = _strip_xml_declaration(svg_markup)
    return nested_svg.replace("<svg ", f'<svg x="{x:.2f}" y="{y:.2f}" ', 1)


def _strip_xml_declaration(svg_markup: str) -> str:
    if svg_markup.startswith("<?xml"):
        return svg_markup.split("\n", 1)[1]
    return svg_markup


def _render_background_texture(width: float, height: float, colors: dict[str, str]) -> str:
    return (
        '<g opacity="0.65">'
        f'<line x1="{REPORT_MARGIN:.2f}" y1="24" x2="{width - REPORT_MARGIN:.2f}" y2="24" stroke="{colors["accent"]}" stroke-width="1.5"/>'
        f'<line x1="{REPORT_MARGIN:.2f}" y1="132" x2="{width - REPORT_MARGIN:.2f}" y2="132" stroke="{colors["border"]}" stroke-width="1"/>'
        f'<line x1="{width - REPORT_MARGIN - 220:.2f}" y1="48" x2="{width - REPORT_MARGIN:.2f}" y2="48" stroke="{colors["accent_border"]}" stroke-width="2"/>'
        f'<circle cx="{width - REPORT_MARGIN - 22:.2f}" cy="48" r="6" fill="{colors["accent"]}"/>'
        '</g>'
    )


def _report_style(font_family: str, heading_font_family: str, colors: dict[str, str]) -> str:
    return (
        "<style>"
        f'.report-title{{font-family:{heading_font_family};font-size:44px;font-weight:700;fill:{colors["text"]};letter-spacing:-1.2px;}}'
        f'.report-subtitle{{font-family:{font_family};font-size:18px;fill:{colors["muted_text"]};}}'
        f'.report-kicker{{font-family:{font_family};font-size:11px;fill:{colors["muted_text"]};letter-spacing:2.6px;text-transform:uppercase;}}'
        f'.report-card{{fill:{colors["surface"]};stroke:{colors["border"]};stroke-width:1;}}'
        f'.report-card-label{{font-family:{font_family};font-size:11px;fill:{colors["muted_text"]};letter-spacing:1.4px;text-transform:uppercase;}}'
        f'.report-card-value{{font-family:{heading_font_family};font-size:28px;fill:{colors["text"]};}}'
        f'.report-side-note{{font-family:{font_family};font-size:12px;fill:{colors["muted_text"]};}}'
        f'.legend-card{{fill:{colors["surface"]};stroke:{colors["border"]};stroke-width:1;}}'
        f'.legend-title{{font-family:{font_family};font-size:12px;fill:{colors["text"]};font-weight:700;}}'
        f'.legend-label{{font-family:{font_family};font-size:11px;fill:{colors["muted_text"]};}}'
        f'.panel-card{{fill:{colors["surface"]};stroke:{colors["border"]};stroke-width:1;}}'
        f'.panel-shadow{{fill:{colors["accent_soft"]};opacity:0.55;}}'
        "</style>"
    )
