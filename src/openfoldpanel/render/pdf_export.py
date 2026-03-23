"""PDF export from SVG markup using CairoSVG when available."""

from __future__ import annotations

from pathlib import Path


def export_pdf(svg_markup: str, output_path: Path) -> tuple[bool, str | None]:
    """Convert SVG markup to PDF."""

    try:
        import cairosvg  # type: ignore
    except ImportError:
        return False, "CairoSVG is not installed; PDF export was skipped."
    cairosvg.svg2pdf(bytestring=svg_markup.encode("utf-8"), write_to=str(output_path))
    return True, None
