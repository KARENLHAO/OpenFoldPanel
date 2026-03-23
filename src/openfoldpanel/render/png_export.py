"""PNG export from SVG using CairoSVG when available."""

from __future__ import annotations

from pathlib import Path


def export_png(svg_path: Path, output_path: Path) -> tuple[bool, str | None]:
    """Convert an SVG file to PNG."""

    try:
        import cairosvg  # type: ignore
    except ImportError:
        return False, "CairoSVG is not installed; PNG export was skipped."
    cairosvg.svg2png(url=str(svg_path), write_to=str(output_path))
    return True, None
