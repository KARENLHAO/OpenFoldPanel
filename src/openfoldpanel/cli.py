"""Command-line interface for openfoldpanel."""

from __future__ import annotations

import argparse
from pathlib import Path

from openfoldpanel.constants import (
    DEFAULT_COLUMNS,
    DEFAULT_CONTACT_CUTOFF,
    DEFAULT_FONT_SIZE,
    DEFAULT_HYDROPATHY_WINDOW,
    DEFAULT_MAX_HITS,
    DEFAULT_MSA_DISPLAY_ROWS,
    DEFAULT_STRONG_CONTACT_CUTOFF,
)
from openfoldpanel.logging_utils import configure_root_logger
from openfoldpanel.models import PipelineConfig
from openfoldpanel.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser."""

    parser = argparse.ArgumentParser(description="Render local FoldScript-like flat representation panels.")
    parser.add_argument("--input", required=True, type=Path, help="Input structure file or archive.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--chain", default="ALL", help="Reference chain selector: ALL or a specific chain ID.")
    parser.add_argument("--columns", type=int, default=DEFAULT_COLUMNS, help="Residue columns per render block.")
    parser.add_argument(
        "--msa-display-rows",
        type=int,
        default=DEFAULT_MSA_DISPLAY_ROWS,
        help="Maximum homolog sequence rows to show in the rendered panel.",
    )
    parser.add_argument("--font-size", type=int, default=DEFAULT_FONT_SIZE, help="Base font size.")
    parser.add_argument("--hyd-window", type=int, default=DEFAULT_HYDROPATHY_WINDOW, help="Hydropathy smoothing window.")
    parser.add_argument("--msa-db", type=Path, default=None, help="Local BLAST/MMseqs database path.")
    parser.add_argument("--max-hits", type=int, default=DEFAULT_MAX_HITS, help="Maximum homolog hits for MSA.")
    parser.add_argument("--disable-msa", action="store_true", help="Disable MSA search and alignment.")
    parser.add_argument("--keep-temp", action="store_true", help="Keep extracted temporary files.")
    parser.add_argument("--contact-cutoff", type=float, default=DEFAULT_CONTACT_CUTOFF, help="Weak contact cutoff in Angstrom.")
    parser.add_argument(
        "--strong-contact-cutoff",
        type=float,
        default=DEFAULT_STRONG_CONTACT_CUTOFF,
        help="Strong contact cutoff in Angstrom.",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""

    parser = build_parser()
    args = parser.parse_args(argv)
    if args.chain.upper() == "AUTO":
        parser.error("--chain AUTO has been deprecated. Use the default ALL mode or provide an explicit chain ID.")
    logger = configure_root_logger(args.verbose)
    config = PipelineConfig(
        input_path=args.input,
        outdir=args.outdir,
        chain=args.chain,
        columns=args.columns,
        msa_display_rows=max(0, args.msa_display_rows),
        font_size=args.font_size,
        hyd_window=args.hyd_window,
        msa_db=args.msa_db,
        max_hits=args.max_hits,
        disable_msa=args.disable_msa,
        keep_temp=args.keep_temp,
        contact_cutoff=args.contact_cutoff,
        strong_contact_cutoff=args.strong_contact_cutoff,
        verbose=args.verbose,
    )
    summary = run_pipeline(config, logger)
    logger.info(
        "Jobs: total=%s success=%s partial_success=%s failed=%s",
        summary["total_jobs"],
        summary["success"],
        summary["partial_success"],
        summary["failed"],
    )
    return 0 if summary["failed"] == 0 else 1
