"""Command-line interface for openfoldpanel."""

from __future__ import annotations

import argparse
from pathlib import Path

from openfoldpanel.constants import (
    ALLOWED_EVALUES,
    DEFAULT_COLUMNS,
    DEFAULT_CONTACT_CUTOFF,
    DEFAULT_EVALUE,
    DEFAULT_FONT_SIZE,
    DEFAULT_HYDROPATHY_WINDOW,
    DEFAULT_MAX_HOMOLOGS_DISPLAYED,
    DEFAULT_STRONG_CONTACT_CUTOFF,
    MAX_HOMOLOGS_DISPLAYED_LIMIT,
    validate_evalue,
)
from openfoldpanel.logging_utils import configure_root_logger
from openfoldpanel.models import PipelineConfig
from openfoldpanel.pipeline import run_pipeline


def _parse_max_homologs_displayed(value: str) -> int:
    """Parse and validate the rendered/search homolog limit."""

    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--max-homologs-displayed must be an integer.") from exc
    if not 0 <= parsed <= MAX_HOMOLOGS_DISPLAYED_LIMIT:
        raise argparse.ArgumentTypeError(
            f"--max-homologs-displayed must be between 0 and {MAX_HOMOLOGS_DISPLAYED_LIMIT}."
        )
    return parsed


def _parse_evalue(value: str) -> str:
    """Parse and validate the fixed BLAST/MMseqs significance threshold."""

    try:
        return validate_evalue(value, parameter_name="--evalue")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser."""

    parser = argparse.ArgumentParser(description="Render local FoldScript-like flat representation panels.")
    parser.add_argument("--input", required=True, type=Path, help="Input structure file or archive.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--chain", default="ALL", help="Reference chain selector: ALL or a specific chain ID.")
    parser.add_argument("--columns", type=int, default=DEFAULT_COLUMNS, help="Residue columns per render block.")
    parser.add_argument("--font-size", type=int, default=DEFAULT_FONT_SIZE, help="Base font size.")
    parser.add_argument("--hyd-window", type=int, default=DEFAULT_HYDROPATHY_WINDOW, help="Hydropathy smoothing window.")
    parser.add_argument(
        "--msa-db",
        type=Path,
        default=None,
        help="Local BLAST/MMseqs database prefix, or a protein FASTA such as uniprot_sprot.fasta.",
    )
    parser.add_argument(
        "--max-homologs-displayed",
        type=_parse_max_homologs_displayed,
        default=DEFAULT_MAX_HOMOLOGS_DISPLAYED,
        help="Maximum homolog sequences to search and render (0-25).",
    )
    parser.add_argument(
        "--evalue",
        type=_parse_evalue,
        default=DEFAULT_EVALUE,
        help=(
            "BLAST/MMseqs hit significance threshold "
            f"({', '.join(ALLOWED_EVALUES)}; default: {DEFAULT_EVALUE})."
        ),
    )
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
        max_homologs_displayed=args.max_homologs_displayed,
        evalue=args.evalue,
        font_size=args.font_size,
        hyd_window=args.hyd_window,
        msa_db=args.msa_db,
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
