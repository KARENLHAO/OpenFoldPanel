"""Pipeline orchestration for openfoldpanel."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

from openfoldpanel.extractors.archive import extract_archive
from openfoldpanel.extractors.discovery import discover_jobs_from_extracted_root, discover_jobs_from_structure
from openfoldpanel.extractors.validators import is_supported_archive, validate_input_path
from openfoldpanel.features.accessibility import build_accessibility_track
from openfoldpanel.features.conservation import compute_conservation
from openfoldpanel.features.contacts import compute_contacts
from openfoldpanel.features.disulfide import infer_disulfides
from openfoldpanel.features.dssp_runner import run_dssp
from openfoldpanel.features.hydropathy import compute_hydropathy
from openfoldpanel.features.msa_align import align_sequences
from openfoldpanel.features.msa_search import search_homologs
from openfoldpanel.features.secondary_structure import build_secondary_structure_track
from openfoldpanel.io.json_dump import write_tracks_json
from openfoldpanel.io.writers import write_summary
from openfoldpanel.logging_utils import attach_file_logger, detach_handler
from openfoldpanel.models import JobPanelData, JobReportData, JobRunResult, MSAData, MSARow, ModelTracks, PipelineConfig
from openfoldpanel.parsers.plddt_reader import residue_plddt
from openfoldpanel.parsers.sequence_mapper import align_chain_to_axis, build_sequence_axis
from openfoldpanel.parsers.structure_parser import (
    collect_reference_chains,
    get_chain_or_best_match,
    parse_structure,
    select_reference_chain,
)
from openfoldpanel.render.layout import build_render_config
from openfoldpanel.render.pdf_export import export_pdf
from openfoldpanel.render.report_renderer import reference_chain_pdf_name, render_html_report, render_reference_chain_report_svg
from openfoldpanel.utils.filesystem import ensure_directory, safe_rmtree
from openfoldpanel.utils.text import humanize_model_name, safe_chain_slug


def run_pipeline(config: PipelineConfig, logger: logging.Logger) -> dict[str, int]:
    """Run the full batch pipeline and return a status summary."""

    validate_input_path(config.input_path)
    ensure_directory(config.outdir)

    temp_dir = Path(tempfile.mkdtemp(prefix="openfoldpanel_"))
    try:
        if is_supported_archive(config.input_path):
            extracted_root = temp_dir / "archive"
            extract_archive(config.input_path, extracted_root)
            jobs = discover_jobs_from_extracted_root(extracted_root, logger)
        else:
            jobs = discover_jobs_from_structure(config.input_path)

        results: list[JobRunResult] = []
        for job in jobs:
            results.append(_run_job(job, config, logger, temp_dir))

        return {
            "total_jobs": len(results),
            "success": sum(result.status == "success" for result in results),
            "partial_success": sum(result.status == "partial_success" for result in results),
            "failed": sum(result.status == "failed" for result in results),
        }
    finally:
        if not config.keep_temp:
            safe_rmtree(temp_dir)


def _run_job(job, config: PipelineConfig, logger: logging.Logger, temp_dir: Path) -> JobRunResult:
    output_dir = ensure_directory(config.outdir / job.name)
    log_path = output_dir / "logs.txt"
    file_handler = attach_file_logger(logger, log_path)
    warnings: list[str] = list(job.ignored_files)
    partial_reasons: list[str] = []
    try:
        if not job.structure_files:
            result = JobRunResult(job_name=job.name, status="failed", output_dir=str(output_dir), error="No structure files found.")
            write_summary(result, None, output_dir / "summary.txt")
            return result

        parsed_structures = []
        for structure_path in job.structure_files:
            logger.info("Parsing %s for job %s", structure_path.name, job.name)
            try:
                parsed_structures.append(parse_structure(structure_path, logger))
            except Exception as exc:
                warning = f"Failed to parse {structure_path.name}: {exc}"
                logger.warning(warning)
                warnings.append(warning)
                partial_reasons.append(warning)

        if not parsed_structures:
            result = JobRunResult(
                job_name=job.name,
                status="failed",
                output_dir=str(output_dir),
                warnings=warnings,
                error="All models failed to parse.",
            )
            write_summary(result, None, output_dir / "summary.txt")
            return result

        reference_structure = parsed_structures[0]
        default_reference_chain = select_reference_chain(reference_structure, "AUTO")
        requested_reference_chains = collect_reference_chains(reference_structure, config.chain)
        render_config = build_render_config(config.columns, config.font_size, config.msa_display_rows)

        chain_panels: list[JobPanelData] = []
        for reference_chain_id in requested_reference_chains:
            panel_data, chain_partial_reasons = _build_panel_data_for_reference_chain(
                parsed_structures=parsed_structures,
                reference_chain_id=reference_chain_id,
                config=config,
                render_config=render_config,
                workdir=temp_dir / job.name / safe_chain_slug(reference_chain_id),
                logger=logger,
            )
            if panel_data is None:
                warning = f"Skipped {reference_chain_id}: no compatible models could be mapped to the reference axis."
                logger.warning(warning)
                warnings.append(warning)
                partial_reasons.append(warning)
                continue

            chain_panels.append(panel_data)
            prefixed_warnings = [f"Chain {reference_chain_id}: {note}" for note in panel_data.warnings]
            warnings.extend(prefixed_warnings)
            partial_reasons.extend(f"Chain {reference_chain_id}: {note}" for note in chain_partial_reasons)

        if not chain_panels:
            result = JobRunResult(
                job_name=job.name,
                status="failed",
                output_dir=str(output_dir),
                warnings=warnings,
                error="No compatible models could be mapped to any protein reference chain.",
            )
            write_summary(result, None, output_dir / "summary.txt")
            return result

        if default_reference_chain not in {panel.reference_chain for panel in chain_panels}:
            default_reference_chain = chain_panels[0].reference_chain

        report_data = JobReportData(
            job_name=job.name,
            default_reference_chain=default_reference_chain,
            chain_panels=chain_panels,
            warnings=warnings,
            status="success",
        )

        artifacts = []
        for panel_data in chain_panels:
            pdf_path = output_dir / reference_chain_pdf_name(panel_data.reference_chain)
            report_svg = render_reference_chain_report_svg(job.name, panel_data, default_reference_chain)
            pdf_ok, pdf_warning = export_pdf(report_svg, pdf_path)
            if pdf_ok:
                artifacts.append(pdf_path.name)
            elif pdf_warning:
                warnings.append(pdf_warning)
                partial_reasons.append(pdf_warning)

        job_status = "partial_success" if partial_reasons else "success"
        report_data.status = job_status

        html_path = output_dir / "report.html"
        html_path.write_text(render_html_report(report_data), encoding="utf-8")
        artifacts.append(html_path.name)

        json_path = output_dir / "tracks.json"
        write_tracks_json(report_data, json_path)
        artifacts.append(json_path.name)
        artifacts.append("summary.txt")

        result = JobRunResult(
            job_name=job.name,
            status=job_status,
            output_dir=str(output_dir),
            warnings=warnings,
            artifacts=artifacts,
        )
        write_summary(result, report_data, output_dir / "summary.txt")
        return result
    except Exception as exc:
        logger.exception("Job %s failed", job.name)
        result = JobRunResult(
            job_name=job.name,
            status="failed",
            output_dir=str(output_dir),
            warnings=warnings,
            error=str(exc),
        )
        write_summary(result, None, output_dir / "summary.txt")
        return result
    finally:
        detach_handler(logger, file_handler)


def _build_panel_data_for_reference_chain(
    *,
    parsed_structures: list,
    reference_chain_id: str,
    config: PipelineConfig,
    render_config,
    workdir: Path,
    logger: logging.Logger,
) -> tuple[JobPanelData | None, list[str]]:
    reference_structure = parsed_structures[0]
    reference_chain = reference_structure.chains[reference_chain_id]
    axis = build_sequence_axis(reference_chain)
    hydropathy = compute_hydropathy(axis, config.hyd_window)

    warnings: list[str] = []
    partial_reasons: list[str] = []
    model_tracks: list[ModelTracks] = []
    for structure in parsed_structures:
        chain = get_chain_or_best_match(structure, reference_chain_id, reference_chain.sequence)
        if chain is None:
            warning = f"Skipping {structure.source_path.name}: no compatible protein chain found."
            logger.warning(warning)
            warnings.append(warning)
            partial_reasons.append(warning)
            continue

        alignment = align_chain_to_axis(axis, chain)
        warnings.extend(alignment.warnings)

        dssp_features, dssp_warnings = run_dssp(structure.source_path, logger)
        warnings.extend(dssp_warnings)

        secondary_structure = build_secondary_structure_track(axis, alignment.residue_by_axis_index, dssp_features)
        accessibility = build_accessibility_track(axis, alignment.residue_by_axis_index, dssp_features)
        contacts = compute_contacts(
            structure,
            chain,
            alignment.residue_by_axis_index,
            axis,
            cutoff=config.contact_cutoff,
            strong_cutoff=config.strong_contact_cutoff,
        )
        plddt = [
            residue_plddt(alignment.residue_by_axis_index[position.residue_index])
            if position.residue_index in alignment.residue_by_axis_index
            else None
            for position in axis
        ]
        disulfides = infer_disulfides(alignment.residue_by_axis_index)
        model_tracks.append(
            ModelTracks(
                name=f"{structure.name}_{chain.chain_id}",
                source_path=str(structure.source_path),
                chain=chain.chain_id,
                secondary_structure=secondary_structure,
                plddt=plddt,
                accessibility=accessibility,
                contacts=contacts,
                disulfides=disulfides,
                display_name=humanize_model_name(structure.name, chain.chain_id),
            )
        )

    if not model_tracks:
        return None, partial_reasons

    msa = _build_msa_data(
        axis=axis,
        config=config,
        workdir=workdir,
        logger=logger,
        warnings=warnings,
    )
    panel_status = "partial_success" if partial_reasons else "success"
    panel_data = JobPanelData(
        job_name=reference_structure.name,
        reference_chain=reference_chain_id,
        sequence_axis=axis,
        models=model_tracks,
        msa=msa,
        hydropathy=hydropathy,
        render_config=render_config,
        warnings=warnings,
        status=panel_status,
    )
    return panel_data, partial_reasons


def _build_msa_data(*, axis, config: PipelineConfig, workdir: Path, logger: logging.Logger, warnings: list[str]) -> MSAData:
    ensure_directory(workdir)
    query_sequence = "".join(position.one_letter for position in axis)
    if config.disable_msa:
        return MSAData(
            enabled=False,
            query=query_sequence,
            rows=[MSARow(identifier="Query Sequence", sequence=query_sequence, is_query=True)],
            warnings=["MSA disabled by user."],
        )

    if config.msa_db is None:
        return MSAData(
            enabled=False,
            query=query_sequence,
            rows=[MSARow(identifier="Query Sequence", sequence=query_sequence, is_query=True)],
            warnings=["No MSA database provided; sequence alignment was skipped."],
        )

    query_fasta = workdir / "query.fasta"
    query_fasta.write_text(f">query\n{query_sequence}\n", encoding="utf-8")
    hits, hit_warnings = search_homologs(
        query_fasta,
        config.msa_db,
        max_hits=config.max_hits,
        workdir=workdir,
        logger=logger,
    )
    if hit_warnings:
        warnings.extend(hit_warnings)
    if not hits:
        return MSAData(
            enabled=False,
            query=query_sequence,
            rows=[MSARow(identifier="Query Sequence", sequence=query_sequence, is_query=True)],
            warnings=hit_warnings,
        )

    alignment_rows, align_warnings = align_sequences([("Query Sequence", query_sequence), *hits], workdir / "alignment.fasta", logger)
    if align_warnings:
        warnings.extend(align_warnings)
    if not alignment_rows:
        return MSAData(
            enabled=False,
            query=query_sequence,
            rows=[MSARow(identifier="Query Sequence", sequence=query_sequence, is_query=True)],
            warnings=hit_warnings + align_warnings,
        )

    projected_rows = _project_alignment_to_query_axis(alignment_rows)
    return MSAData(
        enabled=True,
        query=query_sequence,
        rows=projected_rows,
        conservation=compute_conservation(projected_rows),
        warnings=hit_warnings + align_warnings,
    )


def _project_alignment_to_query_axis(rows: list[MSARow]) -> list[MSARow]:
    query = next((row for row in rows if row.is_query), rows[0])
    keep_indices = [index for index, residue in enumerate(query.sequence) if residue != "-"]
    projected: list[MSARow] = []
    for row in rows:
        sequence = "".join(row.sequence[index] for index in keep_indices if index < len(row.sequence))
        projected.append(MSARow(identifier=row.identifier, sequence=sequence, is_query=row.is_query))
    return projected
