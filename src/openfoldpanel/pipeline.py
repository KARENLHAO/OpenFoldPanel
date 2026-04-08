"""Pipeline orchestration for openfoldpanel."""

from __future__ import annotations

from dataclasses import dataclass, field
import logging
import tempfile
from collections.abc import Iterable
from pathlib import Path

from openfoldpanel.extractors.archive import extract_archive
from openfoldpanel.extractors.discovery import discover_jobs_from_extracted_root, discover_jobs_from_structure
from openfoldpanel.extractors.validators import is_supported_archive, validate_input_path
from openfoldpanel.features.accessibility import build_accessibility_track
from openfoldpanel.features.antibody_annotation import annotate_antibody_chain
from openfoldpanel.features.batch_analysis import build_batch_analysis
from openfoldpanel.features.conservation import compute_conservation
from openfoldpanel.features.contacts import compute_contacts
from openfoldpanel.features.disulfide import infer_disulfides
from openfoldpanel.features.dssp_runner import run_dssp
from openfoldpanel.features.hydropathy import compute_hydropathy
from openfoldpanel.features.msa_align import align_sequences
from openfoldpanel.features.msa_search import search_homologs
from openfoldpanel.features.secondary_structure import build_secondary_structure_track
from openfoldpanel.io.csv_stats import write_statistics_csvs
from openfoldpanel.io.json_dump import write_tracks_json
from openfoldpanel.io.writers import write_summary
from openfoldpanel.logging_utils import attach_file_logger, detach_handler
from openfoldpanel.models import (
    JobDefinition,
    JobPanelData,
    JobReportData,
    JobRunResult,
    MSAData,
    MSARow,
    ModelTracks,
    ParsedStructure,
    PipelineConfig,
    RenderConfig,
    SequenceAxisPosition,
)
from openfoldpanel.parsers.plddt_reader import residue_plddt
from openfoldpanel.parsers.sequence_mapper import align_chain_to_axis, build_sequence_axis
from openfoldpanel.parsers.structure_preprocessor import normalize_structure_to_pdb
from openfoldpanel.parsers.structure_parser import (
    collect_reference_chains,
    get_chain_or_best_match,
    parse_structure,
)
from openfoldpanel.render.layout import build_render_config
from openfoldpanel.render.pdf_export import export_pdf
from openfoldpanel.render.report_renderer import reference_chain_pdf_name, render_html_report, render_reference_chain_report_svg
from openfoldpanel.utils.filesystem import ensure_directory, safe_rmtree
from openfoldpanel.utils.text import humanize_model_name, safe_chain_slug

QUERY_SEQUENCE_IDENTIFIER = "Query Sequence"


@dataclass(slots=True)
class JobIssues:
    warnings: list[str] = field(default_factory=list)
    partial_reasons: list[str] = field(default_factory=list)

    def add_warning(self, message: str) -> None:
        self.warnings.append(message)

    def add_partial(self, message: str) -> None:
        self.warnings.append(message)
        self.partial_reasons.append(message)

    def extend_warnings(self, messages: Iterable[str]) -> None:
        self.warnings.extend(messages)

    def extend_partials(self, messages: Iterable[str]) -> None:
        collected = list(messages)
        self.warnings.extend(collected)
        self.partial_reasons.extend(collected)


def run_pipeline(config: PipelineConfig, logger: logging.Logger) -> dict[str, int]:
    """Run the full batch pipeline and return a status summary."""

    logger.info(
        "Starting OpenFoldPanel run: input=%s outdir=%s chain=%s max_homologs_displayed=%s evalue=%s",
        config.input_path,
        config.outdir,
        config.chain,
        config.max_homologs_displayed,
        config.evalue,
    )
    validate_input_path(config.input_path)
    ensure_directory(config.outdir)

    temp_dir = Path(tempfile.mkdtemp(prefix="openfoldpanel_"))
    try:
        jobs = _discover_jobs(config, temp_dir, logger)
        logger.info("Discovered %s job(s) to process", len(jobs))
        results = [
            _run_job(job, config, logger, temp_dir, job_index=job_index, total_jobs=len(jobs))
            for job_index, job in enumerate(jobs, start=1)
        ]
        return _summarize_job_results(results)
    finally:
        if not config.keep_temp:
            logger.info("Cleaning up temporary workspace %s", temp_dir)
            safe_rmtree(temp_dir)
        else:
            logger.info("Keeping temporary workspace at %s", temp_dir)


def _discover_jobs(config: PipelineConfig, temp_dir: Path, logger: logging.Logger) -> list[JobDefinition]:
    if is_supported_archive(config.input_path):
        extracted_root = temp_dir / "archive"
        logger.info("Extracting archive input to %s", extracted_root)
        extract_archive(config.input_path, extracted_root)
        return discover_jobs_from_extracted_root(extracted_root, logger)

    logger.info("Discovering job from structure input %s", config.input_path)
    return discover_jobs_from_structure(config.input_path)


def _summarize_job_results(results: list[JobRunResult]) -> dict[str, int]:
    return {
        "total_jobs": len(results),
        "success": sum(result.status == "success" for result in results),
        "partial_success": sum(result.status == "partial_success" for result in results),
        "failed": sum(result.status == "failed" for result in results),
    }


def _run_job(
    job: JobDefinition,
    config: PipelineConfig,
    logger: logging.Logger,
    temp_dir: Path,
    *,
    job_index: int,
    total_jobs: int,
) -> JobRunResult:
    output_dir = ensure_directory(config.outdir / job.name)
    log_path = output_dir / "logs.txt"
    file_handler = attach_file_logger(logger, log_path)
    issues = JobIssues(warnings=list(job.ignored_files))
    try:
        logger.info(
            "[Job %s/%s] Starting %s (%s structure file(s))",
            job_index,
            total_jobs,
            job.name,
            len(job.structure_files),
        )
        if not job.structure_files:
            return _write_job_result(
                output_dir,
                JobRunResult(job_name=job.name, status="failed", output_dir=str(output_dir), error="No structure files found."),
                None,
            )

        parsed_structures = _parse_job_structures(
            job,
            logger,
            issues,
            temp_dir=temp_dir,
            job_index=job_index,
            total_jobs=total_jobs,
        )

        if not parsed_structures:
            return _write_job_result(
                output_dir,
                JobRunResult(
                    job_name=job.name,
                    status="failed",
                    output_dir=str(output_dir),
                    warnings=issues.warnings,
                    error="All models failed to parse.",
                ),
                None,
            )

        logger.info(
            "[Job %s/%s] Parsed %s/%s model(s) successfully",
            job_index,
            total_jobs,
            len(parsed_structures),
            len(job.structure_files),
        )
        reference_structure = parsed_structures[0]
        requested_reference_chains = collect_reference_chains(reference_structure, config.chain)
        preferred_default_reference_chain = _preferred_default_reference_chain(requested_reference_chains)
        render_config = build_render_config(config.columns, config.font_size, config.max_homologs_displayed)
        logger.info(
            "[Job %s/%s] Rendering reference chains: %s",
            job_index,
            total_jobs,
            ", ".join(requested_reference_chains),
        )

        chain_panels = _build_chain_panels(
            parsed_structures=parsed_structures,
            requested_reference_chains=requested_reference_chains,
            config=config,
            render_config=render_config,
            temp_dir=temp_dir,
            job_name=job.name,
            logger=logger,
            issues=issues,
            job_index=job_index,
            total_jobs=total_jobs,
        )

        if not chain_panels:
            return _write_job_result(
                output_dir,
                JobRunResult(
                    job_name=job.name,
                    status="failed",
                    output_dir=str(output_dir),
                    warnings=issues.warnings,
                    error="No compatible models could be mapped to any protein reference chain.",
                ),
                None,
            )

        default_reference_chain = _resolve_default_reference_chain(preferred_default_reference_chain, chain_panels)

        report_data = JobReportData(
            job_name=job.name,
            default_reference_chain=default_reference_chain,
            chain_panels=chain_panels,
            warnings=issues.warnings,
            status="success",
        )

        batch_result = build_batch_analysis(
            report_data,
            parsed_structures,
            tm_cluster_cutoff=config.tm_cluster_cutoff,
            disable_tm_clustering=config.disable_tm_clustering,
            logger=logger,
        )
        report_data.batch_analysis = batch_result.analysis
        issues.extend_warnings(batch_result.warnings)
        issues.extend_partials(batch_result.partial_reasons)

        artifacts = _export_reference_chain_pdfs(
            report_data,
            output_dir,
            logger,
            issues,
            job_index=job_index,
            total_jobs=total_jobs,
        )

        job_status = _status_from_partial_reasons(issues.partial_reasons)
        report_data.status = job_status

        artifacts.extend(
            _write_report_artifacts(
                report_data,
                config,
                output_dir,
                logger,
                job_index=job_index,
                total_jobs=total_jobs,
            )
        )
        artifacts.append("summary.txt")

        result = JobRunResult(
            job_name=job.name,
            status=job_status,
            output_dir=str(output_dir),
            warnings=issues.warnings,
            artifacts=artifacts,
        )
        logger.info("[Job %s/%s] Finished %s with status=%s", job_index, total_jobs, job.name, job_status)
        return _write_job_result(output_dir, result, report_data)
    except Exception as exc:
        logger.exception("Job %s failed", job.name)
        return _write_job_result(
            output_dir,
            JobRunResult(
                job_name=job.name,
                status="failed",
                output_dir=str(output_dir),
                warnings=issues.warnings,
                error=str(exc),
            ),
            None,
        )
    finally:
        detach_handler(logger, file_handler)


def _write_job_result(
    output_dir: Path,
    result: JobRunResult,
    report_data: JobReportData | None,
) -> JobRunResult:
    write_summary(result, report_data, output_dir / "summary.txt")
    return result


def _parse_job_structures(
    job: JobDefinition,
    logger: logging.Logger,
    issues: JobIssues,
    *,
    temp_dir: Path,
    job_index: int,
    total_jobs: int,
) -> list[ParsedStructure]:
    parsed_structures: list[ParsedStructure] = []
    for structure_index, structure_path in enumerate(job.structure_files, start=1):
        logger.info(
            "[Job %s/%s] Parsing model %s/%s: %s",
            job_index,
            total_jobs,
            structure_index,
            len(job.structure_files),
            structure_path.name,
        )
        try:
            normalized_dir = _normalized_structure_output_dir(temp_dir, job, structure_path)
            normalized_path = normalize_structure_to_pdb(structure_path, normalized_dir)
        except Exception as exc:
            message = f"Failed to normalize {structure_path.name}: {exc}"
            logger.warning(message)
            issues.add_partial(message)
            continue
        try:
            parsed_structures.append(
                parse_structure(
                    normalized_path,
                    logger,
                    original_source_path=_display_source_path(job, structure_path),
                )
            )
        except Exception as exc:
            message = f"Failed to parse {structure_path.name}: {exc}"
            logger.warning(message)
            issues.add_partial(message)
    return parsed_structures


def _display_source_path(job: JobDefinition, structure_path: Path) -> Path:
    try:
        return structure_path.relative_to(job.root_dir)
    except ValueError:
        return Path(structure_path.name)


def _normalized_structure_output_dir(temp_dir: Path, job: JobDefinition, structure_path: Path) -> Path:
    normalized_root = temp_dir / "normalized_structures" / job.name
    try:
        relative_parent = structure_path.parent.relative_to(job.root_dir)
    except ValueError:
        relative_parent = Path()
    return ensure_directory(normalized_root / relative_parent)


def _build_chain_panels(
    *,
    parsed_structures: list[ParsedStructure],
    requested_reference_chains: list[str],
    config: PipelineConfig,
    render_config: RenderConfig,
    temp_dir: Path,
    job_name: str,
    logger: logging.Logger,
    issues: JobIssues,
    job_index: int,
    total_jobs: int,
) -> list[JobPanelData]:
    chain_panels: list[JobPanelData] = []
    for chain_index, reference_chain_id in enumerate(requested_reference_chains, start=1):
        logger.info(
            "[Job %s/%s] Building chain %s (%s/%s)",
            job_index,
            total_jobs,
            reference_chain_id,
            chain_index,
            len(requested_reference_chains),
        )
        panel_data, chain_partial_reasons = _build_panel_data_for_reference_chain(
            parsed_structures=parsed_structures,
            reference_chain_id=reference_chain_id,
            config=config,
            render_config=render_config,
            workdir=temp_dir / job_name / safe_chain_slug(reference_chain_id),
            logger=logger,
            job_name=job_name,
        )
        if panel_data is None:
            message = f"Skipped {reference_chain_id}: no compatible models could be mapped to the reference axis."
            logger.warning(message)
            issues.add_partial(message)
            continue

        chain_panels.append(panel_data)
        issues.extend_warnings(_prefix_messages(f"Chain {reference_chain_id}", panel_data.warnings))
        issues.extend_partials(_prefix_messages(f"Chain {reference_chain_id}", chain_partial_reasons))
    return chain_panels


def _resolve_default_reference_chain(default_reference_chain: str, chain_panels: list[JobPanelData]) -> str:
    if default_reference_chain in {panel.reference_chain for panel in chain_panels}:
        return default_reference_chain
    return chain_panels[0].reference_chain


def _preferred_default_reference_chain(requested_reference_chains: list[str]) -> str:
    """Pick the preferred default chain before any per-chain build failures."""

    if "A" in requested_reference_chains:
        return "A"
    return requested_reference_chains[0]


def _prefix_messages(prefix: str, messages: Iterable[str]) -> list[str]:
    return [f"{prefix}: {message}" for message in messages]


def _status_from_partial_reasons(partial_reasons: list[str]) -> str:
    return "partial_success" if partial_reasons else "success"


def _export_reference_chain_pdfs(
    report_data: JobReportData,
    output_dir: Path,
    logger: logging.Logger,
    issues: JobIssues,
    *,
    job_index: int,
    total_jobs: int,
) -> list[str]:
    artifacts: list[str] = []
    for panel_data in report_data.chain_panels:
        pdf_path = output_dir / reference_chain_pdf_name(panel_data.reference_chain)
        logger.info(
            "[Job %s/%s] Rendering report SVG and exporting PDF for chain %s",
            job_index,
            total_jobs,
            panel_data.reference_chain,
        )
        report_svg = render_reference_chain_report_svg(panel_data)
        pdf_ok, pdf_warning = export_pdf(report_svg, pdf_path)
        if pdf_ok:
            artifacts.append(pdf_path.name)
        elif pdf_warning:
            issues.add_partial(pdf_warning)
    return artifacts


def _write_report_artifacts(
    report_data: JobReportData,
    config: PipelineConfig,
    output_dir: Path,
    logger: logging.Logger,
    *,
    job_index: int,
    total_jobs: int,
) -> list[str]:
    artifacts: list[str] = []

    html_path = output_dir / "report.html"
    logger.info("[Job %s/%s] Writing HTML report to %s", job_index, total_jobs, html_path)
    html_path.write_text(render_html_report(report_data, config), encoding="utf-8")
    artifacts.append(html_path.name)

    json_path = output_dir / "tracks.json"
    logger.info("[Job %s/%s] Writing tracks JSON to %s", job_index, total_jobs, json_path)
    write_tracks_json(report_data, json_path)
    artifacts.append(json_path.name)

    logger.info("[Job %s/%s] Writing CSV statistics exports to %s", job_index, total_jobs, output_dir / "csv")
    artifacts.extend(write_statistics_csvs(report_data, output_dir))

    return artifacts


def _build_panel_data_for_reference_chain(
    *,
    parsed_structures: list[ParsedStructure],
    reference_chain_id: str,
    config: PipelineConfig,
    render_config: RenderConfig,
    workdir: Path,
    logger: logging.Logger,
    job_name: str,
) -> tuple[JobPanelData | None, list[str]]:
    reference_structure = parsed_structures[0]
    reference_chain = reference_structure.chains[reference_chain_id]
    logger.info("Job %s / Chain %s: building sequence axis (%s residues)", job_name, reference_chain_id, len(reference_chain.residues))
    axis = build_sequence_axis(reference_chain)
    logger.info("Job %s / Chain %s: computing hydropathy track", job_name, reference_chain_id)
    hydropathy = compute_hydropathy(axis, config.hyd_window)
    antibody_numberings, antibody_warnings = annotate_antibody_chain(
        reference_chain.sequence,
        axis,
        chain_id=reference_chain_id,
    )

    warnings: list[str] = []
    warnings.extend(antibody_warnings)
    partial_reasons: list[str] = []
    model_tracks: list[ModelTracks] = []
    for structure_index, structure in enumerate(parsed_structures, start=1):
        display_name = structure.display_source_path.name
        logger.info(
            "Job %s / Chain %s: processing model %s/%s (%s)",
            job_name,
            reference_chain_id,
            structure_index,
            len(parsed_structures),
            display_name,
        )
        chain = get_chain_or_best_match(structure, reference_chain_id, reference_chain.sequence)
        if chain is None:
            warning = f"Skipping {display_name}: no compatible protein chain found."
            logger.warning(warning)
            warnings.append(warning)
            partial_reasons.append(warning)
            continue

        logger.info("Job %s / Chain %s: aligning %s to reference axis", job_name, reference_chain_id, display_name)
        alignment = align_chain_to_axis(axis, chain)
        warnings.extend(alignment.warnings)

        logger.info("Job %s / Chain %s: running DSSP for %s", job_name, reference_chain_id, display_name)
        dssp_features, dssp_warnings = _run_dssp_for_structure(structure, logger)
        warnings.extend(dssp_warnings)

        logger.info("Job %s / Chain %s: computing tracks for %s", job_name, reference_chain_id, display_name)
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
        disulfides = infer_disulfides(
            structure,
            alignment.residue_by_axis_index,
            current_chain_id=chain.chain_id,
        )
        model_tracks.append(
            ModelTracks(
                name=f"{structure.name}_{chain.chain_id}",
                source_path=str(structure.display_source_path),
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
        job_name=job_name,
        reference_chain_id=reference_chain_id,
    )
    warnings.extend(msa.warnings)
    panel_status = "partial_success" if partial_reasons else "success"
    panel_data = JobPanelData(
        job_name=job_name,
        reference_chain=reference_chain_id,
        sequence_axis=axis,
        models=model_tracks,
        msa=msa,
        hydropathy=hydropathy,
        render_config=render_config,
        antibody_numberings=antibody_numberings,
        warnings=warnings,
        status=panel_status,
    )
    return panel_data, partial_reasons


def _run_dssp_for_structure(
    structure: ParsedStructure,
    logger: logging.Logger,
):
    return run_dssp(
        structure.source_path,
        logger,
        display_name=structure.display_source_path.name,
    )


def _build_msa_data(
    *,
    axis: list[SequenceAxisPosition],
    config: PipelineConfig,
    workdir: Path,
    logger: logging.Logger,
    job_name: str,
    reference_chain_id: str,
) -> MSAData:
    ensure_directory(workdir)
    query_sequence = "".join(position.one_letter for position in axis)
    logger.info("Job %s / Chain %s: preparing MSA stage", job_name, reference_chain_id)
    if config.disable_msa:
        logger.info("Job %s / Chain %s: MSA disabled by user", job_name, reference_chain_id)
        return _build_query_only_msa(query_sequence, warnings=["MSA disabled by user."])

    if config.max_homologs_displayed == 0:
        logger.info("Job %s / Chain %s: max_homologs_displayed=0, skipping homolog search", job_name, reference_chain_id)
        return _build_query_only_msa(query_sequence)

    if config.msa_db is None:
        logger.info("Job %s / Chain %s: no MSA database configured, skipping homolog search", job_name, reference_chain_id)
        return _build_query_only_msa(query_sequence, warnings=["No MSA database provided; sequence alignment was skipped."])

    query_fasta = workdir / "query.fasta"
    query_fasta.write_text(f">query\n{query_sequence}\n", encoding="utf-8")
    # max_homologs_displayed caps how many hits we keep, while evalue filters
    # which candidate hits are significant enough to pass the search backend.
    search_limit = config.max_homologs_displayed
    logger.info(
        "Job %s / Chain %s: searching homologs in %s (display limit=%s candidate limit=%s evalue=%s)",
        job_name,
        reference_chain_id,
        config.msa_db,
        config.max_homologs_displayed,
        search_limit,
        config.evalue,
    )
    hits, hit_warnings = search_homologs(
        query_fasta,
        config.msa_db,
        max_homologs_displayed=search_limit,
        evalue=config.evalue,
        workdir=workdir,
        logger=logger,
    )
    if not hits:
        logger.info("Job %s / Chain %s: no homolog hits found", job_name, reference_chain_id)
        return _build_query_only_msa(query_sequence, warnings=hit_warnings)

    logger.info(
        "Job %s / Chain %s: found %s homolog hit(s), running alignment",
        job_name,
        reference_chain_id,
        len(hits),
    )
    alignment_rows, align_warnings = align_sequences(
        [(QUERY_SEQUENCE_IDENTIFIER, query_sequence), *hits],
        workdir / "alignment.fasta",
        logger,
    )
    if not alignment_rows:
        logger.info("Job %s / Chain %s: alignment failed or returned no rows", job_name, reference_chain_id)
        return _build_query_only_msa(query_sequence, warnings=[*hit_warnings, *align_warnings])

    projected_rows = _project_alignment_to_query_axis(alignment_rows)
    leading_display_overrides = _build_leading_display_overrides(alignment_rows, projected_rows)
    displayed_rows, displayed_overrides, filtered_count = _select_display_msa_rows_with_overrides(
        projected_rows,
        leading_display_overrides,
        max_homologs_displayed=config.max_homologs_displayed,
    )
    if filtered_count:
        logger.info(
            "Job %s / Chain %s: trimmed %s homolog row(s) to satisfy display limit",
            job_name,
            reference_chain_id,
            filtered_count,
        )
    if len(displayed_rows) <= 1:
        message = "No homolog rows remained after selection."
        logger.info("Job %s / Chain %s: %s", job_name, reference_chain_id, message)
        return _build_query_only_msa(
            query_sequence,
            warnings=[*hit_warnings, *align_warnings, message],
            rows=displayed_rows[:1],
            leading_display_overrides=displayed_overrides[:1],
        )
    logger.info(
        "Job %s / Chain %s: MSA ready with %s homolog row(s)",
        job_name,
        reference_chain_id,
        max(0, len(displayed_rows) - 1),
    )
    return MSAData(
        enabled=True,
        query=query_sequence,
        rows=displayed_rows,
        conservation=compute_conservation(displayed_rows),
        warnings=hit_warnings + align_warnings,
        leading_display_overrides=displayed_overrides,
    )


def _build_query_row(query_sequence: str) -> MSARow:
    return MSARow(identifier=QUERY_SEQUENCE_IDENTIFIER, sequence=query_sequence, is_query=True)


def _build_query_only_msa(
    query_sequence: str,
    *,
    warnings: list[str] | None = None,
    rows: list[MSARow] | None = None,
    leading_display_overrides: list[str | None] | None = None,
) -> MSAData:
    resolved_rows = rows if rows else [_build_query_row(query_sequence)]
    resolved_overrides = leading_display_overrides if leading_display_overrides else [None] * len(resolved_rows)
    return MSAData(
        enabled=False,
        query=query_sequence,
        rows=resolved_rows,
        warnings=list(warnings or []),
        leading_display_overrides=resolved_overrides,
    )


def _project_alignment_to_query_axis(rows: list[MSARow]) -> list[MSARow]:
    query = next((row for row in rows if row.is_query), rows[0])
    keep_indices = [index for index, residue in enumerate(query.sequence) if residue != "-"]
    projected: list[MSARow] = []
    for row in rows:
        sequence = "".join(row.sequence[index] for index in keep_indices if index < len(row.sequence))
        projected.append(MSARow(identifier=row.identifier, sequence=sequence, is_query=row.is_query))
    return projected


def _build_leading_display_overrides(aligned_rows: list[MSARow], projected_rows: list[MSARow]) -> list[str | None]:
    if not aligned_rows or not projected_rows:
        return []

    query = next((row for row in aligned_rows if row.is_query), aligned_rows[0])
    keep_indices = [index for index, residue in enumerate(query.sequence) if residue != "-"]
    overrides: list[str | None] = []
    for aligned_row, projected_row in zip(aligned_rows, projected_rows, strict=True):
        if aligned_row.is_query:
            overrides.append(None)
            continue
        overrides.append(_leading_gap_display_override(projected_row.sequence, aligned_row.sequence, keep_indices))
    return overrides


def _leading_gap_display_override(projected_sequence: str, aligned_sequence: str, keep_indices: list[int]) -> str | None:
    if not projected_sequence or projected_sequence[0] != "-" or not keep_indices:
        return None

    first_visible_index = next((index for index, residue in enumerate(projected_sequence) if residue != "-"), None)
    if first_visible_index is None or first_visible_index >= len(keep_indices):
        return None

    anchor_index = keep_indices[first_visible_index]
    start_index = min(anchor_index - 1, len(aligned_sequence) - 1)
    for index in range(start_index, -1, -1):
        residue = aligned_sequence[index]
        if residue != "-":
            return residue
    return None


def _select_display_msa_rows(rows: list[MSARow], *, max_homologs_displayed: int) -> tuple[list[MSARow], int]:
    displayed_rows, _, filtered_count = _select_display_msa_rows_with_overrides(
        rows,
        [None] * len(rows),
        max_homologs_displayed=max_homologs_displayed,
    )
    return displayed_rows, filtered_count


def _select_display_msa_rows_with_overrides(
    rows: list[MSARow],
    leading_display_overrides: list[str | None],
    *,
    max_homologs_displayed: int,
) -> tuple[list[MSARow], list[str | None], int]:
    if not rows:
        return [], [], 0

    query_index = next((index for index, row in enumerate(rows) if row.is_query), 0)
    selected: list[MSARow] = [rows[query_index]]
    selected_overrides: list[str | None] = [
        leading_display_overrides[query_index] if query_index < len(leading_display_overrides) else None
    ]
    kept_homologs = 0
    total_homologs = sum(1 for row in rows if not row.is_query)
    for row_index, row in enumerate(rows):
        if row.is_query:
            continue
        if kept_homologs >= max_homologs_displayed:
            continue
        selected.append(row)
        selected_overrides.append(
            leading_display_overrides[row_index] if row_index < len(leading_display_overrides) else None
        )
        kept_homologs += 1
    filtered_count = max(0, total_homologs - kept_homologs)
    return selected, selected_overrides, filtered_count
