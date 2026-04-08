"""Plain-text writers for summary and logs."""

from __future__ import annotations

from pathlib import Path

from openfoldpanel.models import JobReportData, JobRunResult
from openfoldpanel.utils.text import humanize_chain_label, humanize_identifier, humanize_job_status


def write_summary(job_result: JobRunResult, report_data: JobReportData | None, path: Path) -> None:
    """Write a readable summary text file for one job."""

    lines = [
        f"Job Name: {humanize_identifier(job_result.job_name)}",
        f"Job Status: {humanize_job_status(job_result.status)}",
        f"Output Directory: {job_result.output_dir}",
    ]
    if report_data is not None:
        rendered_chains = ", ".join(humanize_chain_label(panel.reference_chain) for panel in report_data.chain_panels)
        lines.extend(
            [
                f"Default Reference Chain: {humanize_chain_label(report_data.default_reference_chain)}",
                f"Rendered Reference Chains: {rendered_chains}",
                f"Rendered Chain Count: {len(report_data.chain_panels)}",
            ]
        )
        if report_data.batch_analysis is not None:
            tm_score = report_data.batch_analysis.tm_score
            if not tm_score.enabled:
                lines.append("TM-score Clustering: Disabled")
            elif not tm_score.available:
                lines.append("TM-score Clustering: Unavailable")
            else:
                lines.append(f"TM-score Cluster Count: {len(tm_score.clusters)}")
                for cluster in tm_score.clusters:
                    lines.append(
                        f"TM-score Cluster {cluster.cluster_id}: center={cluster.center_structure}; size={cluster.size}"
                    )
    if job_result.artifacts:
        lines.append("Artifacts:")
        lines.extend(f"  - {artifact}" for artifact in job_result.artifacts)
    if job_result.warnings:
        lines.append("Warnings:")
        lines.extend(f"  - {warning}" for warning in job_result.warnings)
    if job_result.error:
        lines.append(f"Error: {job_result.error}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
