"""Plain-text writers for summary and logs."""

from __future__ import annotations

from pathlib import Path

from openfoldpanel.models import JobPanelData, JobReportData, JobRunResult
from openfoldpanel.utils.text import humanize_chain_label, humanize_identifier, humanize_job_status


def write_summary(job_result: JobRunResult, panel_data: JobPanelData | JobReportData | None, path: Path) -> None:
    """Write a readable summary text file for one job."""

    lines = [
        f"Job Name: {humanize_identifier(job_result.job_name)}",
        f"Job Status: {humanize_job_status(job_result.status)}",
        f"Output Directory: {job_result.output_dir}",
    ]
    if panel_data is not None:
        if isinstance(panel_data, JobReportData):
            rendered_chains = ", ".join(humanize_chain_label(panel.reference_chain) for panel in panel_data.chain_panels)
            lines.extend(
                [
                    f"Default Reference Chain: {humanize_chain_label(panel_data.default_reference_chain)}",
                    f"Rendered Reference Chains: {rendered_chains}",
                    f"Rendered Chain Count: {len(panel_data.chain_panels)}",
                ]
            )
        else:
            lines.extend(
                [
                    f"Reference Chain: {humanize_chain_label(panel_data.reference_chain)}",
                    f"Sequence Length: {len(panel_data.sequence_axis)}",
                    f"Aligned Models: {len(panel_data.models)}",
                    f"Sequence Alignment: {'Full alignment' if panel_data.msa.enabled else 'Query only'}",
                ]
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
