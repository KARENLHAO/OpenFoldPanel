"""Batch-level contact consensus and TM-score clustering analysis."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

from openfoldpanel.models import (
    BatchAnalysis,
    ContactConsensusAnalysis,
    ContactConsensusResidue,
    ContactConsensusScope,
    JobReportData,
    ParsedStructure,
    TMScoreAnalysis,
    TMScoreCluster,
    TMScoreClusterAssignment,
)
from openfoldpanel.utils.reporting import contact_entry_has_hit, scope_sort_key
from openfoldpanel.utils.subprocess import ExternalToolError, MissingExecutableError, run_command, which

AVERAGE_TM_SCORE_PATTERNS = (
    re.compile(r"Average TM-score=\s*([0-9.]+)\s+\(normalized by average L=", re.IGNORECASE),
    re.compile(r"TM-score=\s*([0-9.]+)\s+\(normalized by average L=", re.IGNORECASE),
    re.compile(r"Average TM-score=\s*([0-9.]+)", re.IGNORECASE),
    re.compile(r"TM-score=\s*([0-9.]+)", re.IGNORECASE),
)


@dataclass(slots=True)
class BatchAnalysisResult:
    analysis: BatchAnalysis
    warnings: list[str]
    partial_reasons: list[str]


def build_batch_analysis(
    report_data: JobReportData,
    parsed_structures: list[ParsedStructure],
    *,
    tm_cluster_cutoff: float,
    disable_tm_clustering: bool,
    logger: logging.Logger,
) -> BatchAnalysisResult:
    """Build job-level TM-score clustering and contact-consensus summaries."""

    tm_score, tm_warnings, tm_partial_reasons = _build_tm_score_analysis(
        parsed_structures,
        cutoff=tm_cluster_cutoff,
        disable_tm_clustering=disable_tm_clustering,
        logger=logger,
    )
    contact_consensus = _build_contact_consensus(report_data, tm_score)
    return BatchAnalysisResult(
        analysis=BatchAnalysis(tm_score=tm_score, contact_consensus=contact_consensus),
        warnings=tm_warnings,
        partial_reasons=tm_partial_reasons,
    )


def _build_tm_score_analysis(
    parsed_structures: list[ParsedStructure],
    *,
    cutoff: float,
    disable_tm_clustering: bool,
    logger: logging.Logger,
) -> tuple[TMScoreAnalysis, list[str], list[str]]:
    structure_names = [str(structure.display_source_path) for structure in parsed_structures]
    if disable_tm_clustering:
        return (
            TMScoreAnalysis(
                enabled=False,
                available=False,
                cutoff=cutoff,
                structure_names=structure_names,
            ),
            [],
            [],
        )

    if not parsed_structures:
        return (
            TMScoreAnalysis(enabled=True, available=False, cutoff=cutoff, structure_names=[]),
            [],
            [],
        )

    if len(parsed_structures) == 1:
        structure_name = structure_names[0]
        cluster = TMScoreCluster(cluster_id=1, size=1, center_structure=structure_name, members=[structure_name], mean_cluster_tm_score=1.0)
        assignment = TMScoreClusterAssignment(
            structure_name=structure_name,
            cluster_id=1,
            cluster_size=1,
            cluster_center=structure_name,
            is_representative=True,
            mean_intra_cluster_tm_score=1.0,
        )
        return (
            TMScoreAnalysis(
                enabled=True,
                available=True,
                cutoff=cutoff,
                structure_names=structure_names,
                matrix=[[1.0]],
                clusters=[cluster],
                assignments=[assignment],
            ),
            [],
            [],
        )

    executable = _find_usalign_executable()
    if executable is None:
        warning = "US-align executable was not found on PATH; TM-score clustering was skipped."
        return (
            TMScoreAnalysis(
                enabled=True,
                available=False,
                cutoff=cutoff,
                structure_names=structure_names,
                warnings=[warning],
            ),
            [warning],
            [warning],
        )

    try:
        matrix = _build_tm_score_matrix(parsed_structures, executable=executable, logger=logger)
    except (ExternalToolError, MissingExecutableError, ValueError) as exc:
        warning = f"US-align TM-score clustering was skipped: {exc}"
        return (
            TMScoreAnalysis(
                enabled=True,
                available=False,
                cutoff=cutoff,
                structure_names=structure_names,
                warnings=[warning],
            ),
            [warning],
            [warning],
        )

    cluster_groups = _cluster_by_average_linkage(matrix, cutoff)
    clusters, assignments = _build_cluster_outputs(structure_names, matrix, cluster_groups)
    return (
        TMScoreAnalysis(
            enabled=True,
            available=True,
            cutoff=cutoff,
            structure_names=structure_names,
            matrix=matrix,
            clusters=clusters,
            assignments=assignments,
        ),
        [],
        [],
    )


def _find_usalign_executable() -> str | None:
    return which("USalign") or which("US-align")


def _build_tm_score_matrix(
    parsed_structures: list[ParsedStructure],
    *,
    executable: str,
    logger: logging.Logger,
) -> list[list[float]]:
    count = len(parsed_structures)
    matrix = [[0.0] * count for _ in range(count)]
    for index in range(count):
        matrix[index][index] = 1.0

    for index_a in range(count):
        for index_b in range(index_a + 1, count):
            structure_a = parsed_structures[index_a]
            structure_b = parsed_structures[index_b]
            tm_score = _run_usalign_pair(structure_a, structure_b, executable=executable, logger=logger)
            matrix[index_a][index_b] = tm_score
            matrix[index_b][index_a] = tm_score
    return matrix


def _run_usalign_pair(
    structure_a: ParsedStructure,
    structure_b: ParsedStructure,
    *,
    executable: str,
    logger: logging.Logger,
) -> float:
    is_multimer = _protein_chain_count(structure_a) > 1 or _protein_chain_count(structure_b) > 1
    command = [
        executable,
        str(structure_a.source_path),
        str(structure_b.source_path),
        "-mol",
        "prot",
        "-a",
        "T",
        "-outfmt",
        "-1",
    ]
    if is_multimer:
        command.extend(["-mm", "1", "-ter", "1"])

    logger.info(
        "Batch analysis: running US-align for %s vs %s%s",
        structure_a.display_source_path.name,
        structure_b.display_source_path.name,
        " (multimer mode)" if is_multimer else "",
    )
    result = run_command(command)
    return _parse_average_tm_score(result.stdout)


def _protein_chain_count(structure: ParsedStructure) -> int:
    return sum(chain.entity_type == "protein" for chain in structure.chains.values())


def _parse_average_tm_score(output_text: str) -> float:
    for pattern in AVERAGE_TM_SCORE_PATTERNS:
        match = pattern.search(output_text)
        if match:
            return round(float(match.group(1)), 5)
    raise ValueError("Could not parse average TM-score from US-align output.")


def _cluster_by_average_linkage(matrix: list[list[float]], cutoff: float) -> list[list[int]]:
    max_distance = 1.0 - cutoff
    clusters = [{index} for index in range(len(matrix))]
    if len(clusters) <= 1:
        return [sorted(cluster) for cluster in clusters]

    while True:
        best_pair: tuple[int, int] | None = None
        best_distance: float | None = None
        for left in range(len(clusters)):
            for right in range(left + 1, len(clusters)):
                distance = _average_linkage_distance(matrix, clusters[left], clusters[right])
                if best_distance is None or distance < best_distance - 1e-12 or (
                    abs(distance - best_distance) <= 1e-12
                    and _cluster_sort_key(clusters[left], clusters[right]) < _cluster_sort_key(clusters[best_pair[0]], clusters[best_pair[1]])
                ):
                    best_distance = distance
                    best_pair = (left, right)
        if best_pair is None or best_distance is None or best_distance > max_distance:
            break

        left, right = best_pair
        merged = set(clusters[left]) | set(clusters[right])
        clusters = [cluster for index, cluster in enumerate(clusters) if index not in {left, right}]
        clusters.append(merged)
        clusters.sort(key=lambda cluster: min(cluster))

    return [sorted(cluster) for cluster in sorted(clusters, key=lambda cluster: min(cluster))]


def _average_linkage_distance(matrix: list[list[float]], cluster_a: set[int], cluster_b: set[int]) -> float:
    similarities = [matrix[index_a][index_b] for index_a in cluster_a for index_b in cluster_b]
    return 1.0 - (sum(similarities) / len(similarities))


def _cluster_sort_key(cluster_a: set[int], cluster_b: set[int]) -> tuple[int, int]:
    return (min(cluster_a), min(cluster_b))


def _build_cluster_outputs(
    structure_names: list[str],
    matrix: list[list[float]],
    cluster_groups: list[list[int]],
) -> tuple[list[TMScoreCluster], list[TMScoreClusterAssignment]]:
    clusters: list[TMScoreCluster] = []
    assignments: list[TMScoreClusterAssignment] = []

    for cluster_id, member_indices in enumerate(cluster_groups, start=1):
        member_scores = {index: _mean_cluster_similarity(matrix, index, member_indices) for index in member_indices}
        center_index = min(member_indices, key=lambda index: (-member_scores[index], index))
        mean_cluster_score = round(sum(member_scores.values()) / len(member_scores), 5) if member_scores else None
        clusters.append(
            TMScoreCluster(
                cluster_id=cluster_id,
                size=len(member_indices),
                center_structure=structure_names[center_index],
                members=[structure_names[index] for index in member_indices],
                mean_cluster_tm_score=mean_cluster_score,
            )
        )
        for member_index in member_indices:
            assignments.append(
                TMScoreClusterAssignment(
                    structure_name=structure_names[member_index],
                    cluster_id=cluster_id,
                    cluster_size=len(member_indices),
                    cluster_center=structure_names[center_index],
                    is_representative=member_index == center_index,
                    mean_intra_cluster_tm_score=member_scores[member_index],
                )
            )

    assignments.sort(key=lambda row: structure_names.index(row.structure_name))
    return clusters, assignments


def _mean_cluster_similarity(matrix: list[list[float]], structure_index: int, member_indices: list[int]) -> float:
    values = [matrix[structure_index][member_index] for member_index in member_indices]
    return round(sum(values) / len(values), 5)


def _build_contact_consensus(report_data: JobReportData, tm_score: TMScoreAnalysis) -> ContactConsensusAnalysis:
    scope_definitions = _build_scope_definitions(tm_score)
    scopes: list[ContactConsensusScope] = []
    residues: list[ContactConsensusResidue] = []

    for panel_data in report_data.chain_panels:
        for scope_name, structure_names, center_structure in scope_definitions:
            if structure_names is None:
                selected_models = list(panel_data.models)
            else:
                selected_models = [model for model in panel_data.models if model.source_path in structure_names]
            if not selected_models:
                continue

            model_count = len(selected_models)
            union_positions: list[int] = []
            intersection_positions: list[int] = []
            selected_structure_names = [model.source_path for model in selected_models]

            for axis_index, position in enumerate(panel_data.sequence_axis):
                occurrence_count = sum(
                    axis_index < len(model.contacts) and contact_entry_has_hit(model.contacts[axis_index])
                    for model in selected_models
                )
                if occurrence_count == 0:
                    continue
                axis_position = axis_index + 1
                in_intersection = occurrence_count == model_count
                union_positions.append(axis_position)
                if in_intersection:
                    intersection_positions.append(axis_position)
                residues.append(
                    ContactConsensusResidue(
                        scope=scope_name,
                        reference_chain=panel_data.reference_chain,
                        cluster_center_structure=center_structure,
                        model_count=model_count,
                        axis_position=axis_position,
                        axis_label=f"{panel_data.reference_chain}{axis_position}",
                        seq_id=position.seq_id,
                        insertion_code=position.insertion_code,
                        uid=position.label,
                        one_letter=position.one_letter,
                        occurrence_count=occurrence_count,
                        occurrence_fraction=round(occurrence_count / model_count, 4),
                        in_intersection=in_intersection,
                    )
                )

            scopes.append(
                ContactConsensusScope(
                    scope=scope_name,
                    reference_chain=panel_data.reference_chain,
                    model_count=model_count,
                    cluster_center_structure=center_structure,
                    structure_names=selected_structure_names,
                    union_count=len(union_positions),
                    intersection_count=len(intersection_positions),
                    union_positions=_compress_axis_positions(panel_data.reference_chain, union_positions),
                    intersection_positions=_compress_axis_positions(panel_data.reference_chain, intersection_positions),
                    union_sequence=_sequence_for_axis_positions(panel_data, union_positions),
                    intersection_sequence=_sequence_for_axis_positions(panel_data, intersection_positions),
                )
            )

    scopes.sort(key=lambda row: (row.reference_chain, scope_sort_key(row.scope)))
    residues.sort(key=lambda row: (row.reference_chain, scope_sort_key(row.scope), row.axis_position))
    return ContactConsensusAnalysis(scopes=scopes, residues=residues)


def _build_scope_definitions(tm_score: TMScoreAnalysis) -> list[tuple[str, set[str] | None, str]]:
    scopes: list[tuple[str, set[str] | None, str]] = [("global", None, "")]
    if not tm_score.available:
        return scopes
    for cluster in tm_score.clusters:
        scopes.append((f"cluster_{cluster.cluster_id}", set(cluster.members), cluster.center_structure))
    return scopes


def _compress_axis_positions(reference_chain: str, positions: list[int]) -> str:
    if not positions:
        return ""
    ordered = sorted(set(positions))
    ranges: list[str] = []
    start = ordered[0]
    end = ordered[0]
    for position in ordered[1:]:
        if position == end + 1:
            end = position
            continue
        ranges.append(_format_axis_range(reference_chain, start, end))
        start = end = position
    ranges.append(_format_axis_range(reference_chain, start, end))
    return ",".join(ranges)


def _format_axis_range(reference_chain: str, start: int, end: int) -> str:
    if start == end:
        return f"{reference_chain}{start}"
    return f"{reference_chain}{start}-{end}"


def _sequence_for_axis_positions(panel_data, positions: list[int]) -> str:
    if not positions:
        return ""
    return "".join(panel_data.sequence_axis[position - 1].one_letter for position in sorted(positions))

