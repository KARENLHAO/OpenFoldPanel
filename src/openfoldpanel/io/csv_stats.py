"""CSV statistics aggregation for job-level structural summaries."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from statistics import mean

from openfoldpanel.constants import SUPPORTED_ANTIBODY_NUMBERING, SUPPORTED_STRUCTURE_SUFFIXES
from openfoldpanel.io.csv_format import write_display_csv, write_raw_csv
from openfoldpanel.models import ContactConsensusScope, JobPanelData, JobReportData
from openfoldpanel.utils.filesystem import ensure_directory
from openfoldpanel.utils.reporting import contact_hits_for_entry, scope_sort_key

SECONDARY_CATEGORY_COLUMNS = (
    "alpha_helix",
    "three_ten_helix",
    "pi_helix",
    "strand",
    "alpha_turn",
    "beta_turn",
    "coil",
)
ACCESSIBILITY_MODE_ORDER = ("buried", "intermediate", "accessible", "highly_exposed")
CDR_REGION_ORDER = ("CDR1", "CDR2", "CDR3")


def write_statistics_csvs(report_data: JobReportData, output_dir: Path) -> list[str]:
    """Write one set of job-level CSV statistics tables."""

    csv_dir = ensure_directory(output_dir / "csv")
    artifacts: list[str] = []

    antibody_rows: list[dict[str, object]] = []
    contact_consensus_rows: list[dict[str, object]] = []
    tm_cluster_rows: list[dict[str, object]] = []
    tm_score_matrix_rows: list[dict[str, object]] = []
    tm_score_matrix_fieldnames: list[str] = []

    for panel_data in report_data.chain_panels:
        panel_residue_rows = _build_residue_summary_rows(panel_data)
        antibody_rows.extend(_build_antibody_summary_rows(panel_data, panel_residue_rows))

    if report_data.batch_analysis is not None:
        contact_consensus_rows = _build_contact_consensus_rows(report_data)
        tm_score_matrix_fieldnames, tm_score_matrix_rows = _build_tm_score_matrix_rows(report_data)
        tm_cluster_rows = _build_tm_cluster_rows(report_data)
    antibody_rows.sort(
        key=lambda row: (
            str(row["reference_chain"]),
            str(row["scheme"]),
            CDR_REGION_ORDER.index(str(row["region_name"])) if str(row["region_name"]) in CDR_REGION_ORDER else 99,
        )
    )
    contact_consensus_rows.sort(key=lambda row: _contact_consensus_cluster_sort_key(row["cluster_id"]))

    if contact_consensus_rows:
        contact_consensus_path = csv_dir / "contact-consensus.csv"
        write_display_csv(contact_consensus_path, _contact_consensus_fieldnames(), contact_consensus_rows)
        artifacts.append(contact_consensus_path.relative_to(output_dir).as_posix())

    if tm_score_matrix_fieldnames and tm_score_matrix_rows:
        tm_score_matrix_path = csv_dir / "tm-score-matrix.csv"
        write_raw_csv(tm_score_matrix_path, tm_score_matrix_fieldnames, tm_score_matrix_rows)
        artifacts.append(tm_score_matrix_path.relative_to(output_dir).as_posix())

    if tm_cluster_rows:
        tm_cluster_path = csv_dir / "tm-clusters.csv"
        write_display_csv(tm_cluster_path, _tm_cluster_fieldnames(), tm_cluster_rows)
        artifacts.append(tm_cluster_path.relative_to(output_dir).as_posix())

    if antibody_rows:
        antibody_path = csv_dir / "antibody-summary.csv"
        write_display_csv(antibody_path, _antibody_summary_fieldnames(), antibody_rows)
        artifacts.append(antibody_path.relative_to(output_dir).as_posix())

    return artifacts


def _build_residue_summary_rows(panel_data: JobPanelData) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    models_total = max(len(panel_data.models), 1)
    region_lookup = _region_lookup_by_scheme(panel_data)
    conservation_by_index = {entry.residue_index: entry for entry in panel_data.msa.conservation}
    chain_type = _residue_summary_chain_type(panel_data)

    for pos, position in enumerate(panel_data.sequence_axis):
        per_model_entries = [model.contacts[pos] for model in panel_data.models if pos < len(model.contacts)]
        per_model_hits = [contact_hits_for_entry(entry) for entry in per_model_entries]
        contact_entries = [entry for entry, hits in zip(per_model_entries, per_model_hits, strict=True) if hits]
        primary_distances = [entry.min_distance for entry in contact_entries if entry.min_distance is not None]
        accessibility_values = [
            model.accessibility[pos].relative
            for model in panel_data.models
            if pos < len(model.accessibility) and model.accessibility[pos].relative is not None
        ]
        accessibility_categories = [
            model.accessibility[pos].category
            for model in panel_data.models
            if pos < len(model.accessibility) and model.accessibility[pos].category is not None
        ]
        plddt_values = [
            value
            for model in panel_data.models
            if pos < len(model.plddt)
            for value in [model.plddt[pos]]
            if value is not None
        ]
        secondary_categories = [
            model.secondary_structure[pos].category
            for model in panel_data.models
            if pos < len(model.secondary_structure)
        ]
        conservation = conservation_by_index.get(pos)

        row: dict[str, object] = {
            "reference_chain": panel_data.reference_chain,
            "chain_type": chain_type,
            "pos": pos,
            "seq_id": position.seq_id,
            "insertion_code": position.insertion_code,
            "uid": position.label,
            "one_letter": position.one_letter,
            "contact_model_count": len(contact_entries),
            "contact_model_fraction": _round(len(contact_entries) / models_total),
            "strong_contact_model_count": sum(
                any(hit.strength_category == "strong" for hit in hits) for hits in per_model_hits
            ),
            "weak_contact_model_count": sum(
                any(hit.strength_category == "weak" for hit in hits) for hits in per_model_hits
            ),
            "protein_contact_model_count": sum(
                any(hit.partner_type == "protein_chain" for hit in hits) for hits in per_model_hits
            ),
            "ligand_contact_model_count": sum(
                any(hit.partner_type != "protein_chain" for hit in hits) for hits in per_model_hits
            ),
            "ion_contact_model_count": sum(any(hit.partner_type == "ion" for hit in hits) for hits in per_model_hits),
            "nucleic_acid_contact_model_count": sum(
                any(hit.partner_type == "nucleic_acid" for hit in hits) for hits in per_model_hits
            ),
            "sugar_contact_model_count": sum(any(hit.partner_type == "sugar" for hit in hits) for hits in per_model_hits),
            "porphyrin_like_contact_model_count": sum(
                any(hit.partner_type == "porphyrin_like" for hit in hits) for hits in per_model_hits
            ),
            "other_ligand_contact_model_count": sum(
                any(hit.partner_type == "other_ligand" for hit in hits) for hits in per_model_hits
            ),
            "multi_contact_model_count": sum(entry.is_multi_contact for entry in per_model_entries),
            "closest_contact_distance_min": _round(min(primary_distances), digits=3) if primary_distances else None,
            "closest_contact_distance_mean": _mean_or_none(primary_distances, digits=4),
            "plddt_mean": _mean_or_none(plddt_values, digits=4),
            "plddt_min": _round(min(plddt_values), digits=4) if plddt_values else None,
            "plddt_max": _round(max(plddt_values), digits=4) if plddt_values else None,
            "accessibility_relative_mean": _mean_or_none(accessibility_values, digits=4),
            "accessibility_category_mode": _mode(accessibility_categories, preferred_order=ACCESSIBILITY_MODE_ORDER),
            "hydropathy_value": panel_data.hydropathy[pos].value if pos < len(panel_data.hydropathy) else None,
            "hydropathy_category": panel_data.hydropathy[pos].category if pos < len(panel_data.hydropathy) else None,
            "conservation_identity_fraction": conservation.identity_fraction if conservation else None,
            "conservation_similarity_fraction": conservation.similarity_fraction if conservation else None,
            "conservation_style": conservation.style if conservation else None,
            "disulfide_model_count": sum(
                any(
                    bond.residue_index_a == pos or bond.residue_index_b == pos
                    for bond in model.disulfides
                )
                for model in panel_data.models
            ),
        }

        secondary_counter = Counter(category for category in secondary_categories if category in SECONDARY_CATEGORY_COLUMNS)
        for category in SECONDARY_CATEGORY_COLUMNS:
            row[f"secondary_{category}_count"] = secondary_counter.get(category, 0)

        for scheme in SUPPORTED_ANTIBODY_NUMBERING:
            row[f"{scheme}_region"] = region_lookup.get(scheme, {}).get(pos, "")

        rows.append(row)
    return rows


def _build_antibody_summary_rows(
    panel_data: JobPanelData,
    residue_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    if not panel_data.antibody_numberings:
        return []

    residue_rows_by_pos = {int(row["pos"]): row for row in residue_rows}
    rows: list[dict[str, object]] = []
    for scheme in SUPPORTED_ANTIBODY_NUMBERING:
        annotation = panel_data.antibody_numberings.get(scheme)
        if annotation is None:
            continue
        region_by_name = {region.name: region for region in annotation.regions}
        for region_name in CDR_REGION_ORDER:
            region = region_by_name.get(region_name)
            if region is None:
                continue
            region_positions = [index for index in range(region.start, min(region.end, len(panel_data.sequence_axis)))]
            if not region_positions:
                continue
            residue_subset = [residue_rows_by_pos[index] for index in region_positions]
            contact_subset = [row for row in residue_subset if int(row["contact_model_count"]) > 0]
            occupancies = [float(row["contact_model_fraction"]) for row in residue_subset]
            rows.append(
                {
                    "reference_chain": panel_data.reference_chain,
                    "scheme": scheme,
                    "chain_type": annotation.chain_type,
                    "region_name": region_name,
                    "region_start_pos": region.start,
                    "region_end_pos": region.end - 1,
                    "region_length": region.end - region.start,
                    "contact_site_count": len(contact_subset),
                    "strong_contact_site_count": sum(int(row["strong_contact_model_count"]) > 0 for row in residue_subset),
                    "mean_site_occupancy_fraction": _mean_or_none(occupancies, digits=4),
                    "max_site_occupancy_fraction": _round(max(occupancies), digits=4) if occupancies else None,
                    "conserved_contact_site_count": sum(
                        float(row["contact_model_fraction"]) >= 0.5
                        and row["conservation_identity_fraction"] is not None
                        and float(row["conservation_identity_fraction"]) >= 0.7
                        for row in residue_subset
                    ),
                }
            )
    return rows


def _build_contact_consensus_rows(report_data: JobReportData) -> list[dict[str, object]]:
    if report_data.batch_analysis is None:
        return []
    grouped_scopes: dict[str, list[ContactConsensusScope]] = {}
    for scope in report_data.batch_analysis.contact_consensus.scopes:
        grouped_scopes.setdefault(scope.scope, []).append(scope)

    rows: list[dict[str, object]] = []
    for scope_name in sorted(grouped_scopes, key=scope_sort_key):
        scope_rows = sorted(grouped_scopes[scope_name], key=lambda item: str(item.reference_chain))
        first_scope = scope_rows[0]
        rows.append(
            {
                "cluster_id": _contact_consensus_cluster_id(scope_name),
                "structure_count": first_scope.model_count,
                "cluster_center_structure": _csv_structure_name(first_scope.cluster_center_structure),
                "union_count": sum(scope.union_count for scope in scope_rows),
                "intersection_count": sum(scope.intersection_count for scope in scope_rows),
                "combine_residue": ",".join(scope.union_positions for scope in scope_rows if scope.union_positions),
                "consensus_residue": ",".join(
                    scope.intersection_positions for scope in scope_rows if scope.intersection_positions
                ),
            }
        )
    return rows


def _contact_consensus_cluster_id(scope_name: str) -> str:
    if scope_name == "global":
        return "all"
    if scope_name.startswith("cluster_"):
        suffix = scope_name.removeprefix("cluster_")
        if suffix:
            return suffix
    return scope_name


def _contact_consensus_cluster_sort_key(cluster_id: object) -> tuple[int, int | str]:
    text = str(cluster_id)
    if text in {"all", "global"}:
        return (0, 0)
    if text.isdigit():
        return (1, int(text))
    return (2, text)


def _build_tm_score_matrix_rows(report_data: JobReportData) -> tuple[list[str], list[dict[str, object]]]:
    if report_data.batch_analysis is None or not report_data.batch_analysis.tm_score.available:
        return [], []
    structure_names = [_csv_structure_name(name) for name in report_data.batch_analysis.tm_score.structure_names]
    first_column = "Structure"
    rows: list[dict[str, object]] = []
    for row_name, values in zip(structure_names, report_data.batch_analysis.tm_score.matrix, strict=True):
        row = {first_column: row_name}
        for column_name, value in zip(structure_names, values, strict=True):
            row[column_name] = value
        rows.append(row)
    return [first_column, *structure_names], rows


def _build_tm_cluster_rows(report_data: JobReportData) -> list[dict[str, object]]:
    if report_data.batch_analysis is None or not report_data.batch_analysis.tm_score.available:
        return []
    rows: list[dict[str, object]] = []
    for assignment in report_data.batch_analysis.tm_score.assignments:
        rows.append(
            {
                "structure_name": _csv_structure_name(assignment.structure_name),
                "cluster_id": assignment.cluster_id,
                "cluster_size": assignment.cluster_size,
                "cluster_center": _csv_structure_name(assignment.cluster_center),
                "is_representative": "1" if assignment.is_representative else "0",
            }
        )
    return rows


def _region_lookup_by_scheme(panel_data: JobPanelData) -> dict[str, dict[int, str]]:
    lookup: dict[str, dict[int, str]] = {}
    axis_length = len(panel_data.sequence_axis)
    for scheme, annotation in panel_data.antibody_numberings.items():
        scheme_lookup: dict[int, str] = {}
        for region in annotation.regions:
            for index in range(max(0, region.start), min(region.end, axis_length)):
                scheme_lookup[index] = region.name
        lookup[scheme] = scheme_lookup
    return lookup


def _mode(values: list[str], *, preferred_order: tuple[str, ...] = ()) -> str | None:
    if not values:
        return None
    counter = Counter(values)
    preferred_rank = {value: index for index, value in enumerate(preferred_order)}
    return min(
        counter,
        key=lambda value: (
            -counter[value],
            preferred_rank.get(value, len(preferred_order)),
            value,
        ),
    )


def _mean_or_none(values: list[float], *, digits: int = 4) -> float | None:
    if not values:
        return None
    return _round(mean(values), digits=digits)


def _round(value: float, *, digits: int = 4) -> float:
    return round(float(value), digits)


def _csv_structure_name(value: str) -> str:
    text = str(value)
    suffix = Path(text).suffix.lower()
    if suffix in SUPPORTED_STRUCTURE_SUFFIXES:
        return str(Path(text).with_suffix(""))
    return text


def _residue_summary_chain_type(panel_data: JobPanelData) -> str:
    for annotation in panel_data.antibody_numberings.values():
        normalized = str(annotation.chain_type).strip().lower()
        if normalized in {"h", "heavy"}:
            return "H"
        if normalized in {"l", "light"}:
            return "L"
    return "-"


def _antibody_summary_fieldnames() -> list[str]:
    return [
        "reference_chain",
        "scheme",
        "chain_type",
        "region_name",
        "region_start_pos",
        "region_end_pos",
        "region_length",
        "contact_site_count",
        "strong_contact_site_count",
        "mean_site_occupancy_fraction",
        "max_site_occupancy_fraction",
        "conserved_contact_site_count",
    ]


def _contact_consensus_fieldnames() -> list[str]:
    return [
        "cluster_id",
        "structure_count",
        "cluster_center_structure",
        "union_count",
        "intersection_count",
        "combine_residue",
        "consensus_residue",
    ]


def _tm_cluster_fieldnames() -> list[str]:
    return [
        "structure_name",
        "cluster_id",
        "cluster_size",
        "cluster_center",
        "is_representative",
    ]
