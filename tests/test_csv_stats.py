from __future__ import annotations

import csv

from openfoldpanel.features.batch_analysis import _build_contact_consensus
from openfoldpanel.features.hydropathy import compute_hydropathy
from openfoldpanel.io.csv_stats import write_statistics_csvs
from openfoldpanel.models import (
    AccessibilityEntry,
    AntibodyAnnotation,
    BatchAnalysis,
    ContactEntry,
    ContactHit,
    ContactConsensusAnalysis,
    ContactConsensusScope,
    ConservationEntry,
    JobPanelData,
    JobReportData,
    MSAData,
    MSARow,
    ModelTracks,
    RegionAnnotation,
    SecondaryStructureEntry,
    SequenceAxisPosition,
    TMScoreAnalysis,
    TMScoreCluster,
    TMScoreClusterAssignment,
)
from openfoldpanel.render.layout import build_render_config


def _build_panel() -> JobPanelData:
    axis = [
        SequenceAxisPosition(index, "A", index + 1, "", resname, residue, str(index + 1))
        for index, (resname, residue) in enumerate(
            [
                ("ALA", "A"),
                ("GLY", "G"),
                ("SER", "S"),
            ]
        )
    ]
    hydropathy = compute_hydropathy(axis, window=3)
    model_a = ModelTracks(
        name="model_a_A",
        source_path="model_a.pdb",
        chain="A",
        secondary_structure=[
            SecondaryStructureEntry(0, "H", "alpha_helix"),
            SecondaryStructureEntry(1, "E", "strand"),
            SecondaryStructureEntry(2, "C", "coil"),
        ],
        plddt=[95.0, 82.0, 60.0],
        accessibility=[
            AccessibilityEntry(0, 120.0, 0.72, "accessible"),
            AccessibilityEntry(1, 18.0, 0.08, "buried"),
            AccessibilityEntry(2, 44.0, 0.31, "intermediate"),
        ],
        contacts=[
            ContactEntry(
                0,
                "protein_chain",
                "B",
                "TYR",
                "1",
                2.8,
                "B",
                "strong",
                all_contacts=[
                    ContactHit("protein_chain", "B", "TYR", "1", 2.8, "B", "strong"),
                ],
            ),
            ContactEntry(1, None, None, None, None, None, None, None),
            ContactEntry(
                2,
                "ion",
                "Z",
                "ZN",
                "1",
                3.5,
                "+",
                "weak",
                all_contacts=[
                    ContactHit("ion", "Z", "ZN", "1", 3.5, "+", "weak"),
                ],
            ),
        ],
        display_name="Model A / Chain A",
    )
    model_b = ModelTracks(
        name="model_b_A",
        source_path="model_b.pdb",
        chain="A",
        secondary_structure=[
            SecondaryStructureEntry(0, "H", "alpha_helix"),
            SecondaryStructureEntry(1, "C", "coil"),
            SecondaryStructureEntry(2, "C", "coil"),
        ],
        plddt=[88.0, 70.0, None],
        accessibility=[
            AccessibilityEntry(0, 110.0, 0.66, "accessible"),
            AccessibilityEntry(1, 52.0, 0.28, "intermediate"),
            AccessibilityEntry(2, None, None, None),
        ],
        contacts=[
            ContactEntry(
                0,
                "protein_chain",
                "C",
                "TYR",
                "5",
                3.0,
                "C",
                "strong",
                is_multi_contact=True,
                all_contacts=[
                    ContactHit("protein_chain", "C", "TYR", "5", 3.0, "C", "strong"),
                    ContactHit("ion", "Z", "ZN", "1", 3.4, "+", "weak"),
                ],
            ),
            ContactEntry(1, None, None, None, None, None, None, None),
            ContactEntry(
                2,
                "nucleic_acid",
                "N",
                "DA",
                "7",
                3.3,
                "*",
                "weak",
                all_contacts=[
                    ContactHit("nucleic_acid", "N", "DA", "7", 3.3, "*", "weak"),
                ],
            ),
        ],
        display_name="Model B / Chain A",
    )

    return JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=[model_a, model_b],
        msa=MSAData(
            enabled=True,
            query="AGS",
            rows=[
                MSARow(identifier="Query Sequence", sequence="AGS", is_query=True),
                MSARow(identifier="hit_1", sequence="AGS"),
            ],
            conservation=[
                ConservationEntry(0, 0.8, 0.8, "similar"),
                ConservationEntry(1, 0.2, 0.2, "default"),
                ConservationEntry(2, 0.9, 0.9, "similar"),
            ],
        ),
        hydropathy=hydropathy,
        render_config=build_render_config(columns=3, font_size=12),
        antibody_numberings={
            "kabat": AntibodyAnnotation(
                scheme="kabat",
                chain_type="heavy",
                regions=[
                    RegionAnnotation("CDR1", 0, 1, "CDR1 - Kabat"),
                    RegionAnnotation("CDR2", 1, 2, "CDR2 - Kabat"),
                    RegionAnnotation("CDR3", 2, 3, "CDR3 - Kabat"),
                ],
            ),
            "imgt": AntibodyAnnotation(
                scheme="imgt",
                chain_type="heavy",
                regions=[
                    RegionAnnotation("CDR1", 0, 1, "CDR1 - IMGT"),
                    RegionAnnotation("CDR2", 1, 2, "CDR2 - IMGT"),
                    RegionAnnotation("CDR3", 2, 3, "CDR3 - IMGT"),
                ],
            ),
            "chothia": AntibodyAnnotation(
                scheme="chothia",
                chain_type="heavy",
                regions=[
                    RegionAnnotation("CDR1", 0, 1, "CDR1 - Chothia"),
                    RegionAnnotation("CDR2", 1, 2, "CDR2 - Chothia"),
                    RegionAnnotation("CDR3", 2, 3, "CDR3 - Chothia"),
                ],
            ),
        },
    )


def _read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            rows.append({_normalize_header(key): value for key, value in row.items()})
        return rows


def _normalize_header(value: str | None) -> str:
    if value is None:
        return ""
    normalized = value.lower().replace(" ", "_")
    normalized = normalized.replace("seq_id", "seq_id")
    normalized = normalized.replace("imgt", "imgt")
    normalized = normalized.replace("plddt", "plddt")
    return normalized


def test_write_statistics_csvs_exports_expected_tables(tmp_path):
    second_panel = _build_panel()
    second_panel.reference_chain = "B"
    second_panel.antibody_numberings = {}
    report = JobReportData(job_name="demo", default_reference_chain="A", chain_panels=[_build_panel(), second_panel])

    artifacts = write_statistics_csvs(report, tmp_path)

    assert artifacts == ["csv/antibody-summary.csv"]
    assert not (tmp_path / "csv" / "residue-summary.csv").exists()
    assert not (tmp_path / "csv" / "contacts.csv").exists()
    assert not (tmp_path / "csv" / "model-summary.csv").exists()
    assert not (tmp_path / "csv" / "tm-cluster-relations.csv").exists()

    antibody_rows = _read_csv(tmp_path / "csv" / "antibody-summary.csv")
    assert len(antibody_rows) == 9
    kabat_cdr1 = next(
        row for row in antibody_rows if row["reference_chain"] == "A" and row["scheme"] == "kabat" and row["region_name"] == "CDR1"
    )
    assert kabat_cdr1["contact_site_count"] == "1"
    assert kabat_cdr1["strong_contact_site_count"] == "1"
    assert kabat_cdr1["mean_site_occupancy_fraction"] == "1.0"
    assert kabat_cdr1["conserved_contact_site_count"] == "1"


def test_write_statistics_csvs_exports_batch_analysis_tables_when_available(tmp_path):
    panel = _build_panel()
    panel.models[1].contacts[2] = ContactEntry(2, None, None, None, None, None, None, None)
    report = JobReportData(job_name="demo", default_reference_chain="A", chain_panels=[panel])
    tm_score = TMScoreAnalysis(
        enabled=True,
        available=True,
        cutoff=0.7,
        structure_names=["model_a.pdb", "model_b.pdb"],
        matrix=[[1.0, 0.82], [0.82, 1.0]],
        clusters=[TMScoreCluster(cluster_id=1, size=2, center_structure="model_a.pdb", members=["model_a.pdb", "model_b.pdb"], mean_cluster_tm_score=0.91)],
        assignments=[
            TMScoreClusterAssignment("model_a.pdb", 1, 2, "model_a.pdb", True, 0.91),
            TMScoreClusterAssignment("model_b.pdb", 1, 2, "model_a.pdb", False, 0.91),
        ],
    )
    report.batch_analysis = BatchAnalysis(
        tm_score=tm_score,
        contact_consensus=_build_contact_consensus(report, tm_score),
    )

    artifacts = write_statistics_csvs(report, tmp_path)

    assert "csv/contact-consensus.csv" in artifacts
    assert "csv/tm-score-matrix.csv" in artifacts
    assert "csv/tm-clusters.csv" in artifacts
    assert "csv/tm-cluster-relations.csv" not in artifacts

    contact_consensus_rows = _read_csv(tmp_path / "csv" / "contact-consensus.csv")
    assert len(contact_consensus_rows) == 2
    contact_consensus_header_line = (tmp_path / "csv" / "contact-consensus.csv").read_text(encoding="utf-8").splitlines()[0]
    assert contact_consensus_header_line == (
        "Cluster Id,Structure Count,Cluster Center Structure,Combine Count,Consensus Count,Combine Residue,Consensus Residue"
    )
    global_row = next(row for row in contact_consensus_rows if row["cluster_id"] == "all")
    assert global_row["combine_residue"] == "A1,A3"
    assert global_row["consensus_residue"] == "A1"
    cluster_row = next(row for row in contact_consensus_rows if row["cluster_id"] == "1")
    assert cluster_row["structure_count"] == "2"
    assert cluster_row["cluster_center_structure"] == "model_a"
    assert "scope" not in cluster_row
    assert "reference_chain" not in cluster_row
    assert "union_sequence" not in cluster_row
    assert "intersection_sequence" not in cluster_row

    tm_score_matrix_lines = (tmp_path / "csv" / "tm-score-matrix.csv").read_text(encoding="utf-8").splitlines()
    assert tm_score_matrix_lines[0] == "Structure,model_a,model_b"
    assert tm_score_matrix_lines[1] == "model_a,1.0,0.82"

    tm_cluster_rows = _read_csv(tmp_path / "csv" / "tm-clusters.csv")
    assert len(tm_cluster_rows) == 2
    assert tm_cluster_rows[0]["structure"] == "model_a"
    assert tm_cluster_rows[0]["cluster_center"] == "model_a"
    assert tm_cluster_rows[0]["is_representative"] == "1"
    assert "mean_intra_cluster_tm_score" not in tm_cluster_rows[0]
    assert not (tmp_path / "csv" / "tm-cluster-relations.csv").exists()


def test_write_statistics_csvs_merges_contact_consensus_rows_by_cluster_id(tmp_path):
    report = JobReportData(
        job_name="demo",
        default_reference_chain="A",
        chain_panels=[],
        batch_analysis=BatchAnalysis(
            tm_score=TMScoreAnalysis(enabled=False, available=False, cutoff=0.7),
            contact_consensus=ContactConsensusAnalysis(
                scopes=[
                    ContactConsensusScope(
                        scope="global",
                        reference_chain="A",
                        model_count=3,
                        union_count=2,
                        intersection_count=1,
                        union_positions="A1,A3",
                        intersection_positions="A1",
                        union_sequence="AS",
                        intersection_sequence="A",
                    ),
                    ContactConsensusScope(
                        scope="global",
                        reference_chain="B",
                        model_count=3,
                        union_count=1,
                        intersection_count=1,
                        union_positions="B2",
                        intersection_positions="B2",
                        union_sequence="G",
                        intersection_sequence="G",
                    ),
                    ContactConsensusScope(
                        scope="cluster_2",
                        reference_chain="A",
                        model_count=2,
                        cluster_center_structure="model_b.pdb",
                        union_count=1,
                        intersection_count=1,
                        union_positions="A3",
                        intersection_positions="A3",
                        union_sequence="S",
                        intersection_sequence="S",
                    ),
                    ContactConsensusScope(
                        scope="cluster_2",
                        reference_chain="B",
                        model_count=2,
                        cluster_center_structure="model_b.pdb",
                        union_count=2,
                        intersection_count=1,
                        union_positions="B2,B4",
                        intersection_positions="B2",
                        union_sequence="GT",
                        intersection_sequence="G",
                    ),
                ]
            ),
        ),
    )

    write_statistics_csvs(report, tmp_path)

    contact_consensus_rows = _read_csv(tmp_path / "csv" / "contact-consensus.csv")

    assert [row["cluster_id"] for row in contact_consensus_rows] == ["all", "2"]
    assert contact_consensus_rows[0]["combine_count"] == "3"
    assert contact_consensus_rows[0]["consensus_count"] == "2"
    assert contact_consensus_rows[0]["structure_count"] == "3"
    assert contact_consensus_rows[0]["combine_residue"] == "A1,A3,B2"
    assert contact_consensus_rows[0]["consensus_residue"] == "A1,B2"
    assert "union_sequence" not in contact_consensus_rows[0]
    assert "intersection_sequence" not in contact_consensus_rows[0]

    cluster_row = contact_consensus_rows[1]
    assert cluster_row["cluster_center_structure"] == "model_b"
    assert cluster_row["combine_residue"] == "A3,B2,B4"
    assert cluster_row["consensus_residue"] == "A3,B2"
