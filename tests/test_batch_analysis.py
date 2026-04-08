from __future__ import annotations

import logging
from pathlib import Path

import openfoldpanel.features.batch_analysis as batch_analysis_module
from openfoldpanel.features.batch_analysis import (
    _build_contact_consensus,
    _build_tm_score_analysis,
    _compress_axis_positions,
    _run_usalign_pair,
)
from openfoldpanel.models import (
    AccessibilityEntry,
    ChainRecord,
    ContactEntry,
    ContactHit,
    JobPanelData,
    JobReportData,
    MSAData,
    ModelTracks,
    ParsedStructure,
    RenderConfig,
    SecondaryStructureEntry,
    SequenceAxisPosition,
    TMScoreAnalysis,
    TMScoreCluster,
)


def _contact_track(length: int, hit_positions: set[int]) -> list[ContactEntry]:
    track: list[ContactEntry] = []
    for index in range(length):
        if index + 1 in hit_positions:
            track.append(
                ContactEntry(
                    residue_index=index,
                    partner_type="protein_chain",
                    partner_chain="B",
                    partner_resname="TYR",
                    partner_resid=str(index + 1),
                    min_distance=3.0,
                    symbol="B",
                    strength_category="strong",
                    all_contacts=[ContactHit("protein_chain", "B", "TYR", str(index + 1), 3.0, "B", "strong")],
                )
            )
        else:
            track.append(ContactEntry(index, None, None, None, None, None, None, None))
    return track


def _panel_for_contact_consensus() -> JobPanelData:
    sequence = "MNGACSPQRTYVLKF"
    axis = [
        SequenceAxisPosition(index, "A", index + 101, "", "ALA", residue, str(index + 101))
        for index, residue in enumerate(sequence)
    ]
    config = RenderConfig(
        columns=80,
        max_homologs_displayed=5,
        font_size=12,
        cell_width=8.0,
        row_height=16.0,
        label_width=80.0,
        margin=16.0,
        colors={},
        font_family="Times",
        heading_font_family="Times",
    )
    models = [
        ModelTracks(
            name="m1_A",
            source_path="m1.pdb",
            chain="A",
            secondary_structure=[SecondaryStructureEntry(index, "C", "coil") for index in range(len(axis))],
            plddt=[90.0] * len(axis),
            accessibility=[AccessibilityEntry(index, None, None, None) for index in range(len(axis))],
            contacts=_contact_track(len(axis), {4, 5, 6, 8, 9, 10, 11, 12, 15}),
        ),
        ModelTracks(
            name="m2_A",
            source_path="m2.pdb",
            chain="A",
            secondary_structure=[SecondaryStructureEntry(index, "C", "coil") for index in range(len(axis))],
            plddt=[90.0] * len(axis),
            accessibility=[AccessibilityEntry(index, None, None, None) for index in range(len(axis))],
            contacts=_contact_track(len(axis), {4, 5, 6, 8, 9, 10, 11, 12, 15}),
        ),
        ModelTracks(
            name="m3_A",
            source_path="m3.pdb",
            chain="A",
            secondary_structure=[SecondaryStructureEntry(index, "C", "coil") for index in range(len(axis))],
            plddt=[90.0] * len(axis),
            accessibility=[AccessibilityEntry(index, None, None, None) for index in range(len(axis))],
            contacts=_contact_track(len(axis), {4}),
        ),
    ]
    return JobPanelData(
        job_name="demo",
        reference_chain="A",
        sequence_axis=axis,
        models=models,
        msa=MSAData(enabled=False, query=sequence),
        hydropathy=[],
        render_config=config,
    )


def _parsed_structure(name: str, protein_chain_count: int = 1) -> ParsedStructure:
    chains = {
        chr(ord("A") + index): ChainRecord(chain_id=chr(ord("A") + index), residues=[], entity_type="protein")
        for index in range(protein_chain_count)
    }
    return ParsedStructure(name=Path(name).stem, source_path=Path(f"/tmp/{name}"), chains=chains, format="pdb", original_source_path=Path(name))


def test_compress_axis_positions_uses_chain_prefixed_one_based_ranges():
    assert _compress_axis_positions("A", [4, 5, 6, 8, 9, 10, 11, 12, 15]) == "A4-6,A8-12,A15"


def test_build_contact_consensus_outputs_global_and_cluster_scopes():
    panel = _panel_for_contact_consensus()
    report = JobReportData(job_name="demo", default_reference_chain="A", chain_panels=[panel])
    tm_score = TMScoreAnalysis(
        enabled=True,
        available=True,
        cutoff=0.7,
        structure_names=["m1.pdb", "m2.pdb", "m3.pdb"],
        matrix=[[1.0, 0.82, 0.55], [0.82, 1.0, 0.56], [0.55, 0.56, 1.0]],
        clusters=[
            TMScoreCluster(cluster_id=1, size=2, center_structure="m1.pdb", members=["m1.pdb", "m2.pdb"], mean_cluster_tm_score=0.91),
            TMScoreCluster(cluster_id=2, size=1, center_structure="m3.pdb", members=["m3.pdb"], mean_cluster_tm_score=1.0),
        ],
    )

    consensus = _build_contact_consensus(report, tm_score)

    global_row = next(row for row in consensus.scopes if row.scope == "global" and row.reference_chain == "A")
    assert global_row.model_count == 3
    assert global_row.union_positions == "A4-6,A8-12,A15"
    assert global_row.intersection_positions == "A4"
    cluster_row = next(row for row in consensus.scopes if row.scope == "cluster_1" and row.reference_chain == "A")
    assert cluster_row.cluster_center_structure == "m1.pdb"
    assert cluster_row.union_positions == "A4-6,A8-12,A15"
    assert cluster_row.intersection_positions == "A4-6,A8-12,A15"

    residue_a15 = next(row for row in consensus.residues if row.scope == "global" and row.axis_label == "A15")
    assert residue_a15.occurrence_count == 2
    assert residue_a15.occurrence_fraction == 0.6667
    assert residue_a15.in_intersection is False
    residue_a4 = next(row for row in consensus.residues if row.scope == "global" and row.axis_label == "A4")
    assert residue_a4.occurrence_count == 3
    assert residue_a4.in_intersection is True


def test_build_tm_score_analysis_clusters_from_symmetric_pair_scores(monkeypatch):
    parsed_structures = [_parsed_structure("a.pdb"), _parsed_structure("b.pdb"), _parsed_structure("c.pdb")]
    scores = {
        ("a.pdb", "b.pdb"): 0.82,
        ("a.pdb", "c.pdb"): 0.58,
        ("b.pdb", "c.pdb"): 0.57,
    }

    monkeypatch.setattr(batch_analysis_module, "_find_usalign_executable", lambda: "USalign")
    monkeypatch.setattr(
        batch_analysis_module,
        "_run_usalign_pair",
        lambda structure_a, structure_b, **kwargs: scores[tuple(sorted((str(structure_a.display_source_path), str(structure_b.display_source_path))))],
    )

    tm_score, warnings, partial_reasons = _build_tm_score_analysis(
        parsed_structures,
        cutoff=0.7,
        disable_tm_clustering=False,
        logger=logging.getLogger("test"),
    )

    assert warnings == []
    assert partial_reasons == []
    assert tm_score.available is True
    assert tm_score.matrix == [
        [1.0, 0.82, 0.58],
        [0.82, 1.0, 0.57],
        [0.58, 0.57, 1.0],
    ]
    assert [cluster.center_structure for cluster in tm_score.clusters] == ["a.pdb", "c.pdb"]
    assert [assignment.cluster_id for assignment in tm_score.assignments] == [1, 1, 2]
    assert tm_score.assignments[0].cluster_center == "a.pdb"


def test_run_usalign_pair_uses_multimer_flags_when_needed(monkeypatch):
    commands: list[list[str]] = []

    def fake_run_command(command):
        commands.append(list(command))
        return type("Result", (), {"stdout": "Average TM-score= 0.88862 (normalized by average L=227.00 and d0=5.59)\n"})()

    monkeypatch.setattr(batch_analysis_module, "run_command", fake_run_command)

    monomer_a = _parsed_structure("mono_a.pdb", protein_chain_count=1)
    monomer_b = _parsed_structure("mono_b.pdb", protein_chain_count=1)
    multimer = _parsed_structure("multi.pdb", protein_chain_count=2)

    assert _run_usalign_pair(monomer_a, monomer_b, executable="USalign", logger=logging.getLogger("test")) == 0.88862
    assert _run_usalign_pair(multimer, monomer_b, executable="USalign", logger=logging.getLogger("test")) == 0.88862

    assert commands[0][:8] == [
        "USalign",
        "/tmp/mono_a.pdb",
        "/tmp/mono_b.pdb",
        "-mol",
        "prot",
        "-a",
        "T",
        "-outfmt",
    ]
    assert "-mm" not in commands[0]
    assert commands[1][-4:] == ["-mm", "1", "-ter", "1"]
