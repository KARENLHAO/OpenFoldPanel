from __future__ import annotations

import json
import tarfile
import zipfile
from importlib import resources
from types import SimpleNamespace

import pytest

import openfoldpanel.pipeline as pipeline_module
from openfoldpanel.cli import build_parser, main
from openfoldpanel.constants import ALLOWED_EVALUES, DEFAULT_EVALUE
from openfoldpanel.models import AntibodyAnnotation, BatchAnalysis, ContactConsensusAnalysis, RegionAnnotation, TMScoreAnalysis
from openfoldpanel.features.dssp_runner import DSSPResidueFeature
from openfoldpanel.utils.text import summarize_msa_database_path
from tests.conftest import write_test_mmcif, write_test_pdb


@pytest.fixture(autouse=True)
def fake_pdf_export_success(monkeypatch):
    def fake_export_pdf(svg_markup, output_path):
        output_path.write_bytes(b"%PDF-1.4\n% OpenFoldPanel test fixture\n")
        return True, None

    monkeypatch.setattr(pipeline_module, "export_pdf", fake_export_pdf)


def test_cli_default_outputs_multi_chain_reports(tmp_path):
    input_path = write_test_pdb(tmp_path / "single_model.pdb")
    outdir = tmp_path / "out"

    exit_code = main(["--input", str(input_path), "--outdir", str(outdir)])
    assert exit_code == 0

    job_dir = outdir / "single_model"
    assert (job_dir / "report.html").exists()
    assert (job_dir / "reference-chain-A.pdf").exists()
    assert (job_dir / "reference-chain-B.pdf").exists()
    assert (job_dir / "tracks.json").exists()
    assert (job_dir / "summary.txt").exists()
    assert (job_dir / "logs.txt").exists()
    assert not (job_dir / "csv" / "residue-summary.csv").exists()
    assert not (job_dir / "csv" / "contacts.csv").exists()
    assert not (job_dir / "csv" / "model-summary.csv").exists()
    assert (job_dir / "csv" / "contact-consensus.csv").exists()
    assert (job_dir / "csv" / "tm-score-matrix.csv").exists()
    assert (job_dir / "csv" / "tm-clusters.csv").exists()
    assert not (job_dir / "csv" / "tm-cluster-relations.csv").exists()
    assert not (job_dir / "csv" / "antibody-summary.csv").exists()
    assert not (job_dir / "panel.svg").exists()
    assert not (job_dir / "panel.png").exists()

    payload = json.loads((job_dir / "tracks.json").read_text())
    assert payload["default_reference_chain"] == "A"
    assert payload["status"] == "success"
    assert {panel["reference_chain"] for panel in payload["chain_panels"]} == {"A", "B"}
    assert {panel["job_name"] for panel in payload["chain_panels"]} == {"single_model"}
    assert "batch_analysis" in payload
    assert payload["batch_analysis"]["tm_score"]["matrix"] == [[1.0]]
    assert payload["batch_analysis"]["tm_score"]["clusters"][0]["center_structure"] == "single_model.pdb"
    assert payload["batch_analysis"]["contact_consensus"]["scopes"]
    assert payload["chain_panels"][0]["models"]
    assert any(model["display_name"] == "Single Model / Chain A" for model in payload["chain_panels"][0]["models"])
    assert "disulfides" in payload["chain_panels"][0]["models"][0]

    html_text = (job_dir / "report.html").read_text()
    assert 'lang="en"' in html_text
    assert "@font-face" in html_text
    assert "OpenFoldPanel Times New Roman" in html_text
    assert "Chain A" in html_text
    assert "Chain B" in html_text
    assert "Query Sequence" in html_text
    assert "Accessibility" in html_text
    assert "Hydropathy" in html_text
    assert "Dark blue indicates very low surface accessibility and mostly buried residues." in html_text
    assert "Grey-blue indicates partial exposure between buried and exposed states." in html_text
    assert "Light blue indicates residues that are more solvent-accessible." in html_text
    assert "Gold indicates residues with strong surface exposure." in html_text
    assert "Orange indicates locally hydrophobic regions, often buried or lipid-facing." in html_text
    assert "Gray indicates locally intermediate physicochemical character." in html_text
    assert "Blue indicates locally hydrophilic regions with stronger aqueous compatibility." in html_text
    assert "Intramolecular Disulfide" in html_text
    assert "Intermolecular Disulfide" in html_text
    assert "Confidence" in html_text
    assert "OpenFoldPanel / ARCHIVE" in html_text
    assert "FoldScript-style panel" not in html_text
    assert "Reference Chain" in html_text
    assert "Legend" in html_text
    assert "single_model_A" in html_text
    assert "Hydropathy Window" in html_text
    assert "E-value Threshold" in html_text
    assert "Weak Contact Cutoff" in html_text
    assert "Strong Contact Cutoff" in html_text
    assert "Homolog Display Limit" in html_text
    assert "Database" in html_text
    assert "3.7 A" in html_text
    assert "3.2 A" in html_text
    assert "Based on the shortest non-hydrogen atom distance: below 3.2 A." in html_text
    assert "Based on the shortest non-hydrogen atom distance: between 3.2 A and 3.7 A, inclusive." in html_text
    assert "5" in html_text
    assert "Not set" in html_text
    assert "--hyd-window" in html_text
    assert "--evalue" in html_text
    assert "--contact-cutoff" in html_text
    assert "--strong-contact-cutoff" in html_text
    assert "--max-homologs-displayed" in html_text
    assert "--msa-db" in html_text
    assert 'data-chain-select' in html_text
    assert 'data-antibody-scheme-select' in html_text
    assert 'data-report-page' in html_text
    assert 'data-active-chain-panel' in html_text
    assert 'data-figure-sheet' in html_text
    assert 'data-summary-grid' in html_text
    assert 'ofp-toolbar-summary' in html_text
    assert 'id="ofp-report-payload"' in html_text
    assert 'data-chain-templates' in html_text
    assert 'template id="ofp-chain-A-kabat"' in html_text
    assert 'template id="ofp-chain-B-kabat"' in html_text
    assert 'data-chain-figure="A"' in html_text
    assert 'data-chain-figure="B"' in html_text
    assert 'data-antibody-scheme="kabat"' in html_text
    assert 'data-panel-width="' in html_text
    assert "--active-panel-width:" in html_text
    assert 'ofp-figure-wrap figure-wrap' in html_text
    assert 'ofp-figure-sheet figure-sheet' in html_text
    assert 'data-legend-deck' in html_text
    assert 'data-legend-card="structure"' in html_text
    assert 'data-legend-card="tracks"' in html_text
    assert 'data-legend-card="contacts"' in html_text
    assert 'data-legend-kind="strand"' in html_text
    assert 'data-legend-kind="alpha-helix"' in html_text
    assert 'data-legend-kind="three-ten-helix"' in html_text
    assert 'data-legend-kind="pi-helix"' in html_text
    assert 'data-legend-kind="helix"' not in html_text
    assert 'data-legend-kind="alpha-turn"' in html_text
    assert 'data-legend-kind="beta-turn"' in html_text
    assert 'data-legend-kind="turn"' not in html_text
    assert 'data-antibody-legend hidden' in html_text
    assert 'data-legend-kind="antibody-numbering"' in html_text
    assert "Alpha Helix" in html_text
    assert "3₁₀ Helix" in html_text
    assert "Pi Helix" in html_text
    assert "Alpha Turn" in html_text
    assert "Beta Turn" in html_text
    assert "Antibody Numbering" in html_text
    assert "Shows CDR / framework region labels using the currently selected antibody numbering scheme." in html_text
    assert 'data-legend-kind="accessibility-buried"' in html_text
    assert 'data-legend-kind="accessibility-intermediate"' in html_text
    assert 'data-legend-kind="accessibility-accessible"' in html_text
    assert 'data-legend-kind="accessibility-highly-exposed"' in html_text
    assert 'data-legend-kind="hydropathy-hydrophobic"' in html_text
    assert 'data-legend-kind="hydropathy-intermediate"' in html_text
    assert 'data-legend-kind="hydropathy-hydrophilic"' in html_text
    assert 'data-legend-kind="confidence-very-high"' in html_text
    assert 'data-legend-kind="contact-disulfide-intra"' in html_text
    assert 'data-legend-kind="contact-disulfide-inter"' in html_text
    assert 'data-legend-kind="contact-multi"' in html_text
    assert 'data-current-chain-label' not in html_text
    assert 'supporting-rail' not in html_text
    assert html_text.index("OpenFoldPanel / ARCHIVE") < html_text.index('data-summary-grid')
    assert html_text.index('data-figure-sheet') < html_text.index('data-legend-deck')
    assert html_text.index('data-legend-deck') < html_text.index('data-chain-warnings')
    assert "fetch(" not in html_text
    assert "Export Current Chain PDF" not in html_text
    assert "Download Tracks JSON" not in html_text
    assert "View Summary" not in html_text
    assert "View Logs" not in html_text
    assert "Terminology" not in html_text
    assert "Research Collaboration Report" not in html_text

    summary_text = (job_dir / "summary.txt").read_text()
    assert "csv/residue-summary.csv" not in summary_text
    assert "csv/contacts.csv" not in summary_text
    assert "csv/model-summary.csv" not in summary_text
    assert "csv/contact-consensus.csv" in summary_text
    assert "csv/tm-score-matrix.csv" in summary_text
    assert "csv/tm-clusters.csv" in summary_text
    assert "csv/tm-cluster-relations.csv" not in summary_text
    assert "TM-score Cluster Count: 1" in summary_text


def test_cli_tracks_json_uses_helix_subtype_categories(tmp_path, monkeypatch):
    input_path = write_test_pdb(tmp_path / "single_model.pdb")
    outdir = tmp_path / "out"

    def fake_run_dssp(structure_path, logger, *, display_name=None):
        assert display_name == "single_model.pdb"
        return {
            ("A", 1, ""): DSSPResidueFeature("A", 1, "", "A", "G", 42.0),
            ("A", 2, ""): DSSPResidueFeature("A", 2, "", "G", "H", 42.0),
            ("A", 3, ""): DSSPResidueFeature("A", 3, "", "S", "I", 42.0),
        }, []

    monkeypatch.setattr(pipeline_module, "run_dssp", fake_run_dssp)

    exit_code = main(["--input", str(input_path), "--outdir", str(outdir), "--chain", "A"])
    assert exit_code == 0

    payload = json.loads((outdir / "single_model" / "tracks.json").read_text())
    categories = [entry["category"] for entry in payload["chain_panels"][0]["models"][0]["secondary_structure"][:3]]
    assert categories == ["three_ten_helix", "alpha_helix", "pi_helix"]
    assert "helix" not in categories


def test_cli_rejects_deprecated_auto_and_supports_specific_reference_chain_selection(tmp_path):
    input_path = write_test_pdb(tmp_path / "single_model.pdb")

    auto_outdir = tmp_path / "deprecated_auto_out"
    with pytest.raises(SystemExit) as exc_info:
        main(["--input", str(input_path), "--outdir", str(auto_outdir), "--chain", "AUTO"])
    assert exc_info.value.code == 2

    chain_b_outdir = tmp_path / "chain_b_out"
    assert main(["--input", str(input_path), "--outdir", str(chain_b_outdir), "--chain", "B"]) == 0
    chain_b_payload = json.loads((chain_b_outdir / "single_model" / "tracks.json").read_text())
    assert chain_b_payload["default_reference_chain"] == "B"
    assert len(chain_b_payload["chain_panels"]) == 1
    assert chain_b_payload["chain_panels"][0]["reference_chain"] == "B"
    assert (chain_b_outdir / "single_model" / "reference-chain-B.pdf").exists()


def test_cli_parser_accepts_max_homologs_displayed_in_valid_range():
    parser = build_parser()
    args = parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--max-homologs-displayed", "5"])
    assert args.max_homologs_displayed == 5

    args = parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--max-homologs-displayed", "0"])
    assert args.max_homologs_displayed == 0

    args = parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--max-homologs-displayed", "25"])
    assert args.max_homologs_displayed == 25


def test_cli_parser_uses_default_evalue_and_accepts_allowed_values():
    parser = build_parser()

    default_args = parser.parse_args(["--input", "demo.pdb", "--outdir", "out"])
    assert default_args.evalue == DEFAULT_EVALUE

    for evalue in ("1e-4", "1e-12"):
        args = parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--evalue", evalue])
        assert args.evalue == evalue


def test_cli_parser_accepts_tm_cluster_cutoff():
    parser = build_parser()
    args = parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--tm-cluster-cutoff", "0.8"])
    assert args.tm_cluster_cutoff == 0.8


def test_cli_parser_rejects_invalid_tm_cluster_cutoff(capsys):
    parser = build_parser()

    for cutoff in ("0", "1.2", "-0.1", "abc"):
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--tm-cluster-cutoff", cutoff])
        assert exc_info.value.code == 2
        assert "--tm-cluster-cutoff" in capsys.readouterr().err


def test_cli_parser_accepts_fasta_msa_database_path():
    parser = build_parser()
    args = parser.parse_args(
        [
            "--input",
            "demo.pdb",
            "--outdir",
            "out",
            "--msa-db",
            "blastdb/swissprot_fasta/uniprot_sprot.fasta",
        ]
    )
    assert str(args.msa_db) == "blastdb/swissprot_fasta/uniprot_sprot.fasta"


def test_cli_report_renders_uppercase_msa_database_name_from_path_tail(tmp_path):
    input_path = write_test_pdb(tmp_path / "single_model.pdb")
    outdir = tmp_path / "out"

    exit_code = main(
        [
            "--input",
            str(input_path),
            "--outdir",
            str(outdir),
            "--msa-db",
            "./blastdb/swissprot/swissprot",
        ]
    )
    assert exit_code == 0

    html_text = (outdir / "single_model" / "report.html").read_text()
    assert "Homolog Display Limit" in html_text
    assert "5" in html_text
    assert "Database" in html_text
    assert "SWISSPROT" in html_text
    assert html_text.index("Strong Contact Cutoff") < html_text.index("Homolog Display Limit")
    assert html_text.index("Homolog Display Limit") < html_text.index("Database")


def test_cli_tracks_json_includes_antibody_annotation_when_numbering_is_enabled(tmp_path, monkeypatch):
    input_path = write_test_pdb(tmp_path / "single_model.pdb")
    outdir = tmp_path / "out"

    def fake_annotate_antibody_chain(sequence, sequence_axis, *, chain_id):
        assert chain_id == "A"
        return (
            {
                "kabat": AntibodyAnnotation(
                    scheme="kabat",
                    chain_type="heavy",
                    regions=[
                        RegionAnnotation(name="CDR1", start=0, end=1, display_label="CDR1 - Kabat"),
                        RegionAnnotation(name="CDR2", start=1, end=2, display_label="CDR2 - Kabat"),
                        RegionAnnotation(name="CDR3", start=2, end=3, display_label="CDR3 - Kabat"),
                    ],
                ),
                "imgt": AntibodyAnnotation(
                    scheme="imgt",
                    chain_type="heavy",
                    regions=[
                        RegionAnnotation(name="CDR1", start=0, end=1, display_label="CDR1 - IMGT"),
                        RegionAnnotation(name="CDR2", start=1, end=2, display_label="CDR2 - IMGT"),
                        RegionAnnotation(name="CDR3", start=2, end=3, display_label="CDR3 - IMGT"),
                    ],
                ),
            },
            [],
        )

    monkeypatch.setattr(pipeline_module, "annotate_antibody_chain", fake_annotate_antibody_chain)

    exit_code = main(
        [
            "--input",
            str(input_path),
            "--outdir",
            str(outdir),
            "--chain",
            "A",
        ]
    )
    assert exit_code == 0

    payload = json.loads((outdir / "single_model" / "tracks.json").read_text())
    panel = payload["chain_panels"][0]
    assert panel["default_antibody_numbering_scheme"] == "kabat"
    assert set(panel["antibody_numberings"]) == {"kabat", "imgt"}
    annotation = panel["antibody_numberings"]["kabat"]
    assert annotation["scheme"] == "kabat"
    assert annotation["chain_type"] == "heavy"
    assert [region["name"] for region in annotation["regions"]] == ["CDR1", "CDR2", "CDR3"]
    html_text = (outdir / "single_model" / "report.html").read_text()
    assert 'data-antibody-scheme="kabat"' in html_text
    assert 'data-antibody-scheme="imgt"' in html_text
    assert 'data-antibody-legend' in html_text
    assert 'data-antibody-legend hidden' not in html_text
    assert 'data-legend-kind="antibody-numbering"' in html_text
    assert "Antibody Numbering" in html_text
    assert "Shows CDR / framework region labels using the currently selected antibody numbering scheme." in html_text


def test_cli_antibody_numbering_warning_does_not_fail_job(tmp_path, monkeypatch):
    input_path = write_test_pdb(tmp_path / "single_model.pdb")
    outdir = tmp_path / "out"

    monkeypatch.setattr(
        pipeline_module,
        "annotate_antibody_chain",
        lambda sequence, sequence_axis, *, chain_id: ({}, ["abnumber unavailable for test"]),
    )

    exit_code = main(
        [
            "--input",
            str(input_path),
            "--outdir",
            str(outdir),
            "--chain",
            "A",
        ]
    )
    assert exit_code == 0

    payload = json.loads((outdir / "single_model" / "tracks.json").read_text())
    panel = payload["chain_panels"][0]
    assert panel["antibody_numberings"] == {}
    assert "abnumber unavailable for test" in panel["warnings"]


def test_cli_report_renders_contact_legend_thresholds_from_config(tmp_path):
    input_path = write_test_pdb(tmp_path / "single_model.pdb")
    outdir = tmp_path / "out"

    exit_code = main(
        [
            "--input",
            str(input_path),
            "--outdir",
            str(outdir),
            "--contact-cutoff",
            "4.1",
            "--strong-contact-cutoff",
            "3.5",
        ]
    )
    assert exit_code == 0

    html_text = (outdir / "single_model" / "report.html").read_text()
    assert "Based on the shortest non-hydrogen atom distance: below 3.5 A." in html_text
    assert "Based on the shortest non-hydrogen atom distance: between 3.5 A and 4.1 A, inclusive." in html_text


def test_cli_parser_rejects_out_of_range_and_legacy_msa_flags():
    parser = build_parser()

    with pytest.raises(SystemExit) as too_high:
        parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--max-homologs-displayed", "26"])
    assert too_high.value.code == 2

    with pytest.raises(SystemExit) as negative:
        parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--max-homologs-displayed", "-1"])
    assert negative.value.code == 2

    with pytest.raises(SystemExit) as legacy_hits:
        parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--max-hits", "3"])
    assert legacy_hits.value.code == 2

    with pytest.raises(SystemExit) as legacy_rows:
        parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--msa-display-rows", "3"])
    assert legacy_rows.value.code == 2


def test_cli_parser_rejects_invalid_evalue_values(capsys):
    parser = build_parser()

    for evalue in ("1e-3", "0.000001", "1E-6"):
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--evalue", evalue])
        assert exc_info.value.code == 2
        error_output = capsys.readouterr().err
        assert "--evalue must be one of:" in error_output
        assert ", ".join(ALLOWED_EVALUES) in error_output


def test_cli_archive_batch_outputs_job_directories(tmp_path):
    archive_path = tmp_path / "batch.tgz"
    job_a = tmp_path / "job_a"
    job_b = tmp_path / "job_b"
    job_a.mkdir()
    job_b.mkdir()
    write_test_pdb(job_a / "ranked_1.pdb")
    write_test_pdb(job_b / "ranked_2.pdb")
    with tarfile.open(archive_path, "w:gz") as tf:
        tf.add(job_a, arcname="job_a")
        tf.add(job_b, arcname="job_b")

    outdir = tmp_path / "archive_out"
    exit_code = main(["--input", str(archive_path), "--outdir", str(outdir)])
    assert exit_code == 0
    assert (outdir / "job_a" / "report.html").exists()
    assert (outdir / "job_a" / "reference-chain-A.pdf").exists()
    assert (outdir / "job_b" / "tracks.json").exists()


def test_cli_marks_partial_success_when_pdf_export_fails(tmp_path, monkeypatch):
    input_path = write_test_pdb(tmp_path / "single_model.pdb")
    outdir = tmp_path / "out"

    monkeypatch.setattr(
        pipeline_module,
        "export_pdf",
        lambda svg_markup, output_path: (False, "CairoSVG is not installed; PDF export was skipped."),
    )

    exit_code = main(["--input", str(input_path), "--outdir", str(outdir)])
    assert exit_code == 0

    job_dir = outdir / "single_model"
    payload = json.loads((job_dir / "tracks.json").read_text())
    assert payload["status"] == "partial_success"
    assert (job_dir / "report.html").exists()
    assert not (job_dir / "reference-chain-A.pdf").exists()
    summary_text = (job_dir / "summary.txt").read_text()
    assert "Job Status: Partial Success" in summary_text
    assert "CairoSVG is not installed; PDF export was skipped." in summary_text


def test_cli_normalizes_cif_input_to_pdb_for_dssp_but_keeps_original_source_name(tmp_path, monkeypatch):
    input_path = write_test_mmcif(tmp_path / "single_model.cif")
    outdir = tmp_path / "out"
    seen_dssp_inputs = []

    def fake_run_dssp(structure_path, logger, *, display_name=None):
        seen_dssp_inputs.append(structure_path)
        assert display_name == "single_model.cif"
        return {}, []

    monkeypatch.setattr(pipeline_module, "run_dssp", fake_run_dssp)

    exit_code = main(["--input", str(input_path), "--outdir", str(outdir)])
    assert exit_code == 0

    assert seen_dssp_inputs
    assert all(path.suffix == ".pdb" for path in seen_dssp_inputs)
    payload = json.loads((outdir / "single_model" / "tracks.json").read_text())
    assert payload["chain_panels"][0]["models"][0]["source_path"].endswith("single_model.cif")
    logs_text = (outdir / "single_model" / "logs.txt").read_text()
    assert "running DSSP for single_model.cif" in logs_text
    assert "running DSSP for single_model.pdb" not in logs_text


def test_cli_marks_partial_success_when_one_archive_model_cannot_be_normalized_to_pdb(tmp_path):
    archive_path = tmp_path / "batch.zip"
    with zipfile.ZipFile(archive_path, "w") as zf:
        zf.writestr("job_a/ranked_1.pdb", write_test_pdb(tmp_path / "ranked_1.pdb").read_text())
        zf.writestr("job_a/ranked_bad.cif", write_test_mmcif(tmp_path / "ranked_bad.cif", chain_id="AA").read_text())

    outdir = tmp_path / "out"
    exit_code = main(["--input", str(archive_path), "--outdir", str(outdir)])
    assert exit_code == 0

    job_dir = outdir / "job_a"
    payload = json.loads((job_dir / "tracks.json").read_text())
    assert payload["status"] == "partial_success"
    assert (job_dir / "report.html").exists()
    summary_text = (job_dir / "summary.txt").read_text()
    assert "Failed to normalize ranked_bad.cif" in summary_text


def test_cli_archive_mixed_pdb_and_cif_inputs_are_normalized_before_dssp(tmp_path, monkeypatch):
    archive_path = tmp_path / "mixed.zip"
    with zipfile.ZipFile(archive_path, "w") as zf:
        zf.writestr("job_a/ranked_1.pdb", write_test_pdb(tmp_path / "ranked_1.pdb").read_text())
        zf.writestr("job_a/ranked_2.cif", write_test_mmcif(tmp_path / "ranked_2.cif").read_text())

    outdir = tmp_path / "out"
    seen_dssp_inputs = []

    def fake_run_dssp(structure_path, logger, *, display_name=None):
        seen_dssp_inputs.append(structure_path)
        assert display_name in {"ranked_1.pdb", "ranked_2.cif"}
        return {}, []

    monkeypatch.setattr(pipeline_module, "run_dssp", fake_run_dssp)
    monkeypatch.setattr(
        pipeline_module,
        "build_batch_analysis",
        lambda report_data, parsed_structures, **kwargs: SimpleNamespace(
            analysis=BatchAnalysis(
                tm_score=TMScoreAnalysis(enabled=False, available=False, cutoff=0.7),
                contact_consensus=ContactConsensusAnalysis(),
            ),
            warnings=[],
            partial_reasons=[],
        ),
    )

    exit_code = main(["--input", str(archive_path), "--outdir", str(outdir)])
    assert exit_code == 0

    assert len({path.name for path in seen_dssp_inputs}) == 2
    assert all(path.suffix == ".pdb" for path in seen_dssp_inputs)
    assert (outdir / "job_a" / "report.html").exists()


def test_cli_marks_partial_success_when_tm_clustering_is_unavailable_for_multi_model_job(tmp_path):
    archive_path = tmp_path / "models.zip"
    with zipfile.ZipFile(archive_path, "w") as zf:
        zf.writestr("job_a/ranked_1.pdb", write_test_pdb(tmp_path / "ranked_1.pdb").read_text())
        zf.writestr("job_a/ranked_2.pdb", write_test_pdb(tmp_path / "ranked_2.pdb").read_text())

    outdir = tmp_path / "out"
    exit_code = main(["--input", str(archive_path), "--outdir", str(outdir)])
    assert exit_code == 0

    job_dir = outdir / "job_a"
    payload = json.loads((job_dir / "tracks.json").read_text())
    assert payload["status"] == "partial_success"
    assert payload["batch_analysis"]["tm_score"]["available"] is False
    assert (job_dir / "csv" / "contact-consensus.csv").exists()
    assert not (job_dir / "csv" / "tm-score-matrix.csv").exists()
    summary_text = (job_dir / "summary.txt").read_text()
    assert "TM-score Clustering: Unavailable" in summary_text
    assert "US-align executable was not found on PATH" in summary_text


def test_cli_disable_tm_clustering_keeps_multi_model_job_success(tmp_path):
    archive_path = tmp_path / "models.zip"
    with zipfile.ZipFile(archive_path, "w") as zf:
        zf.writestr("job_a/ranked_1.pdb", write_test_pdb(tmp_path / "ranked_1.pdb").read_text())
        zf.writestr("job_a/ranked_2.pdb", write_test_pdb(tmp_path / "ranked_2.pdb").read_text())

    outdir = tmp_path / "out"
    exit_code = main(
        [
            "--input",
            str(archive_path),
            "--outdir",
            str(outdir),
            "--disable-tm-clustering",
        ]
    )
    assert exit_code == 0

    job_dir = outdir / "job_a"
    payload = json.loads((job_dir / "tracks.json").read_text())
    assert payload["status"] == "success"
    assert payload["batch_analysis"]["tm_score"]["enabled"] is False
    assert payload["batch_analysis"]["tm_score"]["available"] is False
    assert (job_dir / "csv" / "contact-consensus.csv").exists()
    assert not (job_dir / "csv" / "tm-score-matrix.csv").exists()
    summary_text = (job_dir / "summary.txt").read_text()
    assert "TM-score Clustering: Disabled" in summary_text


def test_ui_resources_are_available_from_package():
    ui_root = resources.files("openfoldpanel.UI")

    assert ui_root.joinpath("report.template.html").is_file()
    assert ui_root.joinpath("fonts/Times New Roman.ttf").is_file()
    assert ui_root.joinpath("styles/tokens.css").is_file()
    assert ui_root.joinpath("styles/figure.css").is_file()
    assert ui_root.joinpath("scripts/report.js").is_file()


def test_summarize_msa_database_path_maps_known_aliases_and_private_paths():
    assert summarize_msa_database_path(None) == "Not set"
    assert summarize_msa_database_path("pdbaa") == "PDBAA"
    assert summarize_msa_database_path("SWISSPROT") == "SWISSPROT"
    assert summarize_msa_database_path("PDBAA50") == "PDBAA50"
    assert summarize_msa_database_path("pdbaa70") == "PDBAA70"
    assert summarize_msa_database_path("pdbaa90") == "PDBAA90"
    assert summarize_msa_database_path("pdbaa95") == "PDBAA95"
    assert summarize_msa_database_path("./blastdb/swissprot_fasta/uniprot_sprot.fasta") == "SWISSPROT"
    assert summarize_msa_database_path("./A") == "A"
