from __future__ import annotations

import json
import tarfile
from importlib import resources

import pytest

import openfoldpanel.pipeline as pipeline_module
from openfoldpanel.cli import build_parser, main
from openfoldpanel.constants import ALLOWED_EVALUES, DEFAULT_EVALUE
from openfoldpanel.utils.text import summarize_msa_database_path
from tests.conftest import write_test_pdb


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
    assert not (job_dir / "panel.svg").exists()
    assert not (job_dir / "panel.png").exists()

    payload = json.loads((job_dir / "tracks.json").read_text())
    assert payload["default_reference_chain"] == "A"
    assert payload["status"] == "success"
    assert {panel["reference_chain"] for panel in payload["chain_panels"]} == {"A", "B"}
    assert {panel["job_name"] for panel in payload["chain_panels"]} == {"single_model"}
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
    assert 'data-report-page' in html_text
    assert 'data-active-chain-panel' in html_text
    assert 'data-figure-sheet' in html_text
    assert 'data-summary-grid' in html_text
    assert 'ofp-toolbar-summary' in html_text
    assert 'id="ofp-report-payload"' in html_text
    assert 'data-chain-templates' in html_text
    assert 'template id="ofp-chain-A"' in html_text
    assert 'template id="ofp-chain-B"' in html_text
    assert 'data-chain-figure="A"' in html_text
    assert 'data-chain-figure="B"' in html_text
    assert 'data-panel-width="' in html_text
    assert "--active-panel-width:" in html_text
    assert 'ofp-figure-wrap figure-wrap' in html_text
    assert 'ofp-figure-sheet figure-sheet' in html_text
    assert 'data-legend-deck' in html_text
    assert 'data-legend-card="structure"' in html_text
    assert 'data-legend-card="tracks"' in html_text
    assert 'data-legend-card="contacts"' in html_text
    assert 'data-legend-kind="strand"' in html_text
    assert 'data-legend-kind="helix"' in html_text
    assert 'data-legend-kind="turn"' in html_text
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
