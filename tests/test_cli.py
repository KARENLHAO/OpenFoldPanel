from __future__ import annotations

import json
import tarfile

import pytest

import openfoldpanel.pipeline as pipeline_module
from openfoldpanel.cli import build_parser, main
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
    assert payload["chain_panels"][0]["models"]
    assert any(model["display_name"] == "Single Model / 链 A" for model in payload["chain_panels"][0]["models"])

    html_text = (job_dir / "report.html").read_text()
    assert "链 A" in html_text
    assert "链 B" in html_text
    assert "查询序列" in html_text
    assert "可及性" in html_text
    assert "疏水性" in html_text
    assert "FoldScript 风格图板" in html_text
    assert "参考链选择" in html_text
    assert "链摘要" in html_text
    assert "图例" in html_text
    assert "single_model_A" in html_text
    assert 'data-chain-select' in html_text
    assert 'data-report-page' in html_text
    assert 'data-panel-width="' in html_text
    assert "--active-panel-width:" in html_text
    assert 'class="figure-wrap"' in html_text
    assert 'class="figure-sheet"' in html_text
    assert "导出当前链 PDF" not in html_text
    assert "下载 Tracks JSON" not in html_text
    assert "查看摘要" not in html_text
    assert "查看日志" not in html_text
    assert "术语说明" not in html_text
    assert "研究协作报告" not in html_text


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


def test_cli_parser_accepts_msa_display_rows():
    parser = build_parser()
    args = parser.parse_args(["--input", "demo.pdb", "--outdir", "out", "--msa-display-rows", "3"])
    assert args.msa_display_rows == 3


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
