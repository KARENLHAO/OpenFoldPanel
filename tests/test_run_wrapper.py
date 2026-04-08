from __future__ import annotations

import importlib.util
import tarfile
from pathlib import Path

from tests.conftest import write_test_pdb


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "run.py"


def load_run_module():
    spec = importlib.util.spec_from_file_location("run_wrapper", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_run_wrapper_smoke_copies_html_and_archives_results(tmp_path, monkeypatch):
    module = load_run_module()
    input_path = write_test_pdb(tmp_path / "single_model.pdb")
    monkeypatch.chdir(tmp_path)

    def fake_run_ext_cmder(command):
        assert "--input" in command
        output_dir = module.current_output_dir()
        job_dir = output_dir / "single_model"
        job_dir.mkdir(parents=True, exist_ok=True)
        (job_dir / "report.html").write_text("<html>report</html>", encoding="utf-8")
        (job_dir / "reference-chain-A.pdf").write_bytes(b"%PDF-1.4\n")
        return "OK"

    monkeypatch.setattr(module, "run_ext_cmder", fake_run_ext_cmder)

    exit_code = module.main(["--input", str(input_path)])

    assert exit_code == 0
    copied_html = tmp_path / module.HTML_OUTPUT_NAME
    archive_path = tmp_path / module.RESULT_ARCHIVE_NAME
    assert copied_html.exists()
    assert copied_html.read_text(encoding="utf-8") == "<html>report</html>"
    assert archive_path.exists()
    assert not module.current_output_dir().exists()

    with tarfile.open(archive_path, "r:gz") as handle:
        assert sorted(handle.getnames()) == [
            "single_model/reference-chain-A.pdf",
            "single_model/report.html",
        ]
