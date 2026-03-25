from __future__ import annotations

import importlib.util
import io
import json
from pathlib import Path

import pytest


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "blastdb" / "Build_pdbaa_clusters.py"


def load_build_pdbaa_clusters_module():
    spec = importlib.util.spec_from_file_location("build_pdbaa_clusters", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FakeUrlopenResponse:
    def __init__(self, payload: dict):
        self._buffer = io.StringIO(json.dumps(payload))

    def __enter__(self):
        return self._buffer

    def __exit__(self, exc_type, exc, tb):
        self._buffer.close()
        return False


def write_cluster_file(path: Path, tokens: list[str]) -> None:
    path.write_text("\n".join(tokens) + "\n", encoding="utf-8")


def write_pdb_seqres(path: Path, records: list[tuple[str, str]]) -> None:
    chunks: list[str] = []
    for header, sequence in records:
        chunks.append(f">{header}")
        chunks.append(sequence)
    path.write_text("\n".join(chunks) + "\n", encoding="utf-8")


@pytest.fixture
def module(tmp_path, monkeypatch):
    module = load_build_pdbaa_clusters_module()
    monkeypatch.setattr(module, "BUILD_DIR", tmp_path / "build")
    monkeypatch.setattr(module, "BUILD_BLASTDB_SCRIPT", tmp_path / "Build_blastdb.sh")
    monkeypatch.setattr(module.time, "sleep", lambda _seconds: None)
    return module


def test_iter_representative_entities_keeps_only_pdb_tokens_and_warns(tmp_path, module, capsys):
    cluster_path = tmp_path / "clusters-by-entity-95.txt"
    write_cluster_file(
        cluster_path,
        [
            "5B8C_2 5DK3_2 5GGS_1",
            "AF_AFP01670F1_1 AF_AFP01671F1_1",
            "MA_1234_1 5B8C_1",
            "NOT_A_VALID_TOKEN",
            "1QFW_3",
        ],
    )

    representatives = module.iter_representative_entities(cluster_path)

    assert representatives == ["5B8C_2", "1QFW_3"]
    stderr = capsys.readouterr().err
    assert "AF_AFP01670F1_1" in stderr
    assert "MA_1234_1" in stderr
    assert "NOT_A_VALID_TOKEN" in stderr


def test_load_pdb_seqres_index_keeps_only_protein_sequences_with_uppercase_entry_ids(tmp_path, module):
    seqres_path = tmp_path / "pdb_seqres.txt"
    write_pdb_seqres(
        seqres_path,
        [
            ("100d_A mol:na length:10 DNA", "CCGGCGCCGG"),
            ("5b8c_A mol:protein length:119 Pembrolizumab light chain variable region", "DIQMTQSPSS"),
            ("5b8c_D mol:protein length:119 Pembrolizumab light chain variable region", "TTTMTQSPSS"),
            ("1qfw_B mol:protein length:120 Example heavy chain", "EVQLVESGGG"),
        ],
    )

    seq_index = module.load_pdb_protein_chain_sequences(seqres_path)

    assert seq_index == {
        ("5B8C", "A"): "DIQMTQSPSS",
        ("5B8C", "D"): "TTTMTQSPSS",
        ("1QFW", "B"): "EVQLVESGGG",
    }


def test_build_identity_fasta_uses_cached_or_fetched_chain_mapping_and_local_seqres_sequences(
    tmp_path, module, monkeypatch
):
    cluster_dir = tmp_path / "pdb_cluster_src"
    cluster_dir.mkdir()
    write_cluster_file(cluster_dir / "clusters-by-entity-95.txt", ["5B8C_2 5DK3_2", "1QFW_3"])
    write_pdb_seqres(
        cluster_dir / "pdb_seqres.txt",
        [
            ("5b8c_A mol:protein length:119 Pembrolizumab light chain variable region", "AAAAAA"),
            ("5b8c_D mol:protein length:119 Pembrolizumab light chain variable region", "DDDDDD"),
            ("1qfw_B mol:protein length:120 Example heavy chain", "BBBBBB"),
            ("100d_A mol:na length:10 DNA", "CCGGCGCCGG"),
        ],
    )

    payloads = {
        "https://data.rcsb.org/rest/v1/core/polymer_entity/5B8C/2": {
            "rcsb_polymer_entity_container_identifiers": {
                "auth_asym_ids": ["A", "D", "G", "J"],
            }
        },
        "https://data.rcsb.org/rest/v1/core/polymer_entity/1QFW/3": {
            "rcsb_polymer_entity_container_identifiers": {
                "auth_asym_ids": ["B"],
            }
        },
    }
    seen_urls: list[str] = []

    def fake_urlopen(request, timeout=30):
        seen_urls.append(request.full_url)
        return FakeUrlopenResponse(payloads[request.full_url])

    monkeypatch.setattr(module, "urlopen", fake_urlopen)

    fasta_path, blast_prefix = module.build_identity_fasta(95, cluster_dir, sleep_seconds=0, force=False)

    assert fasta_path == module.BUILD_DIR / "pdbaa95" / "pdbaa95.fasta"
    assert blast_prefix == module.BUILD_DIR / "pdbaa95" / "pdbaa95"
    assert fasta_path.read_text(encoding="utf-8") == (
        ">pdb|5B8C|A\nAAAAAA\n"
        ">pdb|1QFW|B\nBBBBBB\n"
    )
    assert seen_urls == [
        "https://data.rcsb.org/rest/v1/core/polymer_entity/5B8C/2",
        "https://data.rcsb.org/rest/v1/core/polymer_entity/1QFW/3",
    ]

    cache_path = module.BUILD_DIR / "pdbaa_cache" / "entity_chain_map.json"
    assert json.loads(cache_path.read_text(encoding="utf-8")) == {
        "5B8C_2": {"entry_id": "5B8C", "entity_id": "2", "chain_id": "A"},
        "1QFW_3": {"entry_id": "1QFW", "entity_id": "3", "chain_id": "B"},
    }

    monkeypatch.setattr(module, "urlopen", lambda *_args, **_kwargs: pytest.fail("urlopen should not be called"))

    cached_fasta_path, _cached_blast_prefix = module.build_identity_fasta(95, cluster_dir, sleep_seconds=0, force=True)

    assert cached_fasta_path.read_text(encoding="utf-8") == fasta_path.read_text(encoding="utf-8")


def test_build_identity_fasta_fails_when_chain_mapping_cannot_be_fetched(tmp_path, module, monkeypatch):
    cluster_dir = tmp_path / "pdb_cluster_src"
    cluster_dir.mkdir()
    write_cluster_file(cluster_dir / "clusters-by-entity-50.txt", ["5B8C_2"])
    write_pdb_seqres(
        cluster_dir / "pdb_seqres.txt",
        [("5b8c_A mol:protein length:119 Pembrolizumab light chain variable region", "AAAAAA")],
    )

    def fake_urlopen(_request, timeout=30):
        raise module.URLError("network down")

    monkeypatch.setattr(module, "urlopen", fake_urlopen)

    with pytest.raises(SystemExit):
        module.build_identity_fasta(50, cluster_dir, sleep_seconds=0, force=False)

    assert not (module.BUILD_DIR / "pdbaa50" / "pdbaa50.fasta").exists()


def test_build_identity_fasta_skips_mapped_chain_without_local_protein_sequence(tmp_path, module, monkeypatch, capsys):
    cluster_dir = tmp_path / "pdb_cluster_src"
    cluster_dir.mkdir()
    write_cluster_file(cluster_dir / "clusters-by-entity-70.txt", ["5B8C_2", "1QFW_3"])
    write_pdb_seqres(
        cluster_dir / "pdb_seqres.txt",
        [("1qfw_B mol:protein length:120 Example heavy chain", "BBBBBB")],
    )

    payloads = {
        "https://data.rcsb.org/rest/v1/core/polymer_entity/5B8C/2": {
            "rcsb_polymer_entity_container_identifiers": {
                "auth_asym_ids": ["A", "D"],
            }
        },
        "https://data.rcsb.org/rest/v1/core/polymer_entity/1QFW/3": {
            "rcsb_polymer_entity_container_identifiers": {
                "auth_asym_ids": ["B"],
            }
        },
    }

    def fake_urlopen(request, timeout=30):
        return FakeUrlopenResponse(payloads[request.full_url])

    monkeypatch.setattr(module, "urlopen", fake_urlopen)

    fasta_path, _blast_prefix = module.build_identity_fasta(70, cluster_dir, sleep_seconds=0, force=False)

    assert fasta_path.read_text(encoding="utf-8") == ">pdb|1QFW|B\nBBBBBB\n"
    assert "5B8C_2" in capsys.readouterr().err


def test_build_identity_fasta_resumes_from_existing_temp_without_rewriting_completed_records(
    tmp_path, module, monkeypatch
):
    cluster_dir = tmp_path / "pdb_cluster_src"
    cluster_dir.mkdir()
    write_cluster_file(cluster_dir / "clusters-by-entity-50.txt", ["5B8C_2", "1QFW_3", "2ABC_1"])
    write_pdb_seqres(
        cluster_dir / "pdb_seqres.txt",
        [
            ("5b8c_A mol:protein length:119 Pembrolizumab light chain variable region", "AAAAAA"),
            ("1qfw_B mol:protein length:120 Example heavy chain", "BBBBBB"),
            ("2abc_C mol:protein length:88 Example protein", "CCCCCC"),
        ],
    )

    first_pass_payloads = {
        "https://data.rcsb.org/rest/v1/core/polymer_entity/5B8C/2": {
            "rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": ["A"]}
        },
        "https://data.rcsb.org/rest/v1/core/polymer_entity/1QFW/3": {
            "rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": ["B"]}
        },
    }
    second_pass_payloads = {
        "https://data.rcsb.org/rest/v1/core/polymer_entity/2ABC/1": {
            "rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": ["C"]}
        }
    }
    first_pass_seen: list[str] = []
    second_pass_seen: list[str] = []

    def first_pass_urlopen(request, timeout=30):
        first_pass_seen.append(request.full_url)
        payload = first_pass_payloads.get(request.full_url)
        if payload is None:
            raise module.URLError("connection dropped")
        return FakeUrlopenResponse(payload)

    monkeypatch.setattr(module, "urlopen", first_pass_urlopen)

    with pytest.raises(SystemExit):
        module.build_identity_fasta(50, cluster_dir, sleep_seconds=0, force=False)

    temp_path = module.BUILD_DIR / "pdbaa50" / "pdbaa50.fasta.tmp"
    assert temp_path.read_text(encoding="utf-8") == (
        ">pdb|5B8C|A\nAAAAAA\n"
        ">pdb|1QFW|B\nBBBBBB\n"
    )
    assert first_pass_seen == [
        "https://data.rcsb.org/rest/v1/core/polymer_entity/5B8C/2",
        "https://data.rcsb.org/rest/v1/core/polymer_entity/1QFW/3",
        "https://data.rcsb.org/rest/v1/core/polymer_entity/2ABC/1",
    ]

    def second_pass_urlopen(request, timeout=30):
        second_pass_seen.append(request.full_url)
        payload = second_pass_payloads.get(request.full_url)
        if payload is None:
            pytest.fail(f"unexpected resumed fetch: {request.full_url}")
        return FakeUrlopenResponse(payload)

    monkeypatch.setattr(module, "urlopen", second_pass_urlopen)

    fasta_path, _blast_prefix = module.build_identity_fasta(50, cluster_dir, sleep_seconds=0, force=False)

    assert fasta_path.read_text(encoding="utf-8") == (
        ">pdb|5B8C|A\nAAAAAA\n"
        ">pdb|1QFW|B\nBBBBBB\n"
        ">pdb|2ABC|C\nCCCCCC\n"
    )
    assert second_pass_seen == [
        "https://data.rcsb.org/rest/v1/core/polymer_entity/2ABC/1",
    ]
    assert not temp_path.exists()


def test_build_identity_fasta_requires_pdb_seqres_file(tmp_path, module):
    cluster_dir = tmp_path / "pdb_cluster_src"
    cluster_dir.mkdir()
    write_cluster_file(cluster_dir / "clusters-by-entity-90.txt", ["5B8C_2"])

    with pytest.raises(SystemExit):
        module.build_identity_fasta(90, cluster_dir, sleep_seconds=0, force=False)


def test_build_blast_database_reuses_existing_script_contract(tmp_path, module, monkeypatch):
    fasta_path = tmp_path / "pdbaa95.fasta"
    fasta_path.write_text(">pdb|5B8C|A\nAAAAAA\n", encoding="utf-8")
    blast_prefix = tmp_path / "pdbaa95"

    recorded: dict[str, object] = {}

    def fake_run(command, check):
        recorded["command"] = command
        recorded["check"] = check

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    module.build_blast_database(95, fasta_path, blast_prefix, force=True)

    assert recorded == {
        "command": [
            "bash",
            str(module.BUILD_BLASTDB_SCRIPT),
            "--input",
            str(fasta_path),
            "--out-prefix",
            str(blast_prefix),
            "--title",
            "PDBAA95",
            "--force",
        ],
        "check": True,
    }
