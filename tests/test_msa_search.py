from __future__ import annotations

import logging
from pathlib import Path

import pytest

from openfoldpanel.features import msa_search
from openfoldpanel.utils.subprocess import CommandResult


def test_blast_search_accepts_fasta_msa_db_and_restores_full_uniprot_sequence(monkeypatch, tmp_path, caplog):
    database_fasta = tmp_path / "uniprot_sprot.fasta"
    database_fasta.write_text(
        (
            ">sp|Q15116|PDCD1_HUMAN RecName: Full=Programmed cell death protein 1;\n"
            "PWNPPTFSPALL\n"
            ">sp|Q99999|OTHER_MOUSE RecName: Full=Other protein;\n"
            "LMNOP\n"
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        msa_search,
        "which",
        lambda executable: "/usr/bin/blastp" if executable == "blastp" else "/usr/bin/makeblastdb" if executable == "makeblastdb" else None,
    )

    def fake_run_command(command, **_kwargs):
        if command[0] == "makeblastdb":
            assert command[1:5] == ["-in", str(database_fasta), "-dbtype", "prot"]
            assert "-parse_seqids" in command
            return CommandResult(command=list(command), returncode=0, stdout="", stderr="")
        if command[0] == "blastp":
            assert "-evalue" in command
            assert command[command.index("-evalue") + 1] == "1e-6"
            assert command[-1] == "6 sseqid sallseqid stitle sseq"
            db_prefix = Path(command[4])
            assert db_prefix != database_fasta
            assert db_prefix.parent.name == "_msa_blastdb"
            return CommandResult(
                command=list(command),
                returncode=0,
                stdout="sp|Q15116.3|\tsp|Q15116.3|\tRecName: Full=Programmed cell death protein 1;\tWNPPTFSPALL\n",
                stderr="",
            )
        raise AssertionError(f"unexpected command: {command}")

    monkeypatch.setattr(msa_search, "run_command", fake_run_command)
    caplog.set_level(logging.INFO)

    rows, warnings = msa_search.search_homologs(
        tmp_path / "query.fasta",
        database_fasta,
        max_homologs_displayed=5,
        evalue="1e-6",
        workdir=tmp_path / "chain_a",
        logger=logging.getLogger("test"),
    )

    assert warnings == []
    assert rows == [("sp|Q15116|PDCD1_HUMAN", "PWNPPTFSPALL")]
    assert "BLAST evalue threshold: 1e-6" in caplog.text


def test_search_warns_when_makeblastdb_is_missing_for_fasta_input(monkeypatch, tmp_path):
    database_fasta = tmp_path / "uniprot_sprot.fasta"
    database_fasta.write_text(
        ">sp|P01674|KV3AM_MOUSE RecName: Full=Ig kappa chain V-III region MOPC 63;\nACDEF\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(msa_search, "which", lambda executable: "/usr/bin/blastp" if executable == "blastp" else None)

    rows, warnings = msa_search.search_homologs(
        tmp_path / "query.fasta",
        database_fasta,
        max_homologs_displayed=5,
        evalue="1e-6",
        workdir=tmp_path / "chain_a",
        logger=logging.getLogger("test"),
    )

    assert rows == []
    assert warnings == [
        "makeblastdb was not found; FASTA MSA database input could not be prepared for blastp."
    ]


def test_blast_search_uses_structured_identifier_when_present_for_regular_database(monkeypatch, tmp_path):
    monkeypatch.setattr(msa_search, "which", lambda executable: "/usr/bin/blastp" if executable == "blastp" else None)

    def fake_run_command(command, **_kwargs):
        if command[0] == "blastp":
            assert "-evalue" in command
            assert command[command.index("-evalue") + 1] == "1e-7"
            assert command[-1] == "6 sseqid sallseqid stitle sseq"
            return CommandResult(
                command=list(command),
                returncode=0,
                stdout=(
                    "sp|P01674.1|\tsp|P01674|KV3AM_MOUSE;sp|P01674.1|\tRecName: Full=Ig kappa chain V-III region MOPC 63;\tACD-EF\n"
                    "pdb|5B8C|A\tpdb|5B8C|A\tChain A, Example protein\tGHI-JK\n"
                    "gnl|demo|entry-1\tgnl|demo|entry-1\tExample database protein\tLMN-OP\n"
                ),
                stderr="",
            )
        raise AssertionError("blastdbcmd should not be needed when the structured identifier is already present")

    monkeypatch.setattr(msa_search, "run_command", fake_run_command)

    rows, warnings = msa_search.search_homologs(
        tmp_path / "query.fasta",
        Path("/tmp/demo_db"),
        max_homologs_displayed=5,
        evalue="1e-7",
        workdir=tmp_path,
        logger=logging.getLogger("test"),
    )

    assert warnings == [
        "blastp header recovery failed for 1 homolog hit(s); raw identifiers were used as fallback.",
        "blastp full-length sequence recovery fell back to aligned fragments for 3 homolog hit(s).",
    ]
    assert rows == [
        ("sp|P01674|KV3AM_MOUSE", "ACDEF"),
        ("pdb|5B8C|A", "GHIJK"),
        ("gnl|demo|entry-1", "LMNOP"),
    ]


def test_blast_search_recovers_full_sequence_via_blastdbcmd(monkeypatch, tmp_path):
    monkeypatch.setattr(
        msa_search,
        "which",
        lambda executable: "/usr/bin/blastp"
        if executable == "blastp"
        else "/usr/bin/blastdbcmd"
        if executable == "blastdbcmd"
        else None,
    )

    def fake_run_command(command, **_kwargs):
        if command[0] == "blastp":
            return CommandResult(
                command=list(command),
                returncode=0,
                stdout="sp|Q15116.3|\tsp|Q15116.3|\tRecName: Full=Programmed cell death protein 1;\tWNPPTFSPALL\n",
                stderr="",
            )
        if command[0] == "blastdbcmd":
            return CommandResult(
                command=list(command),
                returncode=0,
                stdout=">sp|Q15116|PDCD1_HUMAN RecName: Full=Programmed cell death protein 1;\nPWNPPTFSPALL\n",
                stderr="",
            )
        raise AssertionError(f"unexpected command: {command}")

    monkeypatch.setattr(msa_search, "run_command", fake_run_command)

    rows, warnings = msa_search.search_homologs(
        tmp_path / "query.fasta",
        Path("/tmp/demo_db"),
        max_homologs_displayed=5,
        evalue="1e-6",
        workdir=tmp_path,
        logger=logging.getLogger("test"),
    )

    assert warnings == []
    assert rows == [("sp|Q15116|PDCD1_HUMAN", "PWNPPTFSPALL")]


def test_mmseqs_search_uses_fasta_header_token_and_full_sequence(monkeypatch, tmp_path, caplog):
    database_fasta = tmp_path / "uniprot_sprot.fasta"
    database_fasta.write_text(
        (
            ">sp|Q15116|PDCD1_HUMAN RecName: Full=Programmed cell death protein 1;\n"
            "PWNPPTFSPALL\n"
            ">pdb|5B8C|A Chain A, Example protein\n"
            "MGHIJK\n"
        ),
        encoding="utf-8",
    )
    output_path = tmp_path / "mmseqs_hits.tsv"
    output_path.write_text(
        (
            "sp|Q15116.3|\tRecName: Full=Programmed cell death protein 1;\tWNPPTFSPALL\n"
            "pdb|5B8C|A\tpdb|5B8C|A Chain A, Example protein\tGHI-JK\n"
            "gnl|demo|entry-1\tExample database protein\tLMN-OP\n"
        ),
        encoding="utf-8",
    )

    def fake_which(executable: str):
        if executable == "blastp":
            return None
        if executable == "mmseqs":
            return "/usr/bin/mmseqs"
        return None

    def fake_run_command(command, **_kwargs):
        assert command[:2] == ["mmseqs", "easy-search"]
        assert "target,theader,alntseq" in command
        assert "-e" in command
        assert command[command.index("-e") + 1] == "1e-8"
        return CommandResult(command=list(command), returncode=0, stdout="", stderr="")

    monkeypatch.setattr(msa_search, "which", fake_which)
    monkeypatch.setattr(msa_search, "run_command", fake_run_command)
    caplog.set_level(logging.INFO)

    rows, warnings = msa_search.search_homologs(
        tmp_path / "query.fasta",
        database_fasta,
        max_homologs_displayed=5,
        evalue="1e-8",
        workdir=tmp_path,
        logger=logging.getLogger("test"),
    )

    assert warnings == [
        "mmseqs header recovery failed for 1 homolog hit(s); raw identifiers were used as fallback.",
        "mmseqs full-length sequence recovery fell back to aligned fragments for 1 homolog hit(s).",
    ]
    assert rows == [
        ("sp|Q15116|PDCD1_HUMAN", "PWNPPTFSPALL"),
        ("pdb|5B8C|A", "MGHIJK"),
        ("gnl|demo|entry-1", "LMNOP"),
    ]
    assert "MMseqs evalue threshold: 1e-8" in caplog.text


def test_search_falls_back_to_raw_identifier_when_header_recovery_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(
        msa_search,
        "which",
        lambda executable: "/usr/bin/blastp" if executable == "blastp" else "/usr/bin/blastdbcmd" if executable == "blastdbcmd" else None,
    )

    def fake_run_command(command, **_kwargs):
        if command[0] == "blastdbcmd":
            return CommandResult(command=list(command), returncode=0, stdout=">gnl|demo|entry-1 Example database protein\nACDEF\n", stderr="")
        return CommandResult(
            command=list(command),
            returncode=0,
            stdout="gnl|demo|entry-1\tRecName: Full=Example protein;\tACD-EF\n",
            stderr="",
        )

    monkeypatch.setattr(msa_search, "run_command", fake_run_command)

    rows, warnings = msa_search.search_homologs(
        tmp_path / "query.fasta",
        Path("/tmp/demo_db"),
        max_homologs_displayed=5,
        evalue="1e-6",
        workdir=tmp_path,
        logger=logging.getLogger("test"),
    )

    assert rows == [("gnl|demo|entry-1", "ACDEF")]
    assert warnings == [
        "blastp header recovery failed for 1 homolog hit(s); raw identifiers were used as fallback."
    ]


def test_search_rejects_invalid_evalue_before_running_tools(tmp_path):
    with pytest.raises(ValueError, match="evalue must be one of:"):
        msa_search.search_homologs(
            tmp_path / "query.fasta",
            Path("/tmp/demo_db"),
            max_homologs_displayed=5,
            evalue="0.001",
            workdir=tmp_path,
            logger=logging.getLogger("test"),
        )
