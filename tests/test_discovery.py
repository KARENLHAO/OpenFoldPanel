from __future__ import annotations

import logging
import tarfile
import zipfile

from openfoldpanel.extractors.archive import extract_archive
from openfoldpanel.extractors.discovery import (
    discover_jobs_from_extracted_root,
    discover_jobs_from_structure,
)
from openfoldpanel.utils.sorting import natural_sort_key
from tests.conftest import write_test_pdb


def test_discover_single_structure_job(tmp_path):
    structure = write_test_pdb(tmp_path / "model_1.pdb")
    jobs = discover_jobs_from_structure(structure)
    assert len(jobs) == 1
    assert jobs[0].name == "model_1"
    assert jobs[0].structure_files == [structure]


def test_discover_archive_root_files_as_single_job(tmp_path):
    archive_path = tmp_path / "models.zip"
    with zipfile.ZipFile(archive_path, "w") as zf:
        zf.writestr("ranked_10.pdb", write_test_pdb(tmp_path / "ranked_10.pdb").read_text())
        zf.writestr("ranked_2.pdb", write_test_pdb(tmp_path / "ranked_2.pdb").read_text())
        zf.writestr("notes.txt", "ignore me")

    extracted = tmp_path / "extracted"
    extract_archive(archive_path, extracted)
    jobs = discover_jobs_from_extracted_root(extracted, logging.getLogger("test"))
    assert len(jobs) == 1
    assert [path.name for path in jobs[0].structure_files] == ["ranked_2.pdb", "ranked_10.pdb"]
    assert jobs[0].ignored_files == ["notes.txt"]


def test_discover_archive_directories_as_jobs(tmp_path):
    archive_path = tmp_path / "batch.tar.gz"
    job_a = tmp_path / "job_a"
    job_b = tmp_path / "job_b"
    job_a.mkdir()
    job_b.mkdir()
    write_test_pdb(job_a / "ranked_1.pdb")
    write_test_pdb(job_b / "ranked_2.pdb")
    with tarfile.open(archive_path, "w:gz") as tf:
        tf.add(job_a, arcname="job_a")
        tf.add(job_b, arcname="job_b")

    extracted = tmp_path / "extracted_tar"
    extract_archive(archive_path, extracted)
    jobs = discover_jobs_from_extracted_root(extracted, logging.getLogger("test"))
    assert [job.name for job in jobs] == ["job_a", "job_b"]
    assert all(job.structure_files for job in jobs)


def test_natural_sort_key_orders_ranked_models():
    names = ["ranked_10.pdb", "ranked_2.pdb", "ranked_1.pdb"]
    assert sorted(names, key=natural_sort_key) == ["ranked_1.pdb", "ranked_2.pdb", "ranked_10.pdb"]
