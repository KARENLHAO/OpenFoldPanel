from __future__ import annotations

import io
import tarfile
import zipfile

import pytest

from openfoldpanel.extractors.archive import extract_archive


def test_extract_zip_blocks_zip_slip(tmp_path):
    archive_path = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive_path, "w") as zf:
        zf.writestr("../escape.pdb", "malicious")

    with pytest.raises(ValueError):
        extract_archive(archive_path, tmp_path / "out")


def test_extract_tar_blocks_path_traversal(tmp_path):
    archive_path = tmp_path / "bad.tar.gz"
    with tarfile.open(archive_path, "w:gz") as tf:
        payload = b"malicious"
        info = tarfile.TarInfo("../escape.pdb")
        info.size = len(payload)
        tf.addfile(info, io.BytesIO(payload))

    with pytest.raises(ValueError):
        extract_archive(archive_path, tmp_path / "out")
