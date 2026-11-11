import gzip
import hashlib
import io
import json
import tarfile

import pytest
from scripts import archive_mac_reconstruction_evidence as evidence


def test_archive_roundtrip_hashes_and_browser_cache_exclusion(tmp_path, monkeypatch):
    source = tmp_path / "tmp" / "study"
    source.mkdir(parents=True)
    (source / "native.pdb").write_bytes(b"original observations\n")
    profile = source / "chrome-live-test" / "profile"
    profile.mkdir(parents=True)
    (profile / "Cookies").write_text("do not publish")
    (source / ".DS_Store").write_bytes(b"host metadata")
    (source / "outside-link").symlink_to(tmp_path)
    monkeypatch.setattr(evidence, "ROOTS", ("study",))
    destination = tmp_path / "archive"
    evidence.archive(tmp_path, destination)
    manifest = json.loads((destination / "manifest.json").read_text())
    archive = manifest["archives"][0]
    assert len(archive["files"]) == 1
    assert len(manifest["excluded"]) == 3
    data = b"".join((destination / p["path"]).read_bytes() for p in archive["parts"])
    assert hashlib.sha256(data).hexdigest() == archive["parts"][0]["sha256"]
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as stream:
        members = stream.getmembers()
        assert [m.name for m in members] == ["tmp/study/native.pdb"]
        assert stream.extractfile(members[0]).read() == b"original observations\n"
    with pytest.raises(FileExistsError):
        evidence.archive(tmp_path, destination)


def test_parts_boundaries_preserve_single_compressed_stream(tmp_path):
    output = evidence.Parts(tmp_path, "small", limit=7)
    with gzip.GzipFile(filename="", fileobj=output, mode="wb", mtime=0) as stream:
        stream.write(b"roundtrip" * 100)
    output.close()
    assert len(output.paths) > 1
    assert all(p.stat().st_size <= 7 for p in output.paths)
    assert gzip.decompress(b"".join(p.read_bytes() for p in output.paths)) == b"roundtrip" * 100
