"""Consented model storage: real file validation, hostile archives and no network by default."""

import hashlib
import io
import tarfile
from types import SimpleNamespace

import pytest

from dvbfixer import install_weights


def args(**kwargs):
    return SimpleNamespace(**{
        "checkpoint": None, "runtime": "protpardelle-mps", "download_weights": False,
        "accept_model_terms": False, "dry_run": False, **kwargs,
    })


def inputs(tmp_path, monkeypatch, *, cc94=True):
    requirements = {"checkpoint": {"name": "cc94_epoch3100.pth", "sha256": hashlib.sha256(b"weights").hexdigest()}}
    if cc94:
        requirements["config"] = hashlib.sha256(b"config").hexdigest()
    monkeypatch.setattr(install_weights, "asset_requirements", lambda *a: requirements)
    return requirements


def archive_bytes(entries):
    data = io.BytesIO()
    with tarfile.open(fileobj=data, mode="w:gz") as stream:
        for name, contents, linked in entries:
            member = tarfile.TarInfo(name)
            if linked:
                member.type = tarfile.SYMTYPE
                member.linkname = "/outside"
            else:
                member.size = len(contents)
            stream.addfile(member, io.BytesIO(contents))
    return data.getvalue()


def good_archive():
    return archive_bytes([
        ("protpardelle-1c/model_params/weights/cc94_epoch3100.pth", b"weights", False),
        ("protpardelle-1c/model_params/configs/cc94.yaml", b"config", False),
        ("../../outside", b"unwanted", False),
        ("protpardelle-1c/other-weights.pth", b"unwanted", False),
    ])


def fake_download(monkeypatch, content):
    calls = []

    def download(url, path, size):
        calls.append(url)
        path.write_bytes(content)

    monkeypatch.setattr(install_weights, "download", download)
    monkeypatch.setattr(install_weights, "CC94_MD5", hashlib.md5(content).hexdigest())
    return calls


def test_verified_cache_is_shared_between_cuda_and_mps_without_redownload(tmp_path, monkeypatch):
    inputs(tmp_path, monkeypatch)
    calls = fake_download(monkeypatch, good_archive())
    checkpoint = install_weights.obtain(tmp_path, "protpardelle-mps")
    assert checkpoint.read_bytes() == b"weights"
    assert checkpoint.parent.parent.joinpath("configs/cc94.yaml").read_bytes() == b"config"
    assert install_weights.obtain(tmp_path, "protpardelle-cuda") == checkpoint
    assert len(calls) == 1
    assert not (tmp_path / "outside").exists()
    assert not list((tmp_path / "models").glob(".weights-*"))


def test_bad_digest_publishes_nothing_and_cleans_partial_download(tmp_path, monkeypatch):
    inputs(tmp_path, monkeypatch)
    fake_download(monkeypatch, good_archive())
    monkeypatch.setattr(install_weights, "CC94_MD5", "0" * 32)
    with pytest.raises(ValueError, match="checksum mismatch"):
        install_weights.obtain(tmp_path, "protpardelle-mps")
    assert not (tmp_path / "models/protpardelle-cc94").exists()
    assert list((tmp_path / "models").iterdir()) == []


def test_bad_checkpoint_sha_never_enters_cache(tmp_path, monkeypatch):
    requirements = inputs(tmp_path, monkeypatch)
    requirements["checkpoint"]["sha256"] = "0" * 64
    fake_download(monkeypatch, good_archive())
    with pytest.raises(ValueError, match="SHA256"):
        install_weights.obtain(tmp_path, "protpardelle-mps")
    assert not (tmp_path / "models/protpardelle-cc94").exists()


def test_interrupted_download_never_becomes_a_cache_hit(tmp_path, monkeypatch):
    inputs(tmp_path, monkeypatch)

    def interrupted(url, path, size):
        path.write_bytes(b"partial")
        raise KeyboardInterrupt()

    monkeypatch.setattr(install_weights, "download", interrupted)
    with pytest.raises(KeyboardInterrupt):
        install_weights.obtain(tmp_path, "protpardelle-mps")
    assert list((tmp_path / "models").iterdir()) == []


@pytest.mark.parametrize("linked, duplicate", [(True, False), (False, True)])
def test_expected_archive_members_must_be_regular_and_unique(tmp_path, linked, duplicate):
    member = ("protpardelle-1c/model_params/weights/cc94_epoch3100.pth", b"weights", linked)
    entries = [member, member] if duplicate else [member]
    archive = tmp_path / "archive.gz"
    archive.write_bytes(archive_bytes(entries))
    with pytest.raises(ValueError, match="Invalid or duplicate"):
        install_weights.extract_cc94(archive, tmp_path / "stage")


def test_corrupt_existing_cache_is_preserved_not_overwritten(tmp_path, monkeypatch):
    inputs(tmp_path, monkeypatch)
    cache = tmp_path / "models/protpardelle-cc94/model_params/weights"
    cache.mkdir(parents=True)
    checkpoint = cache / "cc94_epoch3100.pth"
    checkpoint.write_bytes(b"user bytes")
    monkeypatch.setattr(install_weights, "download", lambda *a: pytest.fail("do not overwrite cache"))
    with pytest.raises(ValueError, match="SHA256"):
        install_weights.obtain(tmp_path, "protpardelle-mps")
    assert checkpoint.read_bytes() == b"user bytes"


def test_linked_cache_root_is_refused_without_writing(tmp_path):
    outside = tmp_path / "outside"
    outside.mkdir()
    (tmp_path / "models").symlink_to(outside)
    with pytest.raises(ValueError, match="symlink"):
        install_weights.obtain(tmp_path, "protpardelle-mps")
    assert list(outside.iterdir()) == []


def test_no_download_without_interactive_consent(tmp_path, monkeypatch):
    monkeypatch.setattr(install_weights, "download", lambda *a: pytest.fail("no consent"))
    install_weights.select(args(), tmp_path, interactive=False)
    with pytest.raises(ValueError, match="accept-model-terms"):
        install_weights.select(args(download_weights=True), tmp_path, interactive=False)
    assert not (tmp_path / "models").exists()


def test_dry_run_does_not_download_or_create_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(install_weights, "download", lambda *a: pytest.fail("preview only"))
    selected = args(download_weights=True, accept_model_terms=True, dry_run=True)
    install_weights.select(selected, tmp_path, interactive=False)
    assert selected.checkpoint is None
    assert not (tmp_path / "models").exists()


def test_interactive_local_file_and_deferral_do_not_download(tmp_path, monkeypatch):
    monkeypatch.setattr(install_weights, "download", lambda *a: pytest.fail("local/deferred only"))
    answers = iter(["2", "/my/model.pth", "3"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    local = args()
    install_weights.select(local, tmp_path, interactive=True)
    assert str(local.checkpoint) == "/my/model.pth"
    deferred = args()
    install_weights.select(deferred, tmp_path, interactive=True)
    assert deferred.checkpoint is None
    assert not (tmp_path / "models").exists()


def test_interactive_download_consent_precedes_transfer(tmp_path, monkeypatch, capsys):
    inputs(tmp_path, monkeypatch)
    content = good_archive()
    monkeypatch.setattr(install_weights, "CC94_MD5", hashlib.md5(content).hexdigest())
    monkeypatch.setattr("builtins.input", lambda _: "1")

    def transfer(url, path, size):
        output = capsys.readouterr().out
        assert "Terms: CC BY 4.0" in output
        assert "4.29 GB" in output
        assert "Accept displayed model terms" in output
        path.write_bytes(content)

    monkeypatch.setattr(install_weights, "download", transfer)
    selected = args()
    install_weights.select(selected, tmp_path, interactive=True)
    assert selected.checkpoint.read_bytes() == b"weights"


def test_protenix_download_uses_its_frozen_sha_and_no_other_assets(tmp_path, monkeypatch):
    requirements = inputs(tmp_path, monkeypatch, cc94=False)
    requirements["checkpoint"]["name"] = "protenix_base_default_v1.0.0.pt"
    calls = fake_download(monkeypatch, b"weights")
    checkpoint = install_weights.obtain(tmp_path, "protenix-cuda")
    assert checkpoint.name == "protenix_base_default_v1.0.0.pt"
    assert calls == [install_weights.PROTENIX_URL]
    assert not (tmp_path / "models/protenix-v1/common").exists()
    assert install_weights.obtain(tmp_path, "protenix-cuda") == checkpoint
    assert len(calls) == 1


def test_verified_cache_reuse_does_not_prompt_or_download(tmp_path, monkeypatch):
    inputs(tmp_path, monkeypatch)
    fake_download(monkeypatch, good_archive())
    checkpoint = install_weights.obtain(tmp_path, "protpardelle-mps")
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("verified weights do not need input"))
    monkeypatch.setattr(install_weights, "download", lambda *a: pytest.fail("verified weights do not need download"))
    selected = args(runtime="protpardelle-cuda")
    install_weights.select(selected, tmp_path, interactive=True)
    assert selected.checkpoint == checkpoint


def test_no_plaintext_redirect_is_accepted():
    with pytest.raises(ValueError, match="non-HTTPS"):
        install_weights.SecureRedirect().redirect_request(None, None, 302, "redirect", {}, "http://example.invalid/model")


@pytest.mark.parametrize("payload,size,match", [(b"abc", 4, "Incomplete"), (b"abcde", 4, "exceeds")])
def test_network_transfer_is_size_bounded(tmp_path, monkeypatch, payload, size, match):
    class Response(io.BytesIO):
        def geturl(self):
            return "https://example.invalid/model"

    opener = SimpleNamespace(open=lambda *a, **kw: Response(payload))
    monkeypatch.setattr(install_weights.urllib.request, "build_opener", lambda *a: opener)
    with pytest.raises(ValueError, match=match):
        install_weights.download("https://example.invalid/model", tmp_path / "download", size)
