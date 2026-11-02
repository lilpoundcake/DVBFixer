"""Read-only evidence collection, explicit unknowns and containment safeguards."""

import ast
import importlib.util
import json
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("license_evidence", ROOT / "scripts/collect_license_evidence.py")
assert SPEC is not None and SPEC.loader is not None
evidence = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evidence)


def test_urls_drop_credentials_and_tokens():
    assert evidence.public_url("https://user:secret@example.org/pkg?token=secret#key") == "https://example.org/pkg"
    assert evidence.public_url("https://example.org/t/secret/linux/pkg") == "redacted-token-channel"
    assert evidence.public_url("file:///private/source") == "local-or-non-http-source"


def test_prefix_metadata_and_license_capture_are_not_approval(tmp_path):
    prefix = tmp_path / "prefix"
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "LICENSE").write_text("example license")
    (prefix / "conda-meta/example.json").write_text(json.dumps({
        "name": "example", "version": "1", "build": "0", "license": "MIT",
        "files": ["LICENSE", "COPYING.missing", "lib/native.dylib"],
        "url": "https://user:password@example.org/pkg?auth=secret",
    }))
    info = prefix / "lib/python3.12/site-packages/example-1.dist-info"
    info.mkdir(parents=True)
    (info / "METADATA").write_text("Name: example\nVersion: 1\nLicense-Expression: MIT\nRequires-Dist: other>=1\n\n")
    (info / "LICENSE").write_text("example license")
    (info / "RECORD").write_text("supplied,not-verified,hashes\n")
    (info / "direct_url.json").write_text('{"url":"file:///private/source"}')
    (prefix / "lib/python3.1").symlink_to("python3.12", target_is_directory=True)
    store = evidence.EvidenceStore(tmp_path / "out")
    result = evidence.collect_prefix(prefix, store)
    conda = result["conda"][0]
    assert conda["compatibility_review"] == "unknown"
    assert conda["license_review"] == "captured-not-reviewed"
    assert conda["missing_legal_files"] == ["COPYING.missing"]
    assert conda["native_files"] == ["lib/native.dylib"]
    assert conda["url"] == "https://example.org/pkg"
    assert result["python"][0]["direct_url"] == "local-or-non-http-source"
    assert result["python"][0]["requires_dist"] == ["other>=1"]
    assert len(result["python"]) == 1
    assert result["site_package_aliases"] == {
        "lib/python3.12/site-packages": ["lib/python3.1/site-packages", "lib/python3.12/site-packages"]}
    assert len(list((tmp_path / "out/blobs").iterdir())) == 3


def test_escaping_paths_fail_closed(tmp_path):
    with pytest.raises(ValueError, match="escapes"):
        evidence.contained_file(tmp_path, "../outside")


def test_authors_modules_and_bytecode_are_not_license_texts():
    assert evidence.legal_name(Path("AUTHORS.txt"))
    assert evidence.legal_name(Path("COPYING.LESSER"))
    assert not evidence.legal_name(Path("authors.py"))
    assert not evidence.legal_name(Path("authors.cpython-313.pyc"))
    assert not evidence.legal_name(Path("license_key.txt"))
    assert evidence.legal_name(Path("GNU_LGPL_v3"))


def test_evidence_conflicts_are_not_overwritten(tmp_path):
    source = tmp_path / "LICENSE"
    source.write_text("text")
    store = evidence.EvidenceStore(tmp_path / "out")
    result = store.capture(source, "test")
    (tmp_path / "out" / result["blob"]).write_text("tampered")
    with pytest.raises(ValueError, match="conflicting"):
        store.capture(source, "test")


def test_reports_are_exclusive(tmp_path):
    store = evidence.EvidenceStore(tmp_path / "out")
    evidence.write_report(store.output, {"status": "unknown"})
    with pytest.raises(FileExistsError):
        evidence.write_report(store.output, {"status": "approved"})


def test_npm_mismatch_and_missing_packages_stay_visible(tmp_path):
    root = tmp_path / "gui"
    package = root / "node_modules/example"
    package.mkdir(parents=True)
    (root / "package-lock.json").write_text(json.dumps({"packages": {
        "": {}, "node_modules/example": {"version": "1"},
        "node_modules/platform-only": {"version": "1"},
    }}))
    (package / "package.json").write_text('{"name":"example","version":"2","license":"MIT"}')
    (package / "LICENSE").write_text("example legal text")
    result = evidence.collect_npm(root, evidence.EvidenceStore(tmp_path / "out"))
    assert result["missing_lock_entries"] == ["node_modules/platform-only"]
    entry = result["packages"][0]
    assert not entry["matches_locked_version"]
    assert entry["compatibility_review"] == "unknown"
    assert entry["distribution_inclusion"] == "unverified"


def test_archive_hashes_payload_without_extracting(tmp_path):
    path = tmp_path / "example.whl"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("example/__init__.py", "# code")
        archive.writestr("share/notices/LICENSE.txt", "notice")
        archive.writestr("weights.ckpt", "not an actual model")
    result = evidence.inspect_archive(path)
    assert len(result["entries"]) == 3
    assert result["legal_paths"] == ["share/notices/LICENSE.txt"]
    assert result["restricted_or_sensitive_name_findings"] == ["weights.ckpt"]
    assert result["compatibility_review"] == "unknown"
    assert not (tmp_path / "example").exists()


def test_archive_rejects_traversal(tmp_path):
    path = tmp_path / "bad.whl"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("../escape", "data")
    with pytest.raises(ValueError, match="Unsafe"):
        evidence.inspect_archive(path)


def test_protpardelle_source_and_checkpoint_terms_are_separate_without_loading_ml():
    tree = ast.parse((ROOT / "deploy/protpardelle-1c/checkpoint_gap_smoke.py").read_text())
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name) and node.func.id == "BackendProvenance"]
    assert len(calls) == 1
    values = {keyword.arg: keyword.value for keyword in calls[0].keywords}
    assert ast.literal_eval(values["source_license"]) == "MIT"
    assert ast.literal_eval(values["checkpoint_license"]) == "CC-BY-4.0"
    assert isinstance(values["checkpoint_sha256"], ast.Name)
    assert values["checkpoint_sha256"].id == "checkpoint_sha"


def test_source_fetch_failures_remain_explicit(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "scripts"))
    import fetch_license_sources

    def fail(*args, **kwargs):
        raise OSError("unavailable")

    monkeypatch.setattr(fetch_license_sources, "urlopen", fail)
    record, data = fetch_license_sources.fetch_source({"id": "test", "url": "https://example.org/LICENSE"})
    assert record["status"] == "unavailable"
    assert record["error_type"] == "OSError"
    assert data is None
    with pytest.raises(ValueError, match="public"):
        fetch_license_sources.fetch_source({"id": "bad", "url": "https://user:secret@example.org/LICENSE"})


def test_registered_source_digest_is_checked(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "scripts"))
    import fetch_license_sources

    class Response:
        url = "https://example.org/source.tar.gz"
        headers = {"Content-Type": "application/gzip"}

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self, limit):
            return b"bounded source"

    monkeypatch.setattr(fetch_license_sources, "urlopen", lambda *a, **kw: Response())
    source = {"id": "test", "url": Response.url, "expected_sha256": evidence.digest(b"bounded source")}
    record, data = fetch_license_sources.fetch_source(source)
    assert record["status"] == "captured-not-reviewed" and data == b"bounded source"
    record, data = fetch_license_sources.fetch_source({**source, "url": Response.url + "?filename=source.tar.gz"})
    assert record["status"] == "captured-not-reviewed"
    with pytest.raises(ValueError, match="public"):
        fetch_license_sources.fetch_source({**source, "url": Response.url + "?token=secret"})
    record, data = fetch_license_sources.fetch_source({**source, "expected_sha256": "0" * 64})
    assert record["status"] == "unavailable" and data is None


def test_python_record_owned_bundled_licenses_are_not_missed(tmp_path):
    prefix = tmp_path / "prefix"
    (prefix / "conda-meta").mkdir(parents=True)
    site = prefix / "lib/python3.12/site-packages"
    info = site / "example-1.dist-info"
    info.mkdir(parents=True)
    (info / "METADATA").write_text("Name: example\nVersion: 1\n\n")
    owned = site / "example/bundled"
    owned.mkdir(parents=True)
    (owned / "LICENSE").write_text("separate bundled-component terms")
    (owned / "license_key.txt").write_text("must not be captured")
    (info / "RECORD").write_text("example/bundled/LICENSE,,\nexample/bundled/license_key.txt,,\nexample/native.dylib,,\n")
    result = evidence.collect_prefix(prefix, evidence.EvidenceStore(tmp_path / "out"))
    package = result["python"][0]
    assert len(package["legal_evidence"]) == 1
    assert package["legal_evidence"][0]["origin"].endswith("example/bundled/LICENSE")
    assert package["native_file_candidates"] == ["lib/python3.12/site-packages/example/native.dylib"]


def test_forcefield_comparison_is_content_evidence_not_clearance(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "scripts"))
    import compare_forcefield_archive

    result = compare_forcefield_archive.compare_manifests(
        {"same": "a", "changed": "b", "local": "c"},
        {"same": "a", "changed": "x", "upstream": "d"},
    )
    assert result["matching"] == ["same"]
    assert result["modified"] == [{"path": "changed", "local_sha256": "b", "upstream_sha256": "x"}]
    assert result["only_local"] == ["local"]
    assert result["only_upstream"] == ["upstream"]
    with pytest.raises(ValueError, match="Both"):
        compare_forcefield_archive.compare_manifests({}, {"a": "b"})


def test_forcefield_archive_prefix_and_identity_are_explicit(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "scripts"))
    import compare_forcefield_archive

    local = tmp_path / "FF/example.ff"
    local.mkdir(parents=True)
    (local / "a.itp").write_bytes(b"parameter")
    (local / "untracked.itp").write_bytes(b"user work")
    archive = tmp_path / "reference.zip"
    with zipfile.ZipFile(archive, "w") as stream:
        stream.writestr("example.ff/a.itp", b"parameter")
        stream.writestr("LICENSE", "reference-only terms")

    class GitFiles:
        stdout = b"FF/example.ff/a.itp\0"

    monkeypatch.setattr(compare_forcefield_archive.subprocess, "run", lambda *a, **kw: GitFiles())
    sha = evidence.digest(archive.read_bytes())
    result = compare_forcefield_archive.build_report(tmp_path, archive, local, "example.ff", sha)
    assert result["matching"] == ["a.itp"]
    assert result["outside_archive_prefix"] == ["LICENSE"]
    assert result["only_local"] == []
    assert result["redistribution_clearance"] == "unknown"
    with pytest.raises(ValueError, match="predeclared"):
        compare_forcefield_archive.build_report(tmp_path, archive, local, "example.ff", "0" * 64)
    with pytest.raises(ValueError, match="safe"):
        compare_forcefield_archive.build_report(tmp_path, archive, local, "../escape", sha)


def test_charmm_comparison_binds_current_tracked_files():
    import hashlib

    report = json.loads((ROOT / "docs/provenance/charmm36-jul2022-comparison.json").read_text())
    assert report["redistribution_clearance"] == "unknown"
    assert len(report["matching"]) == 59
    assert {e["path"] for e in report["modified"]} == {"atomtypes.atp", "ffnonbonded.itp", "ions.itp"}
    assert report["only_upstream"] == ["watermodels.dat"]
    for relative, expected in report["local_manifest"].items():
        assert hashlib.sha256((ROOT / report["local_tree"] / relative).read_bytes()).hexdigest() == expected


def test_evidence_index_refuses_incomplete_or_corrupt_captures(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "scripts"))
    import summarize_license_evidence

    with pytest.raises(ValueError, match="missing"):
        summarize_license_evidence.summarize(tmp_path, tmp_path)
    for label in summarize_license_evidence.PREFIX_LABELS:
        path = tmp_path / "installed-v4" / label
        path.mkdir(parents=True)
        (path / "inventory.json").write_text(json.dumps({
            "kind": "installed-metadata-and-legal-evidence", "prefix_label": label,
            "conda": [], "python": [],
        }))
    for label in ("npm-v3", "upstream-v2", "build/inspection-v3-wheel", "build/inspection-v3-sdist"):
        path = tmp_path / label
        path.mkdir(parents=True)
        (path / "inventory.json").write_text('{"kind":"test-input"}')
    result = summarize_license_evidence.summarize(tmp_path, tmp_path)
    assert len(result["reports"]) == 12
    path = tmp_path / "npm-v3/inventory.json"
    path.write_text(json.dumps({"kind": "test", "evidence": {"blob": "blob", "sha256": "0" * 64}}))
    (path.parent / "blob").write_text("tampered")
    with pytest.raises(ValueError, match="digest mismatch"):
        summarize_license_evidence.summarize(tmp_path, tmp_path)
