"""Declaration inventory stays deterministic and never invents legal approval."""

import importlib.util
import json
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("license_inventory", ROOT / "scripts/gen_license_inventory.py")
assert SPEC is not None and SPEC.loader is not None
inventory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inventory)


@pytest.fixture
def declared_repo(tmp_path):
    (tmp_path / "pyproject.toml").write_text('''
[build-system]
requires = ["setuptools>=64"]
[project]
version = "1.0"
dependencies = ["unknown-license-package>=1"]
[project.optional-dependencies]
dev = ["pytest"]
[tool.setuptools.data-files]
"share/FF" = ["FF/example.ff/*"]
[tool.setuptools.package-data]
sample = ["data/*.xml"]
''')
    catalog = tmp_path / inventory.CATALOG
    catalog.parent.mkdir(parents=True)
    catalog.write_text('''
[profiles.standard]
groups = ["science"]
[profiles.desktop]
inherits = "standard"
groups = ["gui"]
[addons.optional]
groups = ["optional"]
[[packages]]
spec = "example"
group = "science"
[[packages]]
spec = "extra"
group = "optional"
[external.operator]
kind = "operator"
''')
    (tmp_path / "gui").mkdir()
    (tmp_path / "gui/package.json").write_text(json.dumps({"dependencies": {"example": "^1"}}))
    (tmp_path / "gui/package-lock.json").write_text(json.dumps({"packages": {
        "": {"name": "root"}, "node_modules/example": {"version": "1.0", "license": "MIT"},
        "node_modules/unknown": {"version": "2.0"},
    }}))
    for path in ("FF/example.ff/forcefield.doc", "src/sample/data/caps.xml", "deploy/model/environment.yml"):
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("example\n")
    return tmp_path, {p.relative_to(tmp_path).as_posix() for p in tmp_path.rglob("*") if p.is_file()}


def test_inventory_is_deterministic_and_has_no_clearance(declared_repo):
    root, tracked = declared_repo
    result = inventory.build_inventory(root, tracked)
    assert inventory.render(result) == inventory.render(inventory.build_inventory(root, tracked))
    assert result["project_license_declaration"] is None
    assert result["gui_license_declaration"] is None
    assert result["kind"] == "declaration-inventory-not-sbom"
    assert all(r["compatibility_review"] == "unknown" for r in result["records"])
    assert {r["license_review"] for r in result["records"]} == {"unknown", "declared"}
    npm = [r for r in result["records"] if r["id"].startswith("npm-lock:")]
    assert len(npm) == 2
    assert all(r["distribution_inclusion"] == "unverified" for r in npm)
    assert result["record_count"] == len(result["records"])


def test_inherited_membership_and_addons_are_distinct(declared_repo):
    root, tracked = declared_repo
    records = inventory.build_inventory(root, tracked)["records"]
    assert next(r for r in records if r["id"] == "catalog:0")["profiles"] == ["desktop", "standard"]
    addon = next(r for r in records if r["id"] == "catalog:1")
    assert addon["profiles"] == []
    assert addon["addons"] == ["optional"]


def test_drift_is_visible_and_local_untracked_data_is_excluded(declared_repo):
    root, tracked = declared_repo
    before = inventory.build_inventory(root, tracked)
    (root / "FF/example.ff/local.itp").write_text("user work")
    assert inventory.build_inventory(root, tracked) == before
    (root / "FF/example.ff/forcefield.doc").write_text("changed")
    assert inventory.build_inventory(root, tracked) != before
    (root / "deploy/model/environment.yml").write_text("changed runtime")
    after = inventory.build_inventory(root, tracked)
    assert after["input_sha256"] != before["input_sha256"]
    pending = next(r for r in after["records"] if r["source"] == "deploy/model/environment.yml")
    assert pending["dependency_enumeration"] == "pending"


def test_missing_declared_data_fails_closed(declared_repo):
    root, tracked = declared_repo
    tracked.remove("FF/example.ff/forcefield.doc")
    with pytest.raises(ValueError, match="No tracked files"):
        inventory.build_inventory(root, tracked)


def test_project_legal_text_changes_are_bound_even_before_git_add(declared_repo):
    root, tracked = declared_repo
    before = inventory.build_inventory(root, tracked)
    (root / "LICENSE").write_text("approved project terms\n")
    after = inventory.build_inventory(root, tracked)
    assert "LICENSE" not in before["project_legal_files_sha256"]
    assert after["project_legal_files_sha256"]["LICENSE"] == after["input_sha256"]["LICENSE"]
    (root / "LICENSE").write_text("changed terms\n")
    assert inventory.build_inventory(root, tracked) != after


def test_inheritance_cycles_fail(declared_repo):
    root, tracked = declared_repo
    path = root / inventory.CATALOG
    path.write_text(path.read_text().replace('[profiles.standard]', '[profiles.standard]\ninherits = "desktop"'))
    with pytest.raises(ValueError, match="cycle"):
        inventory.build_inventory(root, tracked)


def test_check_mode_reports_drift_without_rewriting(declared_repo, monkeypatch, capsys):
    root, tracked = declared_repo
    output = root / inventory.OUTPUT
    output.parent.mkdir(parents=True)
    output.write_text(inventory.render(inventory.build_inventory(root, tracked)))
    monkeypatch.setattr(inventory, "ROOT", root)
    monkeypatch.setattr(inventory, "tracked_files", lambda _: tracked)
    monkeypatch.setattr(sys, "argv", ["gen_license_inventory.py", "--check"])
    assert inventory.main() == 0
    before = output.read_bytes()
    (root / "gui/package.json").write_text('{"dependencies": {"new": "*"}}')
    assert inventory.main() == 1
    assert "Out of date" in capsys.readouterr().out
    assert output.read_bytes() == before


def test_repository_inventory_is_current():
    subprocess.run([sys.executable, str(ROOT / "scripts/gen_license_inventory.py"), "--check"],
                   check=True, cwd=ROOT)


def test_generated_inventory_covers_python_requirements_and_preserves_unknowns():
    data = json.loads((ROOT / inventory.OUTPUT).read_text())
    runtime = [r for r in data["records"] if r["boundary"] == "python:runtime"]
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert len(runtime) == 20
    assert sorted(r["requirement"] for r in runtime) == sorted(project["project"]["dependencies"])
    assert all(r["compatibility_review"] == "unknown" for r in data["records"])
    assert len({r["id"] for r in data["records"]}) == len(data["records"])


def test_original_code_license_metadata_and_gui_texts_are_synchronized():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert project["project"]["license"] == {"text": "BSD-3-Clause"}
    assert project["tool"]["setuptools"]["license-files"] == ["LICENSE"]
    gui = json.loads((ROOT / "gui/package.json").read_text())
    lock = json.loads((ROOT / "gui/package-lock.json").read_text())
    assert gui["license"] == lock["packages"][""]["license"] == "BSD-3-Clause"
    text = (ROOT / "LICENSE").read_bytes()
    assert b"Copyright (c) 2026, Stanislav Skripkin" in text
    assert b"<year>" not in text and b"<copyright holder>" not in text
    assert (ROOT / "gui/LICENSE").read_bytes() == text
    assert (ROOT / "gui/public/LICENSE").read_bytes() == text
    assert "include LICENSE CONTRIBUTING.md" in (ROOT / "MANIFEST.in").read_text()


def test_licensing_completion_is_scoped_and_does_not_claim_third_party_clearance():
    record = json.loads((ROOT / "docs/research/licensing-completion-2026-11-02.json").read_text())
    assert record["task_status"] == "complete-for-agreed-scope"
    assert record["owner_action_required_for_this_task"] is False
    assert record["original_code"]["license"] == "BSD-3-Clause"
    assert record["delivery"]["primary"] == "source-plus-existing-standard-desktop-installer"
    assert record["delivery"]["prebuilt_environment_container_publication_in_scope"] is False
    assert record["delivery"]["restricted_binary_or_weight_redistribution_authorized"] is False
    assert record["overall_third_party_clearance"] is False
    assert "resolved-gpl-combined-work-interpretation" in record["not_claimed"]


def test_bundled_force_fields_keep_data_and_notice_packaging():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    data = project["tool"]["setuptools"]["data-files"]
    for name in ("amber19sb.ff", "amber99sb-ildn-lipid21.ff", "charmm36_ljpme-jul2022.ff"):
        assert data[f"share/dvbfixer/FF/{name}"] == [f"FF/{name}/*"]
        assert (ROOT / "FF" / name / "forcefield.doc").is_file()
    assert data["share/dvbfixer/notices"] == ["THIRD_PARTY_NOTICES/*"]
    notice = (ROOT / "THIRD_PARTY_NOTICES/FORCE-FIELDS.md").read_text()
    assert "charmm-terminal-caps.xml" in notice
    assert "1bbd3bf" in notice
    assert "not relicense" in notice
