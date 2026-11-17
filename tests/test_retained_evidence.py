import json
from pathlib import Path

import pytest
from scripts.mac_reconstruction_comparison import digest
from scripts.retained_evidence import freeze_code, verify_bindings


def declaration(tmp_path):
    code, source = tmp_path / "runner.py", tmp_path / "source.pdb"
    code.write_text("original = True\n")
    source.write_text("original data\n")
    bindings = {str(p): digest(p) for p in (code, source)}
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"source_code_and_input_sha256": bindings}))
    freeze_code(plan)
    return code, source, bindings, plan


def test_archived_python_must_match_the_exact_declared_bytes(tmp_path):
    code, source, bindings, plan = declaration(tmp_path)
    code.write_text("later = True\n")
    resolved = verify_bindings(bindings, (plan,))
    assert resolved[str(code)] != code
    assert digest(resolved[str(code)]) == bindings[str(code)]
    assert resolved[str(source)] == source


def test_archive_never_excuses_changed_structural_input(tmp_path):
    _, source, bindings, plan = declaration(tmp_path)
    source.write_text("changed data\n")
    with pytest.raises(ValueError, match="changed or disappeared"):
        verify_bindings(bindings, (plan,))


def test_archive_refuses_tampered_snapshot_or_changed_plan(tmp_path):
    code, _, bindings, plan = declaration(tmp_path)
    code.write_text("changed = True\n")
    index = json.loads((tmp_path / "executed-code/index.json").read_text())
    Path(index["snapshots"][0]["snapshot_path"]).write_text("not declared\n")
    with pytest.raises(ValueError, match="changed or disappeared"):
        verify_bindings(bindings, (plan,))
    plan.write_text(plan.read_text() + "\n")
    with pytest.raises(ValueError, match="different declaration"):
        verify_bindings(bindings, (plan,))
