"""Audit controls use reference coordinates, not checkpoint sampling evidence."""
from pathlib import Path

import pytest
from scripts.audit_mac_reconstruction_smoke import audit, inventory
from scripts.prepare_mac_reconstruction_smoke import prepare

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


def test_audit_reference_control_binds_mask_and_retains_observations(tmp_path):
    inputs = tmp_path / "input"
    prepare(FIXTURE, inputs)
    output = tmp_path / "audit"
    report = audit(inputs, inputs / "reference-projection-NOT-CONDITIONING.pdb", output)
    assert report["exact_sequence_and_atom_inventory"]
    assert report["heldout_generated_rmsd_angstrom"] < 1e-10
    assert report["assessments"]["source-restored"]["geometry_passed"]
    assert (output / "source.pdb").read_bytes() == (inputs / "source.pdb").read_bytes()
    original = [line for line in (inputs / "source.pdb").read_text().splitlines() if line.startswith("ATOM  ")]
    restored = (output / "source-restored.pdb").read_text().splitlines()
    assert all(line in restored for line in original)
    assert not report["scientific_qualification"]


def test_audit_refuses_missing_atom_inventory(tmp_path):
    inputs = tmp_path / "input"
    prepare(FIXTURE, inputs)
    with pytest.raises(ValueError, match="inventory"):
        audit(inputs, inputs / "source.pdb", tmp_path / "audit")
    assert not (tmp_path / "audit").exists()


def test_inventory_refuses_duplicate_atoms(tmp_path):
    line = next(line for line in FIXTURE.read_text().splitlines(True) if line.startswith("ATOM  "))
    path = tmp_path / "duplicate.pdb"
    path.write_text(line + line + "END\n")
    with pytest.raises(ValueError, match="duplicate"):
        inventory(path)


def test_audit_refuses_changed_declared_reference_bytes(tmp_path):
    inputs = tmp_path / "input"
    prepare(FIXTURE, inputs)
    reference = inputs / "reference-projection-NOT-CONDITIONING.pdb"
    reference.write_bytes(reference.read_bytes()+b"REMARK changed after declaration\n")
    with pytest.raises(ValueError, match="projection bytes changed"):
        audit(inputs, reference, tmp_path / "audit")
