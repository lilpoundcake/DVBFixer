"""Pre-sampling repair must stay local, bounded, auditable and source-preserving."""

import hashlib
import json
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from Bio import SeqIO

from dvbfixer.diagnose.report import Finding, Severity
from dvbfixer.model.diffusion import input_preparation as preparation
from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.input_geometry import InputGeometryError, assess_input_geometry
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

FIXTURES = Path(__file__).parent / "fixtures/protein_context"


def _request(tmp_path, name="1UBQ", chain="A", start=10, stop=13):
    fixture = FIXTURES / f"{name}.pdb"
    sequences = {r.annotations["chain"]: str(r.seq) for r in SeqIO.parse(fixture, "pdb-seqres")}
    source = "".join(line for line in fixture.read_text().splitlines(keepends=True)
                     if not line.startswith("HETATM") and not
                     (line.startswith("ATOM  ") and line[21] == chain and start <= int(line[22:26]) <= stop)).encode()
    path = tmp_path / "input.pdb"
    path.write_bytes(source)
    return path, source, build_cli_diffusion_request(path, sequences, seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")


def _mock_damage(monkeypatch, source, displacement=0.1):
    # Test the controller separately from the geometry detector / OpenMM kernel.
    atoms = [{"chain": "A", "resid": "20", "icode": "", "atom": "OG"},
             {"chain": "A", "resid": "21", "icode": "", "atom": "O"}]
    finding = Finding(Severity.ERROR, "clash", "A", "20", "SER", "test measured clash", extra={"atoms": atoms})
    monkeypatch.setattr(preparation, "assess_input_geometry", lambda request, data: (finding,) if data == source else ())
    calls = []

    def minimize(text, *, residues, atoms, restraints, max_iterations):
        calls.append((text, residues, atoms, restraints, max_iterations))
        coordinates = {r.atom: np.array(r.reference_angstrom) for r in restraints}
        coordinates[AtomIdentity("A", "20", "", "OG")][0] += displacement
        return SimpleNamespace(coordinates_angstrom=coordinates, initial_energy_kj_mol=10.0,
                               final_energy_kj_mol=1.0, chirality_repairs=())

    monkeypatch.setattr("dvbfixer.model.diffusion.boundary_refinement.refine_observed_region", minimize)
    return calls


def test_repair_rebases_source_without_reclassifying_observed_atoms(tmp_path, monkeypatch, capsys):
    path, source, request = _request(tmp_path)
    calls = _mock_damage(monkeypatch, source)
    prepared_request, prepared_path = preparation.prepare_input(request, source, tmp_path)
    assert len(calls) == 1
    assert prepared_request.fixed_atoms == request.fixed_atoms
    assert prepared_request.generated_atoms == request.generated_atoms
    assert prepared_request.gaps == request.gaps
    assert prepared_request.template_ownership == request.template_ownership
    assert prepared_request.protein_context is None
    assert prepared_request.normalized_pdb.sha256 == hashlib.sha256(prepared_path.read_bytes()).hexdigest()
    assert path.read_bytes() == source
    report_bytes = (prepared_path.parent / "input-repair.json").read_bytes()
    report = json.loads(report_bytes)
    assert report["status"] == "accepted"
    assert all(m["within_limits"] for m in report["movements"])
    payloads = {"input-original": source, "input-prepared": prepared_path.read_bytes(),
                "input-repair-report": report_bytes, "input-repair-request": prepared_request.to_json().encode()}
    preparation.audit_preparation_payloads(prepared_request, payloads)
    assert (prepared_path.parent / "source-original.pdb").read_bytes() == source
    assert "Input repair accepted" in capsys.readouterr().out
    # The missing region is not a spring between two otherwise distant anchors.
    assert "TER\n" in calls[0][0]
    assert report["policy"]["solvent"] is False
    assert report["policy"]["restart_count"] == 0


def test_large_motion_is_rejected_even_when_geometry_assessor_passes(tmp_path, monkeypatch):
    path, source, request = _request(tmp_path)
    _mock_damage(monkeypatch, source, displacement=2.0)
    with pytest.raises(InputGeometryError, match="displacement bounds"):
        preparation.prepare_input(request, source, tmp_path)
    report = json.loads(next(tmp_path.glob("dvbfixer-input-repair-*/input-repair.json")).read_text())
    assert report["status"] == "failed"
    assert path.read_bytes() == source
    assert next(tmp_path.glob("dvbfixer-input-repair-*/prepared.pdb")).exists()


def test_unrelated_coordinate_or_metadata_changes_cannot_pass_audit(tmp_path, monkeypatch):
    _path, source, request = _request(tmp_path)
    _mock_damage(monkeypatch, source)
    request, path = preparation.prepare_input(request, source, tmp_path)
    report = json.loads((path.parent / "input-repair.json").read_text())
    # Even an attacker updating all outer SHA bindings cannot authorize this.
    bad = path.read_bytes().replace(b"HEADER", b"HEADEQ", 1)
    report["prepared_sha256"] = hashlib.sha256(bad).hexdigest()
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=report["prepared_sha256"]))
    with pytest.raises(InputGeometryError, match="outside the authorized"):
        preparation.audit_prepared_input(request, source, bad, report)


def test_no_minimization_or_artifacts_for_healthy_input(tmp_path, monkeypatch):
    path, source, request = _request(tmp_path)
    monkeypatch.setattr("dvbfixer.model.diffusion.boundary_refinement.refine_observed_region",
                        lambda *_args, **_kwargs: pytest.fail("healthy input must not be minimized"))
    assert preparation.prepare_input(request, source, tmp_path) is None
    assert path.read_bytes() == source
    assert not list(tmp_path.glob("dvbfixer-input-repair-*"))


def test_ineligible_damage_keeps_source_and_failure_report(tmp_path, monkeypatch):
    path, source, request = _request(tmp_path)
    finding = Finding(Severity.ERROR, "coincident_atoms", "A", "20", "SER", "ambiguous coincident input")
    monkeypatch.setattr(preparation, "assess_input_geometry", lambda *_args: (finding,))
    with pytest.raises(InputGeometryError, match="cannot safely handle coincident_atoms"):
        preparation.prepare_input(request, source, tmp_path)
    report = json.loads(next(tmp_path.glob("dvbfixer-input-repair-*/input-repair.json")).read_text())
    assert report["status"] == "failed"
    assert path.read_bytes() == source


@pytest.mark.slow
def test_real_1tim_repairs_all_severe_input_findings_with_bounded_local_motion(tmp_path):
    path, source, request = _request(tmp_path, name="1TIM", start=73, stop=78)
    assert sum(f.severity is Severity.ERROR for f in assess_input_geometry(request, source)) == 29
    prepared_request, prepared_path = preparation.prepare_input(request, source, tmp_path)
    report = json.loads((prepared_path.parent / "input-repair.json").read_text())
    assert report["status"] == "accepted"
    assert len(report["selected_residues"]) == 47
    movements, findings = preparation.audit_prepared_input(prepared_request, source, prepared_path.read_bytes(), report)
    assert not [f for f in findings if f.severity is Severity.ERROR]
    assert all(m.within_limits for m in movements)
    assert path.read_bytes() == source
    assert prepared_request.fixed_atoms == request.fixed_atoms
    assert prepared_request.generated_atoms == request.generated_atoms
    assert [asdict(m) for m in movements] == report["movements"]
