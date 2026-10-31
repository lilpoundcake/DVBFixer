"""Pre-diffusion template atoms become fixed context, never generated gap atoms."""

import hashlib
import json
from dataclasses import replace

import pytest

from dvbfixer.model.diffusion import input_preparation, observed_completion
from dvbfixer.model.diffusion.contract import ArtifactReference
from dvbfixer.model.diffusion.input_geometry import InputGeometryError
from tests.test_diffusion_input_preparation import _request


def incomplete(tmp_path, atoms=("NZ",)):
    path, complete, request = _request(tmp_path)
    source = b"".join(line for line in complete.splitlines(keepends=True)
                      if not (line.startswith(b"ATOM  ") and line[22:26].strip() == b"27"
                              and line[12:16].strip().decode() in atoms))
    path.write_bytes(source)
    removed = {a for a in request.fixed_atoms if a.residue_number == "27" and a.atom_name in atoms}
    request = replace(request, normalized_pdb=ArtifactReference(path.name, hashlib.sha256(source).hexdigest()),
                      fixed_atoms=tuple(a for a in request.fixed_atoms if a not in removed))
    return path, complete, source, request


def reference_rebuild(monkeypatch, complete):
    """Reference-derived unit stub, NOT checkpoint/template scientific evidence."""
    from dvbfixer.model.diffusion.atom_completion_materialize import _source_layout
    from dvbfixer.model.diffusion.input_preparation import _records
    from dvbfixer.model.diffusion.pdb_materialize import _generated_lines

    def rebuild(request, source, missing):
        records = _records(complete)
        represented = {a: (records[a][17:20].strip(), a.atom_name[0],
                           tuple(float(records[a][i:i + 8]) for i in (30, 38, 46))) for a in missing}
        _, serials, _ = _source_layout(source.decode())
        return _generated_lines(missing, represented, (s for s in range(1, 100_000) if s not in serials))

    monkeypatch.setattr(observed_completion, "_rebuild", rebuild)


def payloads(request, source, path):
    return {"input-original": source, "input-prepared": path.read_bytes(),
            "input-repair-report": (path.parent / "input-repair.json").read_bytes(),
            "input-repair-request": request.to_json().encode()}


@pytest.mark.parametrize("atoms", [("NZ",), ("CE", "NZ")])
def test_completion_precedes_context_and_preserves_gap_ownership(tmp_path, monkeypatch, atoms):
    path, complete, source, request = incomplete(tmp_path, atoms)
    reference_rebuild(monkeypatch, complete)
    prepared, output = input_preparation.prepare_input(request, source, tmp_path)
    assert path.read_bytes() == source
    assert set(prepared.fixed_atoms) - set(request.fixed_atoms) == set(observed_completion.missing_observed_atoms(request, source))
    assert prepared.generated_atoms == request.generated_atoms
    assert prepared.gaps == request.gaps
    assert prepared.template_ownership == request.template_ownership
    assert prepared.protein_context is None and prepared.gap_contexts is None
    assert not observed_completion.missing_observed_atoms(prepared, output.read_bytes())
    report = json.loads((output.parent / "input-repair.json").read_text())
    assert report["diffusion_generated"] is False
    assert report["whole_gaps_built"] is False
    assert len(report["added_atoms"]) == len(atoms)
    input_preparation.audit_preparation_payloads(prepared, payloads(prepared, source, output))
    assert input_preparation.preparation_options(prepared)["revision"] == observed_completion.REVISION


def test_observed_coordinate_change_cannot_be_hidden_by_updated_hashes(tmp_path, monkeypatch):
    _, complete, source, request = incomplete(tmp_path)
    reference_rebuild(monkeypatch, complete)
    prepared, output = input_preparation.prepare_input(request, source, tmp_path)
    lines = output.read_bytes().splitlines(keepends=True)
    index = next(i for i, line in enumerate(lines) if line.startswith(b"ATOM  "))
    lines[index] = lines[index][:30] + b"  99.000" + lines[index][38:]
    bad = b"".join(lines)
    report = json.loads((output.parent / "input-repair.json").read_text())
    report["prepared_sha256"] = hashlib.sha256(bad).hexdigest()
    prepared = replace(prepared, normalized_pdb=replace(prepared.normalized_pdb, sha256=report["prepared_sha256"]))
    with pytest.raises(InputGeometryError, match="changed observations"):
        input_preparation.audit_prepared_input(prepared, source, bad, report)


def test_partial_occupancy_is_refused_with_retained_original(tmp_path):
    _, _, source, request = incomplete(tmp_path)
    lines = source.splitlines(keepends=True)
    index = next(i for i, line in enumerate(lines) if line.startswith(b"ATOM  ") and line[22:26].strip() == b"27")
    lines[index] = lines[index][:54] + b"  0.50" + lines[index][60:]
    source = b"".join(lines)
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=hashlib.sha256(source).hexdigest()))
    with pytest.raises(InputGeometryError, match="partial occupancy"):
        input_preparation.prepare_input(request, source, tmp_path)
    directory = next(tmp_path.glob("dvbfixer-input-completion-*"))
    assert (directory / "source-original.pdb").read_bytes() == source
    assert json.loads((directory / "input-repair.json").read_text())["status"] == "failed"


def test_failed_rebuild_is_not_sampling_permission(tmp_path, monkeypatch):
    _, _, source, request = incomplete(tmp_path)

    def fail(*args):
        raise RuntimeError("test rebuild refusal")

    monkeypatch.setattr(observed_completion, "_rebuild", fail)
    with pytest.raises(InputGeometryError, match="Sampling was not started"):
        input_preparation.prepare_input(request, source, tmp_path)
    report = next(tmp_path.glob("dvbfixer-input-completion-*/input-repair.json"))
    assert json.loads(report.read_text())["status"] == "failed"


@pytest.mark.parametrize("atoms", [("NZ",), ("CE", "NZ")])
def test_real_missing_atoms_use_seeded_pdbfixer(tmp_path, atoms):
    _, _, source, request = incomplete(tmp_path, atoms)
    prepared, output = input_preparation.prepare_input(request, source, tmp_path)
    input_preparation.audit_preparation_payloads(prepared, payloads(prepared, source, output))
    assert len(prepared.fixed_atoms) == len(request.fixed_atoms) + len(atoms)


def test_request_builder_accepts_partial_sidechain_before_preparation(tmp_path, monkeypatch):
    from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

    path, complete, source, initial = incomplete(tmp_path)
    targets = {target.chain: target.sequence for target in initial.target_sequences}
    request = build_cli_diffusion_request(path, targets, seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    reference_rebuild(monkeypatch, complete)
    prepared, output = input_preparation.prepare_input(request, source, tmp_path)
    assert prepared.gaps == request.gaps
    assert prepared.generated_atoms == request.generated_atoms
    input_preparation.audit_preparation_payloads(prepared, payloads(prepared, source, output))


def test_preparation_evidence_survives_context_selection_and_freezing(tmp_path):
    from dvbfixer.model.diffusion.protein_context import prepare_cc94_request

    _, _, source, request = incomplete(tmp_path)
    prepared, output = input_preparation.prepare_input(request, source, tmp_path)
    # The original fixture has partial occupancy at 73–76; this is a source-bound
    # publication audit control, not a checkpoint sampling admission claim.
    prepared = prepare_cc94_request(prepared, output.read_bytes(), "auto")
    assert prepared.protein_context is not None
    frozen = type(prepared).from_json(prepared.to_json())
    input_preparation.audit_preparation_payloads(frozen, payloads(frozen, source, output))


def test_completion_geometry_failure_stops_with_evidence(tmp_path, monkeypatch):
    from dvbfixer.diagnose.report import Finding, Severity

    _, complete, source, request = incomplete(tmp_path)
    reference_rebuild(monkeypatch, complete)
    monkeypatch.setattr(observed_completion, "assess_input_geometry", lambda *_: (
        Finding(Severity.ERROR, "clash", "A", "27", "LYS", "test severe completion clash"),))
    with pytest.raises(InputGeometryError, match="severe geometry defects"):
        input_preparation.prepare_input(request, source, tmp_path)
    directory = next(tmp_path.glob("dvbfixer-input-completion-*"))
    assert (directory / "prepared.pdb").exists()
    assert json.loads((directory / "input-repair.json").read_text())["status"] == "failed"
