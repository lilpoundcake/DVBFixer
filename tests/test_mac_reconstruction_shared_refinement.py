"""Mask/orchestration controls with a fake refiner, not scientific evidence."""
from pathlib import Path
from types import SimpleNamespace

import pytest
from scripts import refine_mac_reconstruction_smoke as module
from scripts.prepare_mac_reconstruction_smoke import prepare

from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment, RepairIssue
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


def case(tmp_path):
    inputs = tmp_path / "inputs"
    manifest = prepare(FIXTURE, inputs)
    source = "".join(line for line in (inputs / "source.pdb").read_text().splitlines(True)
                     if not (line.startswith("ATOM  ") and 20 <= int(line[22:26]) <= 21))
    work = tmp_path / "audit"
    work.mkdir()
    source_path = work / "source.pdb"
    source_path.write_text(source)
    request = build_cli_diffusion_request(source_path, {"A": manifest["sequence"]}, seeds=(7,),
                                         profile="protpardelle-1c-cc94-mps-v1")
    assert len(request.gaps) == 2
    (work / "evaluation-mask.json").write_text(request.to_json())
    candidate = (inputs / "reference-projection-NOT-CONDITIONING.pdb").read_text()
    (work / "source-restored.pdb").write_text(candidate)
    return work, request


def test_shared_refiner_gets_all_gaps_not_only_the_flagged_gap(tmp_path, monkeypatch):
    work, request = case(tmp_path)
    results = iter((RepairAssessment((RepairIssue("clash", request.gaps[0].generated_residues[:1], "test-only finding"),)),
                    RepairAssessment()))
    monkeypatch.setattr(module, "assessment_adapter", lambda *_: lambda *_: next(results))
    captured = {}
    def fake_refiner(text, **kwargs):
        captured.update(kwargs)
        coordinates = _coordinates(text)[0]
        return SimpleNamespace(coordinates_angstrom={a: coordinates[a] for a in kwargs["generated_atoms"]},
            initial_energy_kj_mol=0, final_energy_kj_mol=0, platform="Reference", search_diagnostics=())
    monkeypatch.setattr(module, "refine_generated_region", fake_refiner)
    module.run(work, tmp_path / "refinement")
    assert set(captured["generated_atoms"]) == set(request.generated_atoms)
    assert set(captured["generated_residues"]) == {r for gap in request.gaps for r in gap.generated_residues}
    assert not set(captured["generated_atoms"]) & set(request.fixed_atoms)
    assert "observed_anchor_restraints" not in captured


def test_clean_assembly_still_runs_mandatory_shared_refinement(tmp_path, monkeypatch):
    work, request = case(tmp_path)
    monkeypatch.setattr(module, "assessment_adapter", lambda *_: lambda *_: RepairAssessment())
    calls = []
    def fake_refiner(text, **kwargs):
        calls.append(kwargs)
        coordinates = _coordinates(text)[0]
        return SimpleNamespace(coordinates_angstrom={a:coordinates[a] for a in kwargs["generated_atoms"]},
            initial_energy_kj_mol=0, final_energy_kj_mol=0, platform="Reference", search_diagnostics=())
    monkeypatch.setattr(module, "refine_generated_region", fake_refiner)
    module.run(work, tmp_path / "refinement", force_shared_refinement=True)
    assert len(calls) == 1
    assert set(calls[0]["generated_atoms"]) == set(request.generated_atoms)


def test_source_admission_failure_never_authorizes_refinement(tmp_path, monkeypatch):
    work, _ = case(tmp_path)
    monkeypatch.setattr(module, "assessment_adapter", lambda *_: lambda *_: RepairAssessment((), ("source-motion",)))
    def never_refine(*_args, **_kwargs):
        pytest.fail("source admission failure reached coordinate repair")
    monkeypatch.setattr(module, "refine_generated_region", never_refine)
    with pytest.raises(ValueError, match="admission"):
        module.run(work, tmp_path / "refinement")
