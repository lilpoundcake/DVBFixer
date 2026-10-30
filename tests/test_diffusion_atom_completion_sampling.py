"""Private completion transport/admission does not enable production profiles."""

import hashlib
import importlib.util
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from dvbfixer.model.diffusion.atom_completion import (
    plan_atom_completion,
    plan_request_atom_completion,
)
from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
from dvbfixer.model.diffusion.atom_completion_sampling import (
    COMPLETION_RESEARCH_PROFILE,
    assess_completion_sampling,
    completion_engine_view,
)
from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    GapRegion,
    ResidueIdentity,
    TargetInterval,
)
from dvbfixer.model.diffusion.protein_atoms import canonical_heavy_atom_identities
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request
from tests.test_diffusion_atom_completion_validation import completion_control as completion_control


def _private(request):
    return replace(request, profile=COMPLETION_RESEARCH_PROFILE)


@pytest.fixture
def sampling_control(completion_control):
    """Explicit 1-72 subset, not an occupancy-normalized original 1UBQ."""
    source, reference, request = completion_control
    source = b"".join(line for line in source.splitlines(keepends=True)
                      if not line.startswith(b"ATOM  ") or int(line[22:26]) <= 72)
    digest = hashlib.sha256(source).hexdigest()
    request = replace(request, profile=COMPLETION_RESEARCH_PROFILE,
                      normalized_pdb=replace(request.normalized_pdb, sha256=digest),
                      target_sequences=(replace(request.target_sequences[0], sequence=request.target_sequences[0].sequence[:72]),),
                      sequence_placements=(replace(request.sequence_placements[0], target_length=72,
                                                   observed_target_indices=tuple(range(72)),
                                                   observed_residues=request.sequence_placements[0].observed_residues[:72]),),
                      masks=replace(request.masks, source_sha256=digest,
                                    fixed_atoms=tuple(a for a in request.masks.fixed_atoms if int(a.residue_number) <= 72)))
    return source, reference, request


def _rebind(request, source):
    digest = hashlib.sha256(source).hexdigest()
    return replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=digest),
                   masks=replace(request.masks, source_sha256=digest))


def _mixed(completion_control, tmp_path):
    source, _, partial = completion_control
    source = b"".join(line for line in source.splitlines(keepends=True)
                      if not (line.startswith(b"ATOM  ") and 10 <= int(line[22:26]) <= 13))
    path = tmp_path / "mixed.pdb"
    path.write_bytes(source)
    gaps = build_cli_diffusion_request(path, {"A": partial.target_sequences[0].sequence},
                                       seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    request = AtomCompletionRequest(1, gaps.normalized_pdb, gaps.target_sequences, gaps.sequence_placements,
                                    gaps.gaps, plan_request_atom_completion(gaps, source), (7,), COMPLETION_RESEARCH_PROFILE)
    return source, request


def test_partial_transport_preserves_masks_and_native_positions(sampling_control):
    source, _, request = sampling_control
    request = _private(request)
    admission = assess_completion_sampling(request, source)
    assert admission.supported, admission.reasons
    view = completion_engine_view(request, source)
    assert view.motif_placement == "A1-72"
    assert view.fixed_positions == tuple(range(72))
    assert view.residue_index == tuple(range(1, 73))
    assert view.chain_index == (0,) * 72
    assert set(request.generated_atoms).isdisjoint(view.represented_fixed_atoms)
    assert set(view.atom_axis) == set(view.represented_fixed_atoms) | set(request.generated_atoms)
    lys = [line for line in view.motif_text.splitlines() if line.startswith("ATOM  ") and line[22:26].strip() == "27"]
    assert {line[12:16].strip() for line in lys} == {"N", "CA", "CB", "C", "O"}
    assert source == sampling_control[0]


def test_mixed_transport_declares_real_gap_and_diagnostic_bridge_exclusion(sampling_control, tmp_path):
    source, request = _mixed(sampling_control, tmp_path)
    admission = assess_completion_sampling(request, source)
    assert admission.supported, admission.reasons
    assert admission.excluded_gap_bridges == ((AtomIdentity("A", "9", "", "C"), AtomIdentity("A", "14", "", "N")),)
    assert admission.source_findings  # Original parser-derived findings remain retained.
    view = completion_engine_view(request, source)
    assert view.motif_placement == "A1-9/4/A14-72"
    assert view.fixed_positions == tuple(range(9)) + tuple(range(13, 72))
    assert len(view.target_residues) == 72
    assert len(request.generated_atoms) == 32


def test_source_case_and_insertion_codes_round_trip_without_engine_truncation(sampling_control):
    source, _, request = sampling_control
    source = b"".join(line[:21] + b"d" + line[22:26] + (b"A" if line[22:26].strip() == b"27" else line[26:27]) + line[27:]
                      if line.startswith(b"ATOM  ") else line for line in source.splitlines(keepends=True))

    def residue(item):
        return ResidueIdentity("d", item.residue_number, "A" if item.residue_number == "27" else item.insertion_code)

    def atom(item):
        return AtomIdentity("d", item.residue_number, "A" if item.residue_number == "27" else item.insertion_code, item.atom_name)

    masks = replace(request.masks, fixed_atoms=tuple(atom(a) for a in request.masks.fixed_atoms),
                    partial_residues=tuple(replace(p, residue=residue(p.residue), observed_atoms=tuple(atom(a) for a in p.observed_atoms),
                                                   missing_atoms=tuple(atom(a) for a in p.missing_atoms)) for p in request.masks.partial_residues))
    request = replace(request, profile=COMPLETION_RESEARCH_PROFILE, target_sequences=(replace(request.target_sequences[0], chain="d"),),
                      sequence_placements=(replace(request.sequence_placements[0], chain="d",
                                                   observed_residues=tuple(residue(r) for r in request.sequence_placements[0].observed_residues)),), masks=masks)
    request = _rebind(request, source)
    view = completion_engine_view(request, source)
    assert view.target_residues[26] == ResidueIdentity("d", "27", "A")
    assert all(atom.chain == "d" for atom in view.atom_axis)
    assert any(atom == AtomIdentity("d", "27", "A", "NZ") for atom in view.atom_axis)
    assert all(line[21] == "A" and line[26] == " " for line in view.motif_text.splitlines() if line.startswith("ATOM  "))
    assert assess_completion_sampling(request, source).supported


def test_profile_and_extra_seed_are_denied_before_sampling(completion_control):
    source, _, request = completion_control
    assert "not-the-private-completion-profile" in assess_completion_sampling(request, source).reasons
    assert "research-invocation-requires-one-predeclared-seed" in assess_completion_sampling(
        replace(_private(request), seeds=(7, 19)), source).reasons


def test_bad_source_geometry_and_occupancy_are_not_repaired(sampling_control):
    source, _, request = sampling_control
    request = _private(request)
    first = source.splitlines(keepends=True)[0]
    damaged = source.replace(first, first[:30] + b"  80.000  80.000  80.000" + first[54:])
    report = assess_completion_sampling(_rebind(request, damaged), damaged)
    assert not report.supported
    assert "source-geometry-errors" in report.reasons
    assert report.source_findings
    occupancy = source.replace(first, first[:54] + b"  0.50" + first[60:])
    assert "partial-occupancy-not-enabled" in assess_completion_sampling(_rebind(request, occupancy), occupancy).reasons


def test_original_ubq_partial_occupancies_remain_refused(completion_control):
    source, _, request = completion_control
    assert "partial-occupancy-not-enabled" in assess_completion_sampling(_private(request), source).reasons


def _runner():
    path = Path(__file__).resolve().parents[1] / "deploy/protpardelle-1c/atom_completion_research.py"
    spec = importlib.util.spec_from_file_location("test_completion_research", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_research_runner_requires_explicit_opt_in_before_output(tmp_path):
    module = _runner()
    output = tmp_path / "not-created"
    with pytest.raises(ValueError, match="acknowledge"):
        module.run(SimpleNamespace(acknowledge_experimental_completion=False, output=output))
    assert not output.exists()


def test_refused_source_never_reaches_preflight_or_weights(completion_control, tmp_path, monkeypatch):
    source, _, request = completion_control
    first = source.splitlines(keepends=True)[0]
    source = source.replace(first, first[:54] + b"  0.50" + first[60:])
    request = _rebind(_private(request), source)
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, path="source.pdb"))
    (tmp_path / "source.pdb").write_bytes(source)
    path = tmp_path / "request.json"
    path.write_text(request.to_json())
    output = tmp_path / "refused"
    module = _runner()
    monkeypatch.setattr(module, "_load", lambda *args: pytest.fail("refused source reached engine preflight"))
    args = SimpleNamespace(acknowledge_experimental_completion=True, request=path, output=output)
    with pytest.raises(ValueError, match="partial-occupancy"):
        module.run(args)
    trace = json.loads((output / "trace.json").read_text())
    assert trace["checkpoint_sampling_executed"] is False
    assert trace["status"] == "failed"
    assert (output / "source.pdb").read_bytes() == source
    with pytest.raises(FileExistsError):
        module.run(args)  # Crash-left or failed reservations are never adopted.


def test_cc94_conditioning_is_native_not_cc89_intervention():
    config = _runner().conditioning_config((0, 1, 3))
    assert config["discontiguous_motif_assignment"]["fixed_motif_pos"] == [0, 1, 3]
    assert config["crop_conditional_guidance"]["strategy"] == "backbone-sidechain"
    assert not config["replacement_guidance"]["enabled"]
    assert not config["reconstruction_guidance"]["enabled"]


def test_known_peptide_cannot_be_reclassified_as_a_gap(sampling_control):
    source, _, request = sampling_control
    target = request.target_sequences[0]
    target = replace(target, sequence=target.sequence[:9] + "AAAA" + target.sequence[9:])
    placement = request.sequence_placements[0]
    placement = replace(placement, target_length=76,
                        observed_target_indices=tuple(index if index < 9 else index + 4 for index in placement.observed_target_indices))
    generated = tuple(ResidueIdentity("A", "9", code) for code in "ABCD")
    left, right = placement.observed_residues[8:10]
    gap = GapRegion("A", TargetInterval(9, 13), left, right, generated, (left, *generated, right))
    atoms = tuple(atom for residue in generated for atom in canonical_heavy_atom_identities(residue, "ALA"))
    masks = plan_atom_completion(source, source_sha256=request.normalized_pdb.sha256,
                                 target_sequences=(target,), sequence_placements=(placement,),
                                 fixed_atoms=request.masks.fixed_atoms, gap_generated_atoms=atoms)
    request = replace(request, target_sequences=(target,), sequence_placements=(placement,), gaps=(gap,), masks=masks)
    admission = assess_completion_sampling(request, source)
    assert not admission.supported
    assert "gap-crosses-observed-peptide-connectivity" in admission.reasons
    assert not admission.excluded_gap_bridges


def test_partial_backbone_is_explicitly_refused(sampling_control):
    source, _, request = sampling_control
    source = b"".join(line for line in source.splitlines(keepends=True)
                      if not (line.startswith(b"ATOM  ") and line[22:26].strip() == b"27" and line[12:16].strip() == b"N"))
    digest = hashlib.sha256(source).hexdigest()
    masks = plan_atom_completion(source, source_sha256=digest, target_sequences=request.target_sequences,
                                 sequence_placements=request.sequence_placements,
                                 fixed_atoms=tuple(a for a in request.masks.fixed_atoms if a != AtomIdentity("A", "27", "", "N")))
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=digest), masks=masks)
    assert "partial-backbone-not-enabled" in assess_completion_sampling(request, source).reasons
