"""Template comparison/scoring controls, not extra diffusion qualification."""

import hashlib
import inspect
from pathlib import Path

import pytest
from pdbfixer import PDBFixer

from dvbfixer.model.diffusion.atom_completion_benchmark import (
    heldout_completion_metrics,
    template_completion_proposal,
)
from dvbfixer.model.diffusion.atom_completion_validation import audit_completion_candidate
from dvbfixer.model.diffusion.contract import DiffusionContractError
from tests.test_diffusion_atom_completion_sampling import sampling_control as sampling_control
from tests.test_diffusion_atom_completion_validation import completion_control as completion_control


def _reference(sampling_control):
    return b"".join(line for line in sampling_control[1].splitlines(keepends=True)
                    if not line.startswith(b"ATOM  ") or int(line[22:26]) <= 72)


@pytest.fixture
def rebuild_resources():
    if not Path(inspect.getfile(PDBFixer)).with_name("soft.xml").is_file():
        pytest.skip("installed PDBFixer soft.xml is missing; actual rebuild controls run separately in the intact CPU runtime")


def test_heldout_identity_metrics_use_fixed_frame_without_alignment(sampling_control):
    source, _, request = sampling_control
    reference = _reference(sampling_control)
    score = heldout_completion_metrics(request, source, reference, reference)
    assert score.atom_count == 4 and score.rmsd_angstrom == 0
    lines = reference.splitlines(keepends=True)
    for index, line in enumerate(lines):
        if line.startswith(b"ATOM  ") and line[22:26].strip() == b"27" and line[12:16].strip() == b"NZ":
            lines[index] = line[:30] + f"{float(line[30:38]) + 2:8.3f}".encode() + line[38:]
    changed = heldout_completion_metrics(request, source, reference, b"".join(lines))
    assert changed.maximum_displacement_angstrom == pytest.approx(2)
    assert changed.rmsd_angstrom == pytest.approx(1)  # One of four atoms moves 2 Å.


def test_reference_cannot_change_known_observations_or_omit_heldout_atoms(sampling_control):
    source, _, request = sampling_control
    reference = _reference(sampling_control)
    with pytest.raises(DiffusionContractError, match="missing or unrequested"):
        heldout_completion_metrics(request, source, source, reference)
    first = reference.splitlines(keepends=True)[0]
    changed = reference.replace(first, first[:30] + b"  80.000" + first[38:])
    with pytest.raises(DiffusionContractError, match="observation records"):
        heldout_completion_metrics(request, source, changed, reference)


def test_real_template_partial_proposal_preserves_known_atoms(sampling_control, rebuild_resources):
    source, _, request = sampling_control
    proposal = template_completion_proposal(request, source)
    audit_completion_candidate(source, proposal.candidate, request)
    assert proposal.temporary_fixed_max_displacement_angstrom < 1e-4
    assert all(atom.atom_name == "OXT" for atom in proposal.temporary_extra_atoms)
    assert "outer seed applies only to refinement" in proposal.seed_policy
    score = heldout_completion_metrics(request, source, _reference(sampling_control), proposal.candidate)
    assert score.atom_count == 4


def test_real_template_mixed_proposal_maps_gap_or_refuses_d(sampling_control, rebuild_resources):
    from dataclasses import replace

    from dvbfixer.ffutils.geometry import ChiralityError
    from dvbfixer.model.diffusion.atom_completion import plan_atom_completion
    from dvbfixer.model.diffusion.contract import GapRegion, ResidueIdentity, TargetInterval
    from dvbfixer.model.diffusion.protein_atoms import (
        _CANONICAL_ONE_TO_THREE,
        canonical_heavy_atom_identities,
    )

    source, _, request = sampling_control
    source = b"".join(line for line in source.splitlines(keepends=True)
                      if not (line.startswith(b"ATOM  ") and 10 <= int(line[22:26]) <= 13))
    placement = request.sequence_placements[0]
    pairs = [(index, residue) for index, residue in zip(placement.observed_target_indices, placement.observed_residues)
             if not 9 <= index < 13]
    placement = replace(placement, observed_target_indices=tuple(pair[0] for pair in pairs),
                        observed_residues=tuple(pair[1] for pair in pairs))
    generated = tuple(ResidueIdentity("A", str(index)) for index in range(10, 14))
    left, right = ResidueIdentity("A", "9"), ResidueIdentity("A", "14")
    gap = GapRegion("A", TargetInterval(9, 13), left, right, generated, (left, *generated, right))
    gap_atoms = tuple(atom for index, residue in zip(range(9, 13), generated)
                      for atom in canonical_heavy_atom_identities(residue, _CANONICAL_ONE_TO_THREE[request.target_sequences[0].sequence[index]]))
    digest = hashlib.sha256(source).hexdigest()
    masks = plan_atom_completion(source, source_sha256=digest, target_sequences=request.target_sequences,
                                 sequence_placements=(placement,), fixed_atoms=tuple(a for a in request.masks.fixed_atoms
                                 if not 10 <= int(a.residue_number) <= 13), gap_generated_atoms=gap_atoms)
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=digest),
                      sequence_placements=(placement,), gaps=(gap,), masks=masks)
    try:
        proposal = template_completion_proposal(request, source)
    except ChiralityError as exc:
        # A real, retained scientific baseline refusal, not permission to reflect
        # or to weaken the same raw-D admission policy used by completion.
        assert exc.residues
        assert all(10 <= int(residue[1]) <= 13 for residue in exc.residues)
        return
    audit_completion_candidate(source, proposal.candidate, request)
    score = heldout_completion_metrics(request, source, _reference(sampling_control), proposal.candidate)
    assert score.atom_count == 32
    assert proposal.temporary_fixed_max_displacement_angstrom < 1e-4


def test_baseline_does_not_bypass_original_occupancy_admission(completion_control):
    from dataclasses import replace

    from dvbfixer.model.diffusion.atom_completion_sampling import COMPLETION_RESEARCH_PROFILE

    source, _, request = completion_control
    with pytest.raises(DiffusionContractError, match="partial-occupancy"):
        template_completion_proposal(replace(request, profile=COMPLETION_RESEARCH_PROFILE), source)
