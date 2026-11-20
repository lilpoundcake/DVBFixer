"""FrameDiPT transport controls do not import Torch, model code or checkpoints."""
from types import SimpleNamespace

import numpy as np
import pytest
from deploy.framedipt.research_worker import project_fixed_frames
from scripts.audit_framedipt_research import compare_frames
from scripts.framedipt_research import encode, validate_transport

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def case():
    targets, placements, gaps, fixed, generated, source = [], [], [], [], [], {}
    for chain in ("D", "d"):
        residues = [ResidueIdentity(chain, str(i+80), "A" if i == 2 else "") for i in range(6)]
        targets.append(SimpleNamespace(chain=chain, sequence="AAAAAA"))
        placements.append(SimpleNamespace(chain=chain, observed_target_indices=[0, 4, 5],
                                           observed_residues=[residues[i] for i in (0, 4, 5)]))
        gaps.append(SimpleNamespace(chain=chain, left_anchor=residues[0], right_anchor=residues[4],
                    generated_residues=residues[1:4], target_interval=SimpleNamespace(start=1, stop=4)))
        for i, r in enumerate(residues):
            for j, name in enumerate(CANONICAL_HEAVY_ATOMS["ALA"]):
                a = AtomIdentity(chain, r.residue_number, r.insertion_code, name)
                if 1 <= i < 4:
                    generated.append(a)
                else:
                    fixed.append(a)
                    source[chain, r.residue_number, r.insertion_code, name] = {
                        "line": "ATOM  observation", "xyz": [i+1., j+1., 1. if chain == "D" else -1.],
                        "element": name[0], "resname": "ALA"}
    request = SimpleNamespace(fixed_atoms=fixed, generated_atoms=generated, target_sequences=targets,
        sequence_placements=placements, gaps=gaps, retained_explicit_links=(), heterogen_contexts=())
    return request, source


def test_complete_case_sensitive_axis_and_gap_masks():
    request, source = case()
    arrays, receipt = encode(request, source)
    assert receipt["target_residues"] == 12
    assert receipt["fixed_residues"] == receipt["generated_residues"] == 6
    assert receipt["chains"] == [{"chain": "D", "length": 6}, {"chain": "d", "length": 6}]
    assert receipt["residues"][2]["insertion_code"] == "A"
    assert arrays["seq_idx"].tolist() == list(range(6)) + list(range(206, 212))
    assert arrays["chain_idx"].tolist() == [0]*6 + [1]*6
    assert np.all(arrays["res_mask"] == 1)
    assert receipt["sidechain_conditioning"] is False


def test_absent_positions_never_contain_reference_coordinates():
    request, source = case()
    arrays, _ = encode(request, source)
    missing = arrays["fixed_mask"] == 0
    assert not arrays["atom_positions"][missing].any()
    assert not arrays["atom_mask"][missing].any()
    arrays["atom_positions"][1, 0] = [1, 2, 3]
    with pytest.raises(ValueError, match="unknown region"):
        validate_transport(arrays)


def test_gap_residues_cannot_be_dropped_by_backbone_mask():
    request, source = case()
    arrays, _ = encode(request, source)
    arrays["res_mask"][1] = 0
    with pytest.raises(ValueError, match="remain real"):
        validate_transport(arrays)


def test_partial_source_and_partial_generated_ownership_refused():
    request, source = case()
    a = request.generated_atoms.pop()
    with pytest.raises(ValueError, match="whole-gap"):
        encode(request, source)
    request, source = case()
    a = request.fixed_atoms.pop()
    del source[a.chain, a.residue_number, a.insertion_code, a.atom_name]
    with pytest.raises(ValueError, match="incomplete"):
        encode(request, source)


def test_context_is_not_silently_cropped_and_links_are_not_cut():
    request, source = case()
    extra = AtomIdentity("B", "1", "", "CA")
    request.fixed_atoms.append(extra)
    source["B", "1", "", "CA"] = {"line": "ATOM  extra"}
    with pytest.raises(ValueError, match="full source context"):
        encode(request, source)
    request, source = case()
    request.retained_explicit_links = ["unhandled"]
    with pytest.raises(ValueError, match="explicit links"):
        encode(request, source)


def test_native_atomic_preservation_is_not_claimed_by_frame_transport():
    request, source = case()
    _, receipt = encode(request, source)
    assert receipt["frame_contributor_atoms"] == 18
    assert receipt["psi_contributor_oxygens"] == 6
    assert receipt["not_directly_conditioned_source_atoms"] == 6
    assert "frames only" in receipt["native_preservation_contract"]


def test_frame_audit_handles_quaternion_sign_but_not_observed_motion():
    reference = np.array([[1., 0, 0, 0, 1, 2, 3], [1., 0, 0, 0, 4, 5, 6]])
    actual = reference.copy()
    actual[0, :4] *= -1
    actual[1, 4:] += 10
    assert compare_frames(reference, actual, np.array([True, False]))["translation_max_A"] == 0
    actual[0, 4] += .001
    with pytest.raises(ValueError, match="before any copying"):
        compare_frames(reference, actual, np.array([True, False]))


def test_projection_changes_only_declared_fixed_frames_without_relaxing_tolerance():
    reference = np.array([[1., 0, 0, 0, 1, 2, 3], [1., 0, 0, 0, 4, 5, 6]])
    current = reference[None].copy()
    current[0, 0, :4] += 1e-5
    current[0, 1, 4:] += 10
    original = current.copy()
    result = project_fixed_frames(current, reference, np.array([True, False]))
    assert np.array_equal(result[0, 0], reference[0])
    assert np.array_equal(result[0, 1], original[0, 1])
    assert np.array_equal(current, original)
    with pytest.raises(ValueError, match="axis"):
        project_fixed_frames(current, reference, np.array([1, 0]))
