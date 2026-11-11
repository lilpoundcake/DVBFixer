"""Diagnostic severity must not decide which neighboring atoms can relax."""
from dataclasses import asdict, replace

import pytest

from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    DiffusionContractError,
    ResidueIdentity,
)
from dvbfixer.model.diffusion.gap_portfolio import digest
from dvbfixer.model.diffusion.joint_gap_research import (
    STRUCTURE_REVISION,
    build_joint_structure_view,
)
from dvbfixer.model.diffusion.protein_context import validate_cc94_request
from dvbfixer.model.diffusion.spatial_relaxation import (
    audit_spatial_relaxation,
    source_gap_components,
    spatial_mobility,
    spatial_restraints,
)
from tests.test_diffusion_anchor_acceptance import shift
from tests.test_diffusion_gap_portfolio import two_gaps
from tests.test_diffusion_joint_gap_research import two_gaps as context_request


def test_clean_geometry_still_includes_neighboring_observations():
    request, source, good = two_gaps()
    mobility = spatial_mobility(request, source, good)
    assert ResidueIdentity("C", "66") in mobility.observed_residues
    assert ResidueIdentity("C", "68") in mobility.observed_residues
    restraints = spatial_restraints(request, source, mobility)
    assert {r.atom for r in restraints} <= set(request.fixed_atoms)
    assert not {r.atom for r in restraints} & set(request.generated_atoms)
    assert asdict(AnchorMovementLimits()) == {
        "backbone_rmsd_angstrom_max": 0.5, "heavy_displacement_angstrom_max": 1.25,
        "ca_displacement_angstrom_max": 0.5}


def test_small_observed_motion_is_audited_without_requiring_a_fatal_defect():
    request, source, good = two_gaps()
    mobility, movements, assessment, passed = audit_spatial_relaxation(request, source, good, shift(good, "66", 0.04))
    assert passed
    assert not assessment.admission_failures
    assert any(m.residue == ResidueIdentity("C", "66") and m.ca_displacement_angstrom > 0 for m in movements)


def test_limits_and_metadata_are_not_relaxed_with_the_shell():
    request, source, good = two_gaps()
    _, movements, _, passed = audit_spatial_relaxation(request, source, good, shift(good, "66", 0.6))
    assert not passed and any(not m.within_limits for m in movements)
    with pytest.raises(DiffusionContractError, match="metadata"):
        audit_spatial_relaxation(request, source, good, good.replace("  1.00", "  0.50", 1))
    with pytest.raises(DiffusionContractError):
        spatial_mobility(request, source + "REMARK altered\n", good)


def test_source_groups_are_connected_components_not_one_pass_neighbor_lists():
    request, source, good = two_gaps()
    groups, edges = source_gap_components(request, source)
    assert groups == ((0, 1),)
    assert edges and edges[0]["minimum_flank_distance_angstrom"] <= 6.0


def test_structure_research_budget_does_not_change_production_admission(tmp_path):
    parent, source = context_request(tmp_path)
    original = parent.to_json()
    view = build_joint_structure_view(parent, source, seed=7)
    assert parent.to_json() == original
    assert view.request.protein_context.policy.name == STRUCTURE_REVISION
    assert view.request.protein_context.policy.residue_budget == 1024
    assert view.gap_indices == tuple(range(len(parent.gaps)))
    assert view.request.generated_atoms == parent.generated_atoms
    with pytest.raises(ValueError):
        validate_cc94_request(view.request)


def test_outside_observed_shell_is_exact_even_below_old_drift_tolerance(tmp_path):
    from tests.test_diffusion_input_preparation import FIXTURES, _request

    _, source, request = _request(tmp_path)
    # Synthetic ATOM-only test projection: deposited MASTER counts include waters
    # excluded by the helper, which would be an unrelated metadata admission error.
    source = ("".join(line for line in source.decode().splitlines(keepends=True)
                     if line.startswith("ATOM  ")) + "END\n").encode()
    request = replace(request, normalized_pdb=ArtifactReference("input.pdb", digest(source.decode())))
    good = "".join(line for line in (FIXTURES / "1UBQ.pdb").read_text().splitlines(keepends=True)
                   if line.startswith("ATOM  ")) + "END\n"
    mobility = spatial_mobility(request, source.decode(), good)
    selected = set(mobility.observed_residues)
    outside = next(a for a in request.fixed_atoms if a.atom_name == "CA"
        and ResidueIdentity(a.chain, a.residue_number, a.insertion_code) not in selected)
    changed = "".join(line[:30] + f"{float(line[30:38]) + 0.001:8.3f}" + line[38:]
        if line.startswith("ATOM  ") and line[21] == outside.chain
        and line[22:26].strip() == outside.residue_number and line[26].strip() == outside.insertion_code
        and line[12:16].strip() == "CA" else line for line in good.splitlines(keepends=True))
    with pytest.raises(DiffusionContractError, match="nonmobile"):
        audit_spatial_relaxation(request, source.decode(), good, changed)
