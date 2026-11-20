import numpy as np
import pytest
from openmm import Context, HarmonicBondForce, OpenMMException, Platform, System, VerletIntegrator
from openmm.app import Topology, element
from scripts.anchor_refinement_guards import (
    immutable_context_admission,
    returned_state_guard,
    xh_equilibria,
)


def model():
    topology = Topology()
    chain = topology.addChain("B")
    ser = topology.addResidue("SER", chain, "161", insertionCode="A")
    water = topology.addResidue("HOH", chain, "12")
    oxygen = topology.addAtom("OG", element.oxygen, ser)
    hydrogen = topology.addAtom("HG", element.hydrogen, ser)
    topology.addAtom("O", element.oxygen, water)
    topology.addBond(oxygen, hydrogen)
    system = System()
    for mass in (0, 1, 0):
        system.addParticle(mass)
    force = HarmonicBondForce()
    force.addBond(0, 1, .096, 400000)
    system.addForce(force)
    points = np.array([[0., 0, 0], [.096, 0, 0], [.5, 0, 0]])
    return topology, system, points


def test_xh_uses_ff_equilibrium_and_reports_fixed_parent_blocker():
    topology, system, points = model()
    points[1] = [.3, 0, 0]
    rows = xh_equilibria(topology, system)
    assert len(rows) == 1 and rows[0]["equilibrium_nm"] == .096
    assert rows[0]["massless_indices"] == [0]
    assert rows[0]["parent"]["icode"] == "A"
    assert not rows[0]["openmm_hard_constraint_supported"]


def test_massless_constraint_is_not_supported_by_actual_engine():
    system = System()
    system.addParticle(0)
    system.addParticle(1)
    system.addConstraint(0, 1, .096)
    with pytest.raises(OpenMMException, match="massless"):
        Context(system, VerletIntegrator(.001), Platform.getPlatformByName("CPU"), {"Threads": "1"})


def test_immutable_contact_refuses_without_moving_or_deleting_water():
    topology, system, points = model()
    points[2] = [.1658382, 0, 0]
    before = points.copy()
    report = immutable_context_admission(topology, system, points)
    assert not report["admitted"]
    assert report["reasons"] == ["incompatible_immutable_context"]
    assert len(report["immutable_severe_clashes"]) == 1
    assert np.array_equal(points, before)
    system.setParticleMass(2, 16)
    assert immutable_context_admission(topology, system, points)["admitted"]


def test_unknown_context_is_partial_not_admitted():
    topology, system, points = model()
    residue = topology.addResidue("UNK", list(topology.chains())[0], "15")
    topology.addAtom("C1", element.carbon, residue)
    topology.addAtom("C2", element.carbon, residue)
    system.addParticle(0)
    system.addParticle(0)
    report = immutable_context_admission(topology, system, np.vstack([points, [1, 0, 0], [1.15, 0, 0]]))
    assert "unresolved_steric_connectivity" in report["reasons"]
    assert report["steric"]["status"] == "partial"


def test_fixed_reaction_forces_are_not_mobile_gradients():
    topology, system, points = model()
    forces = np.array([[1.e10]*3, [1.]*3, [-1.e10]*3])
    report = returned_state_guard(topology, system, points, points, forces, set())
    assert report["passed"]
    assert report["active_cartesian_component_rms"] == 1
    assert report["generated_heavy_component_rms"] is None


def test_fixed_particles_cannot_dilute_active_stationarity_failure():
    topology, system, points = model()
    forces = np.zeros_like(points)
    forces[1] = 12.
    report = returned_state_guard(topology, system, points, points, forces, set())
    assert report["active_cartesian_component_rms"] == 12
    assert report["massless_zeroed_all_component_rms"] < 10
    assert "nonstationary_active_cartesian_state" in report["reasons"]


def test_h_collapse_detected_across_residues_even_with_valid_parent_length():
    topology, system, points = model()
    points[2] = [.0965, 0, 0]
    report = returned_state_guard(topology, system, points, points, np.zeros_like(points), set())
    assert not report["xh_geometry_failures"]
    assert "hydrogen_atom_collapse" in report["reasons"]
    assert report["hydrogen_collapses"][0]["distance_nm"] == pytest.approx(.0005)


@pytest.mark.parametrize("length", [.01, .3])
def test_both_compressed_and_stretched_xh_bonds_refuse(length):
    topology, system, points = model()
    final = points.copy()
    final[1] = [length, 0, 0]
    report = returned_state_guard(topology, system, points, final, np.zeros_like(points), set())
    assert "covalent_xh_geometry_failure" in report["reasons"]


def test_nonfinite_and_fixed_motion_fail_closed():
    topology, system, points = model()
    final = points.copy()
    final[2, 0] += .001
    assert "native_fixed_observation_motion" in returned_state_guard(
        topology, system, points, final, np.zeros_like(points), set())["reasons"]
    final[1, 0] = np.nan
    assert returned_state_guard(topology, system, points, final, np.zeros_like(points), set())["reasons"] == ["nonfinite_returned_state"]


def test_missing_xh_force_parameters_and_projected_constraint_mode_are_not_guessed():
    topology, system, points = model()
    system.removeForce(0)
    with pytest.raises(ValueError, match="FF harmonic equilibrium"):
        xh_equilibria(topology, system)
    system.addConstraint(0, 1, .096)
    with pytest.raises(ValueError, match="constraint-projected"):
        returned_state_guard(topology, system, points, points, np.zeros_like(points), set())


def test_generated_thr_beta_is_checked_separately_from_ca_and_not_repaired():
    topology = Topology()
    residue = topology.addResidue("THR", topology.addChain("d"), "82", insertionCode="A")
    names = ("CB", "CA", "OG1", "CG2", "N", "C")
    points = np.array([(0, 0, 0), (.15, 0, 0), (-.045, .13, 0),
                       (-.05, -.07, -.11), (.15, 0, .145), (.15, .15, 0)])
    system = System()
    for name in names:
        topology.addAtom(name, element.get_by_symbol(name[0]), residue)
        system.addParticle(12)
    original = points.copy()
    report = returned_state_guard(topology, system, points, points, np.zeros_like(points), set(range(6)))
    assert "generated_stereocentre_geometry" in report["reasons"]
    assert not report["d_ca"]
    assert report["generated_stereocentre_findings"][0]["category"] == "thr_beta_chirality"
    assert np.array_equal(points, original)
