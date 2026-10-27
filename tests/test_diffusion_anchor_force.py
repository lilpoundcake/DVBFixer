"""Explicit observed-anchor force must preserve identities through topology edits."""

import numpy as np
import pytest

openmm = pytest.importorskip("openmm")
pytest.importorskip("pdbfixer")

from openmm.app import Topology, element
from openmm.unit import kilojoule_per_mole, nanometer

from dvbfixer.model.diffusion.anchor_policy import AnchorPositionRestraint, NumericalRefinementError
from dvbfixer.model.diffusion.boundary_refinement import (
    BoundaryNumericalError,
    BoundaryRefinementError,
    _build_anchor_position_force,
    refine_generated_region,
)
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def test_final_topology_indices_case_and_icode_with_physical_energy():
    topology = Topology()
    system = openmm.System()
    for chain_id, code in (("w", ""), ("W", ""), ("W", "A")):
        chain = topology.addChain(chain_id)
        r = topology.addResidue("GLY", chain, id="29", insertionCode=code)
        topology.addAtom("CA", element.carbon, r)
        system.addParticle(12)
    chosen = AtomIdentity("W", "29", "A", "CA")
    restraint = AnchorPositionRestraint(chosen, (10, 20, 30), 1000)
    force = _build_anchor_position_force(topology, (restraint,))
    index, params = force.getParticleParameters(0)
    assert index == 2 and list(params) == pytest.approx([1000, 1, 2, 3])
    system.addForce(force)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(np.array([[50, 50, 50], [90, 90, 90], [1, 2, 3]]) * nanometer)
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) == pytest.approx(0)
    context.setPositions(np.array([[50, 50, 50], [90, 90, 90], [1.1, 2, 3]]) * nanometer)
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) == pytest.approx(5)


def test_missing_and_duplicate_final_restraint_identities_rejected():
    topology = Topology()
    restraint = AnchorPositionRestraint(AtomIdentity("W", "29", "", "CA"), (0, 0, 0), 1000)
    with pytest.raises(BoundaryRefinementError, match="missing"):
        _build_anchor_position_force(topology, (restraint,))
    with pytest.raises(BoundaryRefinementError, match="duplicate"):
        _build_anchor_position_force(topology, (restraint, restraint))


def test_generated_residues_cannot_be_reclassified_as_observed_anchors():
    atom = AtomIdentity("W", "29", "", "CA")
    with pytest.raises(BoundaryRefinementError, match="cannot target generated"):
        refine_generated_region("", generated_residues=(ResidueIdentity("W", "29"),), generated_atoms=(atom,),
                                observed_anchor_restraints=(AnchorPositionRestraint(atom, (0, 0, 0), 1000),))


def test_numerical_error_is_backward_compatible_with_existing_refiner_handlers():
    error = BoundaryNumericalError("No finite all-L minimum")
    assert isinstance(error, BoundaryRefinementError)
    assert isinstance(error, NumericalRefinementError)
