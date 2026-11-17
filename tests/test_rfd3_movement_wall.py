"""Actual Reference force control, identity re-resolution and fail-closed masks."""
import numpy as np
import pytest
from openmm import Context, Platform, System, VerletIntegrator
from openmm.app import Topology, element
from openmm.unit import kilojoule_per_mole, nanometer

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError
from dvbfixer.model.diffusion.rfd3_movement_wall import build_missing_atom_wall


def fixture():
    top = Topology()
    chain = top.addChain("a")
    res = top.addResidue("GLY", chain, "1", insertionCode="A")
    top.addAtom("H", element.hydrogen, res)
    top.addAtom("CA", element.carbon, res)
    identity = AtomIdentity("a", "1", "A", "CA")
    return top, identity


def test_actual_wall_has_zero_interior_energy_and_inward_exterior_force():
    top, atom = fixture()
    force = build_missing_atom_wall(top, {atom: (0., 0., 0.)}, mobile=frozenset({atom}),
        outer_limit_angstrom=3., buffer_angstrom=.2, strength_kj_mol_nm2=1e6)
    assert force.getParticleParameters(0)[0] == 1
    system = System()
    for _ in range(2):
        system.addParticle(12)
    system.addForce(force)
    integrator = VerletIntegrator(.001)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPositions([(0., 0., 0.), (.27, 0., 0.)])
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) == 0
    context.setPositions([(0., 0., 0.), (.31, 0., 0.)])
    state = context.getState(getEnergy=True, getForces=True)
    assert state.getPotentialEnergy().value_in_unit(kilojoule_per_mole) == pytest.approx(450.)
    forces = state.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole / nanometer)
    assert forces[1, 0] < 0 and np.all(forces[0] == 0)


@pytest.mark.parametrize("problem", ["observed-reference", "wrong-chain-case", "no-buffer", "excess-buffer", "strength", "hydrogen"])
def test_invalid_or_inexact_wall_domain_is_refused(problem):
    top, atom = fixture()
    references = {atom: (0., 0., 0.)}
    mobile, buffer, strength = frozenset({atom}), .2, 1e6
    if problem == "observed-reference":
        references[AtomIdentity("A", "1", "A", "CA")] = (0., 0., 0.)
    elif problem == "wrong-chain-case":
        atom = AtomIdentity("A", "1", "A", "CA")
        references, mobile = {atom: (0., 0., 0.)}, frozenset({atom})
    elif problem == "no-buffer":
        buffer = 0.
    elif problem == "excess-buffer":
        buffer = 3.
    elif problem == "strength":
        strength = True
    else:
        atom = AtomIdentity("a", "1", "A", "H")
        references, mobile = {atom: (0., 0., 0.)}, frozenset({atom})
    with pytest.raises(CustodyError):
        build_missing_atom_wall(top, references, mobile=mobile, outer_limit_angstrom=3.,
                                buffer_angstrom=buffer, strength_kj_mol_nm2=strength)


def test_wall_survives_real_final_topology_and_all_geometry_conditioning_phases(monkeypatch):
    from pathlib import Path

    from dvbfixer.model.diffusion import boundary_refinement as kernel
    from dvbfixer.model.diffusion.contract import ResidueIdentity

    source = Path(__file__).parent / "fixtures/8cz8/8cz8_a_u.pdb"
    lines = [line for line in source.read_text().splitlines(keepends=True)
             if line.startswith("ATOM  ") and line[21] == "C" and 6 <= int(line[22:26]) <= 8
             and line[76:78].strip() != "H"]
    mobile = frozenset(AtomIdentity("C", "7", "", line[12:16].strip())
                       for line in lines if int(line[22:26]) == 7)
    freeze = kernel._freeze_observed_heavy_atoms
    recorded = []

    def add_wall(system, topology, original, requested):
        freeze(system, topology, original, requested)
        force = build_missing_atom_wall(topology, {a: tuple(original[a]) for a in mobile}, mobile=mobile,
            outer_limit_angstrom=3., buffer_angstrom=.2, strength_kj_mol_nm2=1e6)
        system.addForce(force)
        recorded.append(force.getNumParticles())

    monkeypatch.setattr(kernel, "_freeze_observed_heavy_atoms", add_wall)
    result = kernel.refine_generated_region("".join(lines) + "END\n",
        generated_residues=(ResidueIdentity("C", "7"),), generated_atoms=tuple(sorted(mobile)),
        platform_name="CPU", max_iterations=10, restart_count=0, random_seed=101, geometry_preconditioning=True)
    assert result.coordinates_angstrom.keys() == mobile
    assert recorded == [len(mobile)]
    assert {row["phase"] for row in result.search_diagnostics} >= {
        "bonded-geometry", "soft-steric-1", "soft-steric-2", "full-forcefield"}
