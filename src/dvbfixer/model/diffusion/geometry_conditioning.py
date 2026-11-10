"""Versioned geometry-first forces for new diffusion refinement requests.

Force-field equilibrium bonded terms close rings and junctions with observations
fixed. Bounded soft-core repulsion is then introduced before the unchanged full
force field. No raw bond length is promoted to an equilibrium constraint.
"""
from __future__ import annotations

from collections import defaultdict

from openmm import (
    CustomCompoundBondForce,
    CustomExternalForce,
    CustomNonbondedForce,
    NonbondedForce,
)

from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits, atom_residue

GEOMETRY_CONDITIONING_REVISION = "bonded-softcore-full-v1"
BONDED_GROUP = 0
SOFTCORE_GROUP = 30
FULL_NONBONDED_GROUP = 31
BOUNDARY_WALL_K = 1_000_000.0


def add_geometry_conditioning(system):
    """Return soft-core force or refuse unsupported nonbonded compositions."""
    forces = [f for f in system.getForces() if isinstance(f, NonbondedForce)]
    if len(forces) != 1:
        raise ValueError("geometry conditioning requires one canonical nonbonded force")
    if any(f.getForceGroup() != BONDED_GROUP for f in system.getForces()):
        raise ValueError("geometry conditioning cannot reinterpret existing force groups")
    full = forces[0]
    if full.getNonbondedMethod() != NonbondedForce.NoCutoff:
        raise ValueError("geometry conditioning requires the no-solvent NoCutoff refiner")
    full.setForceGroup(FULL_NONBONDED_GROUP)
    soft = CustomNonbondedForce(
        "soft_scale*epsilon*(sigma^2/(r^2+soft_alpha*sigma^2))^6;"
        "sigma=0.5*(sigma1+sigma2);epsilon=sqrt(epsilon1*epsilon2)")
    soft.addGlobalParameter("soft_scale", 0.0)
    soft.addGlobalParameter("soft_alpha", 1.0)
    soft.addPerParticleParameter("sigma")
    soft.addPerParticleParameter("epsilon")
    for index in range(full.getNumParticles()):
        _, sigma, epsilon = full.getParticleParameters(index)
        soft.addParticle([sigma, epsilon])
    for index in range(full.getNumExceptions()):
        first, second, *_ = full.getExceptionParameters(index)
        soft.addExclusion(first, second)
    soft.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    soft.setForceGroup(SOFTCORE_GROUP)
    system.addForce(soft)
    return soft


def add_anchor_boundary_walls(system, topology, restraints, limits: AnchorMovementLimits):
    """One-sided coordinate walls; final independent bounds remain authoritative."""
    if not restraints:
        return
    atoms = {(a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.name): a.index
             for a in topology.atoms()}
    force = CustomExternalForce(
        "0.5*k*max(0,sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2+1e-16)-radius)^2")
    for name in ("k", "x0", "y0", "z0", "radius"):
        force.addPerParticleParameter(name)
    grouped = defaultdict(dict)
    for restraint in restraints:
        atom = restraint.atom
        index = atoms[(atom.chain, atom.residue_number, atom.insertion_code, atom.atom_name)]
        radius = limits.heavy_displacement_angstrom_max
        if atom.atom_name == "CA":
            radius = min(radius, limits.ca_displacement_angstrom_max)
        force.addParticle(index, [BOUNDARY_WALL_K, *(x / 10 for x in restraint.reference_angstrom), radius / 10])
        grouped[atom_residue(atom)][atom.atom_name] = (index, restraint.reference_angstrom)
    system.addForce(force)
    distances = "+".join(f"((x{i}-a{i})^2+(y{i}-b{i})^2+(z{i}-c{i})^2)" for i in range(1, 5))
    backbone = CustomCompoundBondForce(4, f"0.5*k*max(0,sqrt(({distances})/4+1e-16)-radius)^2")
    backbone.addPerBondParameter("k")
    backbone.addPerBondParameter("radius")
    for i in range(1, 5):
        for axis in ("a", "b", "c"):
            backbone.addPerBondParameter(f"{axis}{i}")
    for group in grouped.values():
        if not {"N", "CA", "C", "O"} <= group.keys():
            raise ValueError("anchor boundary wall requires complete observed backbone")
        ordered = [group[name] for name in ("N", "CA", "C", "O")]
        backbone.addBond([index for index, _ in ordered],
            [BOUNDARY_WALL_K, limits.backbone_rmsd_angstrom_max / 10,
             *(x / 10 for _, xyz in ordered for x in xyz)])
    system.addForce(backbone)
