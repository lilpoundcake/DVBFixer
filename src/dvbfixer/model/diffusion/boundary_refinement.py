"""Localized post-sampling refinement for generated protein residues."""

from __future__ import annotations

import hashlib
import math
import random
from dataclasses import dataclass
from io import StringIO
from typing import Any

import numpy as np
from openmm import (
    CustomAngleForce,
    CustomExternalForce,
    CustomTorsionForce,
    OpenMMException,
    Platform,
    VerletIntegrator,
)
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile, Simulation
from openmm.unit import kilojoule_per_mole, nanometer, picosecond
from pdbfixer import PDBFixer

from dvbfixer.diagnose.chemistry import (
    check_backbone_bond_angles,
    check_bond_lengths,
    check_peptide_omegas,
)
from dvbfixer.diagnose.report import Severity
from dvbfixer.diagnose.steric import clashes_python
from dvbfixer.ffutils import FF_ALIASES
from dvbfixer.ffutils.geometry import (
    assert_all_l,
    build_ca_chirality_force,
    find_d_residues,
    fix_ca_chirality,
    rebuild_missing_atoms_with_retry,
    repair_misplaced_hydrogens,
)
from dvbfixer.model.diffusion.anchor_policy import AnchorPositionRestraint, NumericalRefinementError
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity

BOUNDARY_REFINEMENT_REVISION = "dvbfixer-openmm-boundary-refinement-v8"
BOUNDARY_REFINEMENT_FORCEFIELD = tuple(FF_ALIASES["amber"])
BOUNDARY_REFINEMENT_MAX_ITERATIONS = 250
BOUNDARY_REFINEMENT_TOLERANCE_KJ_MOL_NM = 10.0
BOUNDARY_REFINEMENT_RESTART_COUNT = 8
BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM = 0.75
BOUNDARY_REFINEMENT_OMEGA_KJ_MOL = 500.0
BOUNDARY_REFINEMENT_PEPTIDE_ANGLE_KJ_MOL_RAD2 = 1000.0
_REFINEMENT_BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})
_FIXED_ATOM_INTERNAL_TOLERANCE_ANGSTROM = 1.0e-4


class BoundaryRefinementError(RuntimeError):
    """Raised when localized refinement cannot preserve its hard invariants."""


class BoundaryNumericalError(BoundaryRefinementError, NumericalRefinementError):
    """No finite all-L minimum was found; eligible for typed numerical recovery."""


@dataclass(frozen=True, slots=True)
class BoundaryRefinementResult:
    coordinates_angstrom: dict[AtomIdentity, np.ndarray]
    chirality_repairs: tuple[ResidueIdentity, ...]
    pre_coordinate_sha256: str
    post_coordinate_sha256: str
    initial_energy_kj_mol: float
    final_energy_kj_mol: float
    platform: str


def refine_generated_region(
    pdb_text: str,
    *,
    generated_residues: tuple[ResidueIdentity, ...],
    generated_atoms: tuple[AtomIdentity, ...],
    max_iterations: int = BOUNDARY_REFINEMENT_MAX_ITERATIONS,
    restart_count: int = BOUNDARY_REFINEMENT_RESTART_COUNT,
    perturbation_angstrom: float = BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM,
    platform_name: str = "Reference",
    random_seed: int = 1,
    restrain_amides: bool = True,
    observed_anchor_restraints: tuple[AnchorPositionRestraint, ...] = (),
    mobile_position_restraints: tuple[AnchorPositionRestraint, ...] = (),
) -> BoundaryRefinementResult:
    """Minimize generated atoms and explicitly restrained observed anchor atoms.

    The default remains generated-only. Observed mobility is a separate explicit
    scientific argument, not a reclassification of observed atoms as generated.
    The caller must independently validate anchor eligibility, drift and geometry;
    this low-level method grants no publication exception. Returned coordinates
    cover generated atoms AND authorized observed anchor atoms, nothing else.
    """
    if max_iterations <= 0:
        raise BoundaryRefinementError("refinement max_iterations must be positive")
    if restart_count < 0:
        raise BoundaryRefinementError("refinement restart_count must be non-negative")
    if not math.isfinite(perturbation_angstrom) or perturbation_angstrom <= 0:
        raise BoundaryRefinementError(
            "refinement perturbation_angstrom must be finite and positive"
        )
    if random_seed < 0:
        raise BoundaryRefinementError("refinement random_seed must be non-negative")
    generated_residue_set = set(generated_residues)
    if not generated_residue_set or len(generated_residue_set) != len(generated_residues):
        raise BoundaryRefinementError(
            "refinement generated residues must be non-empty and unique"
        )
    generated_atom_set = set(generated_atoms)
    if not generated_atom_set or len(generated_atom_set) != len(generated_atoms):
        raise BoundaryRefinementError(
            "refinement generated atoms must be non-empty and unique"
        )
    if any(_atom_residue(atom) not in generated_residue_set for atom in generated_atoms):
        raise BoundaryRefinementError(
            "refinement generated atoms must belong to generated residues"
        )

    anchor_atoms = tuple(restraint.atom for restraint in observed_anchor_restraints)
    if len(set(anchor_atoms)) != len(anchor_atoms):
        raise BoundaryRefinementError("duplicate observed anchor restraint")
    if any(_atom_residue(atom) in generated_residue_set for atom in anchor_atoms):
        raise BoundaryRefinementError("observed anchor restraints cannot target generated residues")
    mobile_atoms_requested = (*generated_atoms, *anchor_atoms)
    mobile_atom_set = set(mobile_atoms_requested)
    mobile_residue_set = generated_residue_set | {_atom_residue(atom) for atom in anchor_atoms}
    if (len({r.atom for r in mobile_position_restraints}) != len(mobile_position_restraints)
            or any(r.atom not in generated_atom_set for r in mobile_position_restraints)):
        raise BoundaryRefinementError("mobile positional restraints must uniquely target mobile atoms")

    pdb = PDBFile(StringIO(pdb_text))
    original_coordinates = _coordinates_by_identity(pdb.topology, pdb.positions)
    original_heavy = {AtomIdentity(atom.residue.chain.id, str(atom.residue.id), atom.residue.insertionCode.strip(), atom.name) for atom in pdb.topology.atoms()
        if atom.element is None or atom.element.symbol != "H"}
    original_coordinates = {atom: xyz for atom, xyz in original_coordinates.items() if atom in original_heavy}
    missing = mobile_atom_set - set(original_coordinates)
    if missing:
        raise BoundaryRefinementError(
            "refinement input is missing requested generated atoms: "
            + _atom_labels(missing)
        )

    fixed_d = _d_residue_identities(pdb.topology, pdb.positions) - generated_residue_set
    if fixed_d:
        raise BoundaryRefinementError(
            "refinement refuses to alter D-Cα geometry outside generated residues: "
            + _residue_labels(fixed_d)
        )

    try:
        platform = Platform.getPlatformByName(platform_name)
    except Exception as exc:
        raise BoundaryRefinementError(
            f"OpenMM platform {platform_name!r} is unavailable"
        ) from exc

    # A normalized source may contain unrelated incomplete side chains. Complete
    # them only in this temporary force-field topology; publication still copies
    # coordinates exclusively for the requested generated atoms.
    # Hydrogens are disposable only in the temporary FF system. Publication
    # still rewrites generated coordinates exclusively, preserving source lines.
    temporary = Modeller(pdb.topology, pdb.positions)
    temporary.delete([atom for atom in temporary.topology.atoms()
        if atom.element is not None and atom.element.symbol == "H"])
    heavy_text = StringIO()
    PDBFile.writeFile(temporary.topology, temporary.positions, heavy_text, keepIds=True)
    fixer = PDBFixer(pdbfile=StringIO(heavy_text.getvalue()))
    fixer.platform = platform
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    rebuild_missing_atoms_with_retry(fixer)
    topology = fixer.topology
    positions = fixer.positions
    rebuilt_d = _d_residue_identities(topology, positions)
    rebuilt_fixed_d = rebuilt_d - generated_residue_set
    unsafe_fixed_d = {
        identity
        for identity in rebuilt_fixed_d
        if _has_observed_sidechain(identity, topology, original_coordinates)
    }
    if unsafe_fixed_d:
        raise BoundaryRefinementError(
            "temporary completion produced D-Cα geometry whose repair would move "
            "observed atoms outside generated residues: "
            + _residue_labels(unsafe_fixed_d)
        )
    generated_d = rebuilt_d & generated_residue_set
    repairable_d = generated_d | rebuilt_fixed_d
    repairs = (
        fix_ca_chirality(
            topology,
            positions,
            residue_identities={
                (identity.chain, identity.residue_number, identity.insertion_code)
                for identity in repairable_d
            },
        )
        if repairable_d
        else 0
    )
    if repairs != len(repairable_d):
        raise BoundaryRefinementError(
            "temporary D-Cα repair did not repair every eligible offender"
        )
    assert_all_l(topology, positions)

    forcefield = ForceField(*BOUNDARY_REFINEMENT_FORCEFIELD)
    modeller = Modeller(topology, positions)
    random_state = random.getstate()
    random.seed(random_seed)
    try:
        modeller.addHydrogens(forcefield, pH=7.0, platform=platform)
    finally:
        random.setstate(random_state)
    repair_misplaced_hydrogens(modeller.topology, modeller.positions)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
        rigidWater=False,
    )
    # Use the force field's equilibrium bond lengths. Constraining generated
    # bonds to their raw sampled lengths freezes the very defects being repaired.
    _freeze_observed_heavy_atoms(
        system, modeller.topology, original_coordinates, mobile_atom_set,
    )
    position_restraints = (*observed_anchor_restraints, *mobile_position_restraints)
    if position_restraints:
        system.addForce(_build_anchor_position_force(modeller.topology, position_restraints))
    chirality_force, protected_centres = build_ca_chirality_force(
        modeller.topology,
        modeller.positions,
    )
    if protected_centres:
        system.addForce(chirality_force)
    if restrain_amides:
        omega_force, restrained_amides = _build_trans_amide_force(
            modeller.topology,
            mobile_residue_set,
        )
        if restrained_amides:
            system.addForce(omega_force)
    peptide_angle_force, restrained_angles = _build_peptide_angle_force(
        modeller.topology,
        mobile_residue_set,
    )
    if restrained_angles:
        system.addForce(peptide_angle_force)

    if platform_name == "CPU":
        properties = {"Threads": "1"}
    elif platform_name == "CUDA":
        properties = {"Precision": "mixed", "DeterministicForces": "true"}
    else:
        properties = {}
    integrator = VerletIntegrator(0.001 * picosecond)
    simulation = Simulation(
        modeller.topology,
        system,
        integrator,
        platform,
        properties,
    )
    simulation.context.setPositions(modeller.positions)
    initial_state = simulation.context.getState(getEnergy=True)
    try:
        simulation.minimizeEnergy(
            tolerance=BOUNDARY_REFINEMENT_TOLERANCE_KJ_MOL_NM
            * kilojoule_per_mole / nanometer, maxIterations=max_iterations,
        )
        first_state = simulation.context.getState(getEnergy=True, getPositions=True)
    except OpenMMException:
        first_state = None
    best_state = first_state if _acceptable_minimum(modeller.topology, first_state) else None
    best_score = (_local_geometry_error_count(modeller.topology, best_state.getPositions(), mobile_residue_set)
        if best_state is not None else float("inf"))
    restart_positions = np.asarray(
        (best_state.getPositions(asNumpy=True) if best_state is not None else modeller.positions).value_in_unit(nanometer),
        dtype=np.float64,
    )
    mobile_atoms = [
        atom
        for atom in modeller.topology.atoms()
        if (_topology_residue(atom.residue) in mobile_residue_set
            and (atom.element is not None and atom.element.symbol == "H"
                 or AtomIdentity(atom.residue.chain.id, str(atom.residue.id),
                                 atom.residue.insertionCode.strip(), atom.name) in mobile_atom_set))
    ]
    if restart_count:
        for cycle in range(restart_count):
            rng = np.random.default_rng(random_seed + cycle)
            perturbed = restart_positions.copy()
            for atom in mobile_atoms:
                scale_nm = perturbation_angstrom / 10.0
                if atom.name in {"N", "CA", "C", "O"}:
                    scale_nm *= 0.25
                perturbed[atom.index] += rng.normal(0.0, scale_nm, size=3)
            simulation.context.setPositions(perturbed * nanometer)
            try:
                simulation.minimizeEnergy(
                    tolerance=BOUNDARY_REFINEMENT_TOLERANCE_KJ_MOL_NM
                    * kilojoule_per_mole
                    / nanometer,
                    maxIterations=max_iterations,
                )
            except OpenMMException:
                continue
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            if not _acceptable_minimum(modeller.topology, state):
                continue
            score = _local_geometry_error_count(
                modeller.topology,
                state.getPositions(),
                mobile_residue_set,
            )
            state_energy = float(
                state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            )
            best_energy = float(
                best_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            ) if best_state is not None else float("inf")
            if (score, state_energy) < (best_score, best_energy):
                best_state = state
                best_score = score
            if best_score == 0:
                break
    if best_state is None:
        raise BoundaryNumericalError("localized refinement found no finite all-L minimum after seeded retries")
    final_state = best_state
    refined_positions = final_state.getPositions()
    initial_energy = float(
        initial_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    )
    final_energy = float(
        final_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    )
    if not math.isfinite(initial_energy) or not math.isfinite(final_energy):
        raise BoundaryNumericalError("localized refinement produced non-finite energy")
    assert_all_l(modeller.topology, refined_positions)

    refined_coordinates = _coordinates_by_identity(
        modeller.topology,
        refined_positions,
    )
    for identity, original in original_coordinates.items():
        if identity in mobile_atom_set:
            continue
        refined = refined_coordinates.get(identity)
        # CUDA's mixed-precision coordinate round-trip perturbs massless fixed
        # particles below PDB precision. Publication still copies their source
        # records unchanged; reject any displacement beyond numerical noise.
        if refined is None or not np.allclose(
            refined,
            original,
            rtol=0.0,
            atol=_FIXED_ATOM_INTERNAL_TOLERANCE_ANGSTROM,
        ):
            raise BoundaryRefinementError(
                "localized refinement moved or removed a fixed atom: "
                f"{_atom_label(identity)}"
            )

    result_coordinates = {
        identity: refined_coordinates[identity]
        for identity in mobile_atoms_requested
        if identity in refined_coordinates
    }
    missing = mobile_atom_set - set(result_coordinates)
    if missing:
        raise BoundaryRefinementError(
            "localized refinement lost generated atoms: " + _atom_labels(missing)
        )
    return BoundaryRefinementResult(
        coordinates_angstrom=result_coordinates,
        chirality_repairs=tuple(sorted(repairable_d)),
        pre_coordinate_sha256=_coordinate_digest(
            {identity: original_coordinates[identity] for identity in mobile_atoms_requested}
        ),
        post_coordinate_sha256=_coordinate_digest(result_coordinates),
        initial_energy_kj_mol=initial_energy,
        final_energy_kj_mol=final_energy,
        platform=platform.getName(),
    )


def refine_observed_region(
    pdb_text: str,
    *,
    residues: tuple[ResidueIdentity, ...],
    atoms: tuple[AtomIdentity, ...],
    restraints: tuple[AnchorPositionRestraint, ...],
    max_iterations: int = 1000,
) -> BoundaryRefinementResult:
    """Use the same local minimizer *before* freezing a diffusion request.

    The low-level generated_* arguments describe mobility, not provenance. This
    adapter never modifies a request, generated masks or template ownership.
    All measured atoms outside ``atoms`` are massless/fixed; all selected heavy
    atoms are restrained. No solvent, stochastic restarts or unrestrained repair.
    """
    if {restraint.atom for restraint in restraints} != set(atoms):
        raise BoundaryRefinementError("observed repair requires restraints on every selected atom")
    return refine_generated_region(
        pdb_text, generated_residues=residues, generated_atoms=atoms,
        mobile_position_restraints=restraints, max_iterations=max_iterations,
        restart_count=0, platform_name="CPU", random_seed=1, restrain_amides=False,
    )


def _build_anchor_position_force(
    topology: Any,
    restraints: tuple[AnchorPositionRestraint, ...],
) -> CustomExternalForce:
    """Resolve stable identities only after final topology rebuilding; Å -> nm."""
    expected = {restraint.atom: restraint for restraint in restraints}
    if len(expected) != len(restraints):
        raise BoundaryRefinementError("duplicate observed anchor restraint")
    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    for name in ("k", "x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    seen: set[AtomIdentity] = set()
    for atom in topology.atoms():
        identity = AtomIdentity(atom.residue.chain.id, str(atom.residue.id),
                                atom.residue.insertionCode.strip(), atom.name)
        if identity not in expected:
            continue
        if identity in seen:
            raise BoundaryRefinementError("ambiguous observed anchor topology identity")
        seen.add(identity)
        restraint = expected[identity]
        force.addParticle(atom.index, [restraint.k_kj_mol_nm2,
                                      *(x / 10 for x in restraint.reference_angstrom)])
    if seen != expected.keys():
        raise BoundaryRefinementError("observed anchor missing from final topology")
    return force


def _freeze_observed_heavy_atoms(
    system: Any,
    topology: Any,
    original_coordinates: dict[AtomIdentity, np.ndarray],
    generated_atoms: set[AtomIdentity],
) -> None:
    """Freeze observed heavy atoms, not disposable temporary FF completions.

    All input hydrogens were stripped and rebuilt. Freezing those new positions
    can lock enormous steric energies into the objective. Temporary missing heavy
    atoms likewise need to relax against the immutable observed scaffold. Neither
    kind of completion is copied into the published source records.
    """
    fixed_atoms = original_coordinates.keys() - generated_atoms
    for atom in topology.atoms():
        identity = AtomIdentity(
            atom.residue.chain.id, str(atom.residue.id),
            atom.residue.insertionCode.strip(), atom.name,
        )
        if identity in fixed_atoms:
            system.setParticleMass(atom.index, 0.0)


def _coordinates_by_identity(topology: Any, positions: Any) -> dict[AtomIdentity, np.ndarray]:
    coordinates: dict[AtomIdentity, np.ndarray] = {}
    for atom in topology.atoms():
        identity = AtomIdentity(
            atom.residue.chain.id,
            str(atom.residue.id),
            atom.residue.insertionCode.strip(),
            atom.name,
        )
        if identity in coordinates:
            raise BoundaryRefinementError(
                "refinement topology contains duplicate atom identity: "
                f"{_atom_label(identity)}"
            )
        coordinates[identity] = np.asarray(
            positions[atom.index].value_in_unit(nanometer),
            dtype=np.float64,
        ) * 10.0
    return coordinates


def _state_is_finite(state: Any) -> bool:
    positions = np.asarray(
        state.getPositions(asNumpy=True).value_in_unit(nanometer),
        dtype=np.float64,
    )
    energy = float(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
    return bool(np.isfinite(positions).all() and math.isfinite(energy))


def _acceptable_minimum(topology: Any, state: Any) -> bool:
    """Chirality is an admission gate, never a tie-breaker against energy."""
    return state is not None and _state_is_finite(state) and not find_d_residues(topology, state.getPositions())


def _build_trans_amide_force(
    topology: Any,
    generated_residues: set[ResidueIdentity],
) -> tuple[CustomTorsionForce, int]:
    force = CustomTorsionForce("omega_k*(1+cos(theta))")
    force.addGlobalParameter("omega_k", BOUNDARY_REFINEMENT_OMEGA_KJ_MOL)
    count = 0
    for chain in topology.chains():
        residues = list(chain.residues())
        for left, right in zip(residues, residues[1:]):
            if not ({_topology_residue(left), _topology_residue(right)} & generated_residues):
                continue
            left_atoms = {atom.name: atom for atom in left.atoms()}
            right_atoms = {atom.name: atom for atom in right.atoms()}
            if not {"CA", "C"} <= left_atoms.keys() or not {"N", "CA"} <= right_atoms.keys():
                continue
            force.addTorsion(
                left_atoms["CA"].index,
                left_atoms["C"].index,
                right_atoms["N"].index,
                right_atoms["CA"].index,
                (),
            )
            count += 1
    return force, count


def _build_peptide_angle_force(
    topology: Any,
    generated_residues: set[ResidueIdentity],
) -> tuple[CustomAngleForce, int]:
    force = CustomAngleForce("0.5*angle_k*(theta-theta0)^2")
    force.addGlobalParameter(
        "angle_k",
        BOUNDARY_REFINEMENT_PEPTIDE_ANGLE_KJ_MOL_RAD2,
    )
    force.addPerAngleParameter("theta0")
    count = 0
    targets = (
        ("CA", "C", "N", math.radians(117.5)),
        ("O", "C", "N", math.radians(122.5)),
        ("C", "N", "CA", math.radians(122.5)),
    )
    for chain in topology.chains():
        residues = list(chain.residues())
        for left, right in zip(residues, residues[1:]):
            if not ({_topology_residue(left), _topology_residue(right)} & generated_residues):
                continue
            left_atoms = {atom.name: atom for atom in left.atoms()}
            right_atoms = {atom.name: atom for atom in right.atoms()}
            for first, centre, third, target in targets:
                if third == "N":
                    atom1 = left_atoms.get(first)
                    atom2 = left_atoms.get(centre)
                    atom3 = right_atoms.get(third)
                else:
                    atom1 = left_atoms.get(first)
                    atom2 = right_atoms.get(centre)
                    atom3 = right_atoms.get(third)
                if atom1 is None or atom2 is None or atom3 is None:
                    continue
                force.addAngle(atom1.index, atom2.index, atom3.index, (target,))
                count += 1
    return force, count


def _constrain_generated_heavy_bonds(
    system: Any,
    topology: Any,
    positions: Any,
    generated_residues: set[ResidueIdentity],
) -> None:
    for atom1, atom2 in topology.bonds():
        residue1 = _topology_residue(atom1.residue)
        residue2 = _topology_residue(atom2.residue)
        if residue1 != residue2 or residue1 not in generated_residues:
            continue
        if atom1.element.symbol == "H" or atom2.element.symbol == "H":
            continue
        first = np.asarray(
            positions[atom1.index].value_in_unit(nanometer),
            dtype=np.float64,
        )
        second = np.asarray(
            positions[atom2.index].value_in_unit(nanometer),
            dtype=np.float64,
        )
        system.addConstraint(
            atom1.index,
            atom2.index,
            float(np.linalg.norm(first - second)) * nanometer,
        )


def _local_geometry_error_count(
    topology: Any,
    positions: Any,
    generated_residues: set[ResidueIdentity],
) -> int:
    local = set(generated_residues)
    for chain in topology.chains():
        residues = list(chain.residues())
        for index, residue in enumerate(residues):
            if _topology_residue(residue) not in generated_residues:
                continue
            if index:
                local.add(_topology_residue(residues[index - 1]))
            if index + 1 < len(residues):
                local.add(_topology_residue(residues[index + 1]))
    keys = {(residue.chain, residue.residue_number + residue.insertion_code) for residue in local}
    errors = sum(
        finding.severity is Severity.ERROR
        and (finding.chain, finding.resid) in keys
        for finding in check_bond_lengths(topology, positions)
    )
    errors += sum(
        finding.severity is Severity.ERROR
        and (finding.chain, finding.resid) in keys
        for finding in check_backbone_bond_angles(topology, positions)
    )
    errors += sum(
        finding.category == "non_planar_amide"
        and (finding.chain, finding.resid) in keys
        for finding in check_peptide_omegas(topology, positions)
    )
    for finding in clashes_python(
        topology,
        positions,
        clash_warn_a=0.9,
        clash_error_a=0.9,
    ):
        if finding.severity is not Severity.ERROR:
            continue
        partner = finding.extra.get("clash_partner", {})
        if (finding.chain, finding.resid) in keys or (
            partner.get("chain"),
            partner.get("resid"),
        ) in keys:
            errors += 1
    return errors


def _d_residue_identities(topology: Any, positions: Any) -> set[ResidueIdentity]:
    offenders = find_d_residues(topology, positions)
    identities: set[ResidueIdentity] = set()
    for chain, residue_number, residue_name, _triple in offenders:
        matches = [
            residue
            for residue in topology.residues()
            if residue.chain.id == chain
            and str(residue.id) == residue_number
            and residue.name == residue_name
        ]
        if len(matches) != 1:
            raise BoundaryRefinementError(
                "cannot resolve insertion-code-aware identity for D-Cα residue "
                f"{chain}/{residue_name}{residue_number}"
            )
        identities.add(_topology_residue(matches[0]))
    return identities


def _topology_residue(residue: Any) -> ResidueIdentity:
    return ResidueIdentity(
        residue.chain.id,
        str(residue.id),
        residue.insertionCode.strip(),
    )


def _atom_residue(atom: AtomIdentity) -> ResidueIdentity:
    return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)


def _has_observed_sidechain(
    residue: ResidueIdentity,
    topology: Any,
    original_coordinates: dict[AtomIdentity, np.ndarray],
) -> bool:
    matches = [item for item in topology.residues() if _topology_residue(item) == residue]
    if len(matches) != 1:
        raise BoundaryRefinementError(
            "cannot resolve fixed residue after temporary completion: "
            + _residue_labels({residue})
        )
    return any(
        atom.atom_name not in _REFINEMENT_BACKBONE_ATOMS
        for atom in original_coordinates
        if _atom_residue(atom) == residue
    )


def _coordinate_digest(coordinates: dict[AtomIdentity, np.ndarray]) -> str:
    digest = hashlib.sha256()
    for identity in sorted(coordinates):
        digest.update(
            "\0".join(
                (
                    identity.chain,
                    identity.residue_number,
                    identity.insertion_code,
                    identity.atom_name,
                )
            ).encode("utf-8")
        )
        digest.update(np.asarray(coordinates[identity], dtype=">f8").tobytes())
    return digest.hexdigest()


def _atom_label(identity: AtomIdentity) -> str:
    return (
        f"{identity.chain}/{identity.residue_number}{identity.insertion_code}/"
        f"{identity.atom_name}"
    )


def _atom_labels(identities: set[AtomIdentity]) -> str:
    return ", ".join(_atom_label(identity) for identity in sorted(identities))


def _residue_labels(identities: set[ResidueIdentity]) -> str:
    return ", ".join(
        f"{identity.chain}/{identity.residue_number}{identity.insertion_code}"
        for identity in sorted(identities)
    )
