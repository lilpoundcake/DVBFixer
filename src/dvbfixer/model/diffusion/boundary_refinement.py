"""Localized post-sampling refinement for generated protein residues."""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from io import StringIO

import numpy as np
from openmm import (
    CustomAngleForce,
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
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity

BOUNDARY_REFINEMENT_REVISION = "dvbfixer-openmm-boundary-refinement-v2"
BOUNDARY_REFINEMENT_FORCEFIELD = tuple(FF_ALIASES["amber"])
BOUNDARY_REFINEMENT_MAX_ITERATIONS = 2000
BOUNDARY_REFINEMENT_TOLERANCE_KJ_MOL_NM = 10.0
BOUNDARY_REFINEMENT_RESTART_COUNT = 8
BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM = 0.75
BOUNDARY_REFINEMENT_OMEGA_KJ_MOL = 500.0
BOUNDARY_REFINEMENT_PEPTIDE_ANGLE_KJ_MOL_RAD2 = 1000.0


class BoundaryRefinementError(RuntimeError):
    """Raised when localized refinement cannot preserve its hard invariants."""


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
    platform_name: str = "CPU",
    random_seed: int = 1,
    restrain_amides: bool = True,
) -> BoundaryRefinementResult:
    """Minimize generated residues while keeping every source atom exact."""
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

    pdb = PDBFile(StringIO(pdb_text))
    original_coordinates = _coordinates_by_identity(pdb.topology, pdb.positions)
    missing = generated_atom_set - set(original_coordinates)
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

    # A normalized source may contain unrelated incomplete side chains. Complete
    # them only in this temporary force-field topology; publication still copies
    # coordinates exclusively for the requested generated atoms.
    fixer = PDBFixer(pdbfile=StringIO(pdb_text))
    fixer.platform = Platform.getPlatformByName("CPU")
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    rebuild_missing_atoms_with_retry(fixer)
    topology = fixer.topology
    positions = fixer.positions
    rebuilt_d = _d_residue_identities(topology, positions)
    rebuilt_fixed_d = rebuilt_d - generated_residue_set
    if rebuilt_fixed_d:
        raise BoundaryRefinementError(
            "temporary completion produced D-Cα geometry outside generated residues: "
            + _residue_labels(rebuilt_fixed_d)
        )
    generated_d = rebuilt_d & generated_residue_set
    repairs = fix_ca_chirality(topology, positions) if generated_d else 0
    if repairs != len(generated_d):
        raise BoundaryRefinementError("generated D-Cα repair did not repair every offender")
    assert_all_l(topology, positions)

    forcefield = ForceField(*BOUNDARY_REFINEMENT_FORCEFIELD)
    modeller = Modeller(topology, positions)
    modeller.addHydrogens(forcefield, pH=7.0)
    repair_misplaced_hydrogens(modeller.topology, modeller.positions)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
        rigidWater=False,
    )
    _constrain_generated_heavy_bonds(
        system,
        modeller.topology,
        modeller.positions,
        generated_residue_set,
    )
    for atom in modeller.topology.atoms():
        if _topology_residue(atom.residue) not in generated_residue_set:
            system.setParticleMass(atom.index, 0.0)
    chirality_force, protected_centres = build_ca_chirality_force(
        modeller.topology,
        modeller.positions,
    )
    if protected_centres:
        system.addForce(chirality_force)
    if restrain_amides:
        omega_force, restrained_amides = _build_trans_amide_force(
            modeller.topology,
            generated_residue_set,
        )
        if restrained_amides:
            system.addForce(omega_force)
    peptide_angle_force, restrained_angles = _build_peptide_angle_force(
        modeller.topology,
        generated_residue_set,
    )
    if restrained_angles:
        system.addForce(peptide_angle_force)

    try:
        platform = Platform.getPlatformByName(platform_name)
    except Exception as exc:
        raise BoundaryRefinementError(
            f"OpenMM platform {platform_name!r} is unavailable"
        ) from exc
    properties = {"Threads": "1"} if platform_name == "CPU" else {}
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
    simulation.minimizeEnergy(
        tolerance=BOUNDARY_REFINEMENT_TOLERANCE_KJ_MOL_NM
        * kilojoule_per_mole
        / nanometer,
        maxIterations=max_iterations,
    )
    best_state = simulation.context.getState(getEnergy=True, getPositions=True)
    best_score = _local_geometry_error_count(
        modeller.topology,
        best_state.getPositions(),
        generated_residue_set,
    )
    restart_positions = np.asarray(
        best_state.getPositions(asNumpy=True).value_in_unit(nanometer),
        dtype=np.float64,
    )
    mobile_atoms = [
        atom
        for atom in modeller.topology.atoms()
        if _topology_residue(atom.residue) in generated_residue_set
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
            score = _local_geometry_error_count(
                modeller.topology,
                state.getPositions(),
                generated_residue_set,
            )
            state_energy = float(
                state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            )
            best_energy = float(
                best_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            )
            if (score, state_energy) < (best_score, best_energy):
                best_state = state
                best_score = score
            if best_score == 0:
                break
    final_state = best_state
    refined_positions = final_state.getPositions()
    initial_energy = float(
        initial_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    )
    final_energy = float(
        final_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    )
    if not math.isfinite(initial_energy) or not math.isfinite(final_energy):
        raise BoundaryRefinementError("localized refinement produced non-finite energy")
    assert_all_l(modeller.topology, refined_positions)

    refined_coordinates = _coordinates_by_identity(
        modeller.topology,
        refined_positions,
    )
    for identity, original in original_coordinates.items():
        if _atom_residue(identity) in generated_residue_set:
            continue
        refined = refined_coordinates.get(identity)
        if refined is None or not np.array_equal(refined, original):
            raise BoundaryRefinementError(
                "localized refinement moved or removed a fixed atom: "
                f"{_atom_label(identity)}"
            )

    result_coordinates = {
        identity: refined_coordinates[identity]
        for identity in generated_atoms
        if identity in refined_coordinates
    }
    missing = generated_atom_set - set(result_coordinates)
    if missing:
        raise BoundaryRefinementError(
            "localized refinement lost generated atoms: " + _atom_labels(missing)
        )
    return BoundaryRefinementResult(
        coordinates_angstrom=result_coordinates,
        chirality_repairs=tuple(sorted(generated_d)),
        pre_coordinate_sha256=_coordinate_digest(
            {identity: original_coordinates[identity] for identity in generated_atoms}
        ),
        post_coordinate_sha256=_coordinate_digest(result_coordinates),
        initial_energy_kj_mol=initial_energy,
        final_energy_kj_mol=final_energy,
        platform=platform.getName(),
    )


def _coordinates_by_identity(topology: object, positions: object) -> dict[AtomIdentity, np.ndarray]:
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


def _build_trans_amide_force(
    topology: object,
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
    topology: object,
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
    system: object,
    topology: object,
    positions: object,
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
    topology: object,
    positions: object,
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


def _d_residue_identities(topology: object, positions: object) -> set[ResidueIdentity]:
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


def _topology_residue(residue: object) -> ResidueIdentity:
    return ResidueIdentity(
        residue.chain.id,
        str(residue.id),
        residue.insertionCode.strip(),
    )


def _atom_residue(atom: AtomIdentity) -> ResidueIdentity:
    return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)


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
