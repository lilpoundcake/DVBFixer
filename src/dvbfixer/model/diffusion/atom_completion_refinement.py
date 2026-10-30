"""Generated-atom-only refinement adapter for private completion controls.

No observed mobility, runner enablement or checkpoint qualification is granted.
"""

from __future__ import annotations

from dataclasses import dataclass
from io import StringIO

from openmm import OpenMMException
from openmm.app import PDBFile

from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
from dvbfixer.model.diffusion.atom_completion_materialize import materialize_atom_completion
from dvbfixer.model.diffusion.atom_completion_validation import (
    CompletionValidation,
    audit_completion_candidate,
    validate_completion_candidate,
)
from dvbfixer.model.diffusion.boundary_refinement import (
    BOUNDARY_REFINEMENT_MAX_ITERATIONS,
    BOUNDARY_REFINEMENT_RESTART_COUNT,
    BOUNDARY_REFINEMENT_REVISION,
    BoundaryRefinementError,
    refine_generated_region,
)
from dvbfixer.model.diffusion.contract import DiffusionContractError, ResidueIdentity

COMPLETION_REFINEMENT_REVISION = "atom-completion-generated-only-refinement-v1"


class CompletionRefinementError(BoundaryRefinementError):
    """A failed attempt retains its raw diagnostics and untouched raw bytes."""

    def __init__(self, message: str, raw_candidate: bytes, raw_validation: CompletionValidation) -> None:
        super().__init__(message)
        self.raw_candidate = raw_candidate
        self.raw_validation = raw_validation


@dataclass(frozen=True, slots=True)
class CompletionRefinement:
    raw_candidate: bytes
    refined_candidate: bytes
    raw_validation: CompletionValidation
    refined_validation: CompletionValidation
    initial_energy_kj_mol: float
    final_energy_kj_mol: float
    platform: str
    seed: int
    kernel_revision: str = BOUNDARY_REFINEMENT_REVISION
    revision: str = COMPLETION_REFINEMENT_REVISION


def refine_completion_candidate(
    source: bytes,
    candidate: bytes,
    request: AtomCompletionRequest,
    *,
    seed: int,
    max_iterations: int = BOUNDARY_REFINEMENT_MAX_ITERATIONS,
    restart_count: int = BOUNDARY_REFINEMENT_RESTART_COUNT,
    platform_name: str = "Reference",
) -> CompletionRefinement:
    """Refine exactly the requested missing atoms, never a whole observed residue.

    Raw and refined validation are both retained, including rejected refinement
    outcomes. Malformed identity/ownership is refused before the force field runs.
    This initial adapter is protein-only and refuses raw D-Cα rather than risk
    reflecting observed atoms of a mixed observed/generated side chain.
    """
    audit_completion_candidate(source, candidate, request)
    if type(seed) is not int or seed not in request.seeds:
        raise DiffusionContractError("completion refinement seed was not declared by the request")
    if any(line[:6].strip() == b"HETATM" for line in candidate.splitlines()):
        raise DiffusionContractError("private completion refinement currently supports protein-only inputs")
    raw_validation = validate_completion_candidate(source, candidate, request)
    pdb = PDBFile(StringIO(candidate.decode("utf-8")))
    try:
        assert_all_l(pdb.topology, pdb.positions)
    except ChiralityError as exc:
        raise CompletionRefinementError(
            "completion refinement refuses raw D-Cα; observed sidechain reflection is forbidden",
            candidate, raw_validation,
        ) from exc
    residues = tuple(dict.fromkeys(ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)
                                   for atom in request.generated_atoms))
    try:
        refined = refine_generated_region(
            candidate.decode("utf-8"), generated_residues=residues, generated_atoms=request.generated_atoms,
            max_iterations=max_iterations, restart_count=restart_count, platform_name=platform_name, random_seed=seed,
        )
    except (BoundaryRefinementError, ChiralityError, OpenMMException, ValueError) as exc:
        raise CompletionRefinementError(f"completion refinement failed: {exc}", candidate, raw_validation) from exc
    if set(refined.coordinates_angstrom) != set(request.generated_atoms):
        raise CompletionRefinementError("completion refinement returned coordinates outside the exact generated mask",
                                         candidate, raw_validation)
    output = materialize_atom_completion(source, request, request.generated_atoms,
                                         [tuple(float(value) for value in refined.coordinates_angstrom[atom])
                                          for atom in request.generated_atoms])
    final_validation = validate_completion_candidate(source, output, request)
    return CompletionRefinement(candidate, output, raw_validation, final_validation,
                                 refined.initial_energy_kj_mol, refined.final_energy_kj_mol, refined.platform, seed)
