"""Private template baseline and held-out metrics, never a cc94 sampler trace.

Reference coordinates are used only for scoring, not proposal/refinement input.
"""

from __future__ import annotations

import inspect
import math
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from openmm import Platform
from openmm.unit import angstrom
from pdbfixer import PDBFixer

from dvbfixer.ffutils.geometry import assert_all_l, rebuild_missing_atoms_with_retry
from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
from dvbfixer.model.diffusion.atom_completion_materialize import materialize_atom_completion
from dvbfixer.model.diffusion.atom_completion_sampling import assess_completion_sampling
from dvbfixer.model.diffusion.atom_completion_validation import (
    _coordinates,
    _target_positions,
    audit_completion_candidate,
)
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionContractError, ResidueIdentity
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE


class TemplateBaselineUnavailable(RuntimeError):
    """An incomplete environment is not a rejected scientific proposal."""


@dataclass(frozen=True, slots=True)
class TemplateCompletionProposal:
    candidate: bytes
    temporary_extra_atoms: tuple[AtomIdentity, ...]
    temporary_fixed_max_displacement_angstrom: float
    platform: str
    revision: str = "seeded-pdbfixer-completion-baseline-v1"
    seed_policy: str = "shared-rebuild-helper-seeds-1-to-5; outer seed applies only to refinement"


@dataclass(frozen=True, slots=True)
class HeldOutCompletionMetrics:
    atom_count: int
    rmsd_angstrom: float
    maximum_displacement_angstrom: float
    revision: str = "exact-name-heldout-fixed-frame-v1"


def template_completion_proposal(
    request: AtomCompletionRequest, source: bytes, *, platform_name: str = "CPU",
) -> TemplateCompletionProposal:
    """Use the ordinary seeded heavy-atom helper, then retain only missing atoms.

    Whole-gap insertion is driven by the existing source-bound target placement,
    not a new sequence matcher, guessed numbering, or withheld coordinates.
    This baseline is explicitly distinct from checkpoint-generated evidence.
    """
    admission = assess_completion_sampling(request, source)
    if not admission.supported:
        raise DiffusionContractError("template comparison source failed the shared admission: " + ", ".join(admission.reasons))
    resource = Path(inspect.getfile(PDBFixer)).with_name("soft.xml")
    if not resource.is_file():
        raise TemplateBaselineUnavailable(f"PDBFixer resource is missing (not a scientific template rejection): {resource}")
    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    fixer = PDBFixer(pdbfile=StringIO(source.decode("utf-8")))
    fixer.platform = Platform.getPlatformByName(platform_name)
    chains = list(fixer.topology.chains())
    if len(chains) != 1 or chains[0].id != target.chain:
        raise DiffusionContractError("template baseline changed physical source chain identity")
    observed = tuple(ResidueIdentity(r.chain.id, r.id, r.insertionCode.strip()) for r in chains[0].residues())
    if observed != placement.observed_residues:
        raise DiffusionContractError("template baseline topology does not match observed target placement")
    fixer.findMissingResidues()
    fixer.missingResidues = {
        (chains[0].index, sum(index < gap.target_interval.start for index in placement.observed_target_indices)):
            [_CANONICAL_ONE_TO_THREE[target.sequence[index]]
             for index in range(gap.target_interval.start, gap.target_interval.stop)]
        for gap in request.gaps
    }
    fixer.findMissingAtoms()
    rebuild_missing_atoms_with_retry(fixer, verbose=True, log_prefix="template comparison: ")
    assert_all_l(fixer.topology, fixer.positions)
    rebuilt_residues = tuple(fixer.topology.residues())
    if (len(rebuilt_residues) != len(target.sequence)
            or tuple(r.name for r in rebuilt_residues) != tuple(_CANONICAL_ONE_TO_THREE[letter] for letter in target.sequence)):
        raise DiffusionContractError("template baseline reconstructed the wrong target sequence/length")
    mapping = _target_positions(request)
    positions = fixer.positions.value_in_unit(angstrom)
    rebuilt: dict[AtomIdentity, tuple[float, float, float]] = {}
    for index, residue in enumerate(rebuilt_residues):
        original = mapping[target.chain, index]
        for atom in residue.atoms():
            key = AtomIdentity(original.chain, original.residue_number, original.insertion_code, atom.name)
            if key in rebuilt:
                raise DiffusionContractError("template baseline produced duplicate mapped atom identities")
            point = positions[atom.index]
            rebuilt[key] = (float(point[0]), float(point[1]), float(point[2]))
    source_points = _coordinates(source.decode("utf-8").splitlines())
    drift = max(math.dist(rebuilt[atom], source_points[atom]) for atom in request.masks.fixed_atoms)
    if drift > 1e-4:
        raise DiffusionContractError("template baseline moved known source atoms in its temporary topology")
    expected = set(source_points) | set(request.generated_atoms)
    extras = tuple(atom for atom in rebuilt if atom not in expected)
    if any(atom.atom_name != "OXT" or atom.residue_number != mapping[target.chain, len(target.sequence) - 1].residue_number
           or atom.insertion_code != mapping[target.chain, len(target.sequence) - 1].insertion_code for atom in extras):
        raise DiffusionContractError("template baseline introduced unrequested nonterminal atoms")
    candidate = materialize_atom_completion(source, request, request.generated_atoms,
                                            [rebuilt[atom] for atom in request.generated_atoms])
    return TemplateCompletionProposal(candidate, extras, drift, platform_name)


def heldout_completion_metrics(
    request: AtomCompletionRequest, source: bytes, reference: bytes, candidate: bytes,
) -> HeldOutCompletionMetrics:
    """Exact-name RMSD in the unchanged source frame, with no fitting/leakage.

    The reference must restore exactly the masked observations and preserve every
    known record. Symmetry normalization is intentionally not claimed; callers
    must declare any chemically equivalent atom names when interpreting RMSD.
    """
    truth = audit_completion_candidate(source, reference, request)
    actual = audit_completion_candidate(source, candidate, request)
    distances = [math.dist(left, right) for left, right in zip(
        truth.generated_coordinates, actual.generated_coordinates, strict=True)]
    return HeldOutCompletionMetrics(len(distances), math.sqrt(sum(value * value for value in distances) / len(distances)),
                                    max(distances))
