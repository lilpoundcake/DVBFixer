"""Private cc94 completion scope and lossless source-to-engine transport.

Not registered in the public profile catalog. Wire support is not qualification.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from io import StringIO

from openmm.app import PDBFile
from openmm.unit import angstrom

from dvbfixer.diagnose.report import Finding, Severity
from dvbfixer.model.diffusion.atom_completion import CompletionKind
from dvbfixer.model.diffusion.atom_completion_contract import (
    AtomCompletionCapabilities,
    AtomCompletionRequest,
    audit_completion_request_source,
    require_completion_capability,
)
from dvbfixer.model.diffusion.atom_completion_materialize import _source_layout
from dvbfixer.model.diffusion.atom_completion_validation import (
    _geometry_findings,
    _target_positions,
)
from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    DiffusionContractError,
    GapKind,
    ResidueIdentity,
)
from dvbfixer.model.diffusion.protein_atoms import (
    _CANONICAL_ONE_TO_THREE,
    canonical_heavy_atom_identities,
)

COMPLETION_RESEARCH_PROFILE = "cc94-atom-completion-research-v1"
COMPLETION_TRANSPORT_REVISION = "cc94-completion-full-monomer-v1"


@dataclass(frozen=True, slots=True)
class CompletionSamplingAdmission:
    supported: bool
    reasons: tuple[str, ...]
    source_findings: tuple[Finding, ...]
    scoped_findings: tuple[Finding, ...]
    excluded_gap_bridges: tuple[tuple[AtomIdentity, AtomIdentity], ...]


@dataclass(frozen=True, slots=True)
class CompletionEngineView:
    target_residues: tuple[ResidueIdentity, ...]
    atom_axis: tuple[AtomIdentity, ...]
    atom_target_indices: tuple[int, ...]
    represented_fixed_atoms: tuple[AtomIdentity, ...]
    preserved_only_atoms: tuple[AtomIdentity, ...]
    motif_text: str
    motif_placement: str
    fixed_positions: tuple[int, ...]
    residue_index: tuple[int, ...]
    chain_index: tuple[int, ...]
    revision: str = COMPLETION_TRANSPORT_REVISION


def assess_completion_sampling(request: AtomCompletionRequest, source: bytes) -> CompletionSamplingAdmission:
    """Narrow research admission before creating motif files or loading weights.

    Only one <=128-residue heavy-atom protein chain, <=4 partial sidechains,
    and optionally one internal 3-12-residue gap. Existing bad geometry is not
    repaired. Explicit bonds, heterogens and backbone defects remain unsupported.
    """
    audit_completion_request_source(request, source)
    lines, _, _ = _source_layout(source.decode("utf-8"))
    reasons: list[str] = []
    if request.profile != COMPLETION_RESEARCH_PROFILE:
        reasons.append("not-the-private-completion-profile")
    if len(request.target_sequences) != 1 or sum(len(t.sequence) for t in request.target_sequences) > 128:
        reasons.append("research-context-requires-one-chain-at-most-128-residues")
    if len(request.seeds) != 1:
        reasons.append("research-invocation-requires-one-predeclared-seed")
    if len(request.masks.partial_residues) > 4 or len(request.generated_atoms) > 256:
        reasons.append("research-completion-resource-limit")
    if any(partial.kind is CompletionKind.BACKBONE for partial in request.masks.partial_residues):
        reasons.append("partial-backbone-not-enabled")
    if len(request.gaps) > 1 or any(gap.gap_kind is not GapKind.INTERNAL
                                  or not 3 <= len(gap.generated_residues) <= 12 for gap in request.gaps):
        reasons.append("research-gap-requires-one-internal-region-of-3-to-12-residues")
    if any(line[:6].strip() == "HETATM" for line in lines):
        reasons.append("heterogen-context-not-enabled")
    if any(line[:6].strip() in {"LINK", "SSBOND", "CONECT"} for line in lines):
        reasons.append("explicit-covalent-context-not-enabled")
    for line in lines:
        if line[:6].strip() != "ATOM":
            continue
        if line[76:78].strip().upper() in {"H", "D"} or line[12:16].strip().lstrip("0123456789").startswith("H"):
            reasons.append("research-source-must-be-heavy-atom-only")
        try:
            if float(line[54:60]) != 1.0:
                reasons.append("partial-occupancy-not-enabled")
        except ValueError:
            reasons.append("invalid-source-occupancy")
    if reasons:
        return CompletionSamplingAdmission(False, tuple(dict.fromkeys(reasons)), (), (), ())
    require_completion_capability(request, AtomCompletionCapabilities(
        1, COMPLETION_RESEARCH_PROFILE, sidechain=True, whole_gaps=True, mixed_gaps=True,
        max_target_residues=128, max_generated_atoms=256, max_candidates=1,
    ))
    pdb = PDBFile(StringIO(source.decode("utf-8")))
    raw_findings = _geometry_findings(pdb)
    coordinates = {AtomIdentity(atom.residue.chain.id, atom.residue.id, atom.residue.insertionCode.strip(), atom.name):
                   pdb.positions[atom.index].value_in_unit(angstrom) for atom in pdb.topology.atoms()}
    bridges: list[tuple[AtomIdentity, AtomIdentity]] = []
    for gap in request.gaps:
        assert gap.left_anchor is not None and gap.right_anchor is not None
        c = AtomIdentity(gap.chain, gap.left_anchor.residue_number, gap.left_anchor.insertion_code, "C")
        n = AtomIdentity(gap.chain, gap.right_anchor.residue_number, gap.right_anchor.insertion_code, "N")
        if math.dist(coordinates[c], coordinates[n]) <= 1.8:
            reasons.append("gap-crosses-observed-peptide-connectivity")
        else:
            bridges.append((c, n))
    # PDBFile may infer a peptide across missing target residues. Exclude ONLY
    # those declared, spatially disconnected C-N bridges from this temporary
    # diagnostic topology. Keep raw findings and the exact exclusions in evidence.
    excluded = {frozenset(pair) for pair in bridges}
    pdb.topology._bonds = [bond for bond in pdb.topology.bonds() if frozenset(
        AtomIdentity(atom.residue.chain.id, atom.residue.id, atom.residue.insertionCode.strip(), atom.name)
        for atom in bond) not in excluded]
    scoped = _geometry_findings(pdb)
    if any(finding.severity is Severity.ERROR or finding.category == "non_planar_amide" for finding in scoped):
        reasons.append("source-geometry-errors")
    if any("truncated" in finding.message.lower() for finding in scoped):
        reasons.append("source-diagnostics-incomplete")
    positions = _target_positions(request)
    observed = set(request.sequence_placements[0].observed_residues)
    for index in range(len(request.target_sequences[0].sequence) - 1):
        left, right = positions[request.target_sequences[0].chain, index], positions[request.target_sequences[0].chain, index + 1]
        if left in observed and right in observed:
            c = AtomIdentity(left.chain, left.residue_number, left.insertion_code, "C")
            n = AtomIdentity(right.chain, right.residue_number, right.insertion_code, "N")
            if not 1.20 <= math.dist(coordinates[c], coordinates[n]) <= 1.45:
                reasons.append("observed-peptide-geometry-error")
    return CompletionSamplingAdmission(not reasons, tuple(dict.fromkeys(reasons)), raw_findings, scoped, tuple(bridges))


def completion_engine_view(request: AtomCompletionRequest, source: bytes) -> CompletionEngineView:
    """Full monomer transport; planning alone grants no execution admission."""
    audit_completion_request_source(request, source)
    if len(request.target_sequences) != 1:
        raise DiffusionContractError("completion engine view currently requires one physical source chain")
    positions = _target_positions(request)
    target = request.target_sequences[0]
    residues = tuple(positions[target.chain, index] for index in range(len(target.sequence)))
    atom_axis: list[AtomIdentity] = []
    indices: list[int] = []
    for index, residue in enumerate(residues):
        atoms = canonical_heavy_atom_identities(residue, _CANONICAL_ONE_TO_THREE[target.sequence[index]])
        atom_axis.extend(atoms)
        indices.extend([index] * len(atoms))
    canonical = set(atom_axis)
    fixed = set(request.masks.fixed_atoms) & canonical
    preserved = tuple(atom for atom in request.masks.fixed_atoms if atom not in canonical)
    observed_indices = dict(zip(request.sequence_placements[0].observed_residues,
                                request.sequence_placements[0].observed_target_indices, strict=True))
    motif: list[str] = []
    seen: set[AtomIdentity] = set()
    for line in source.decode("utf-8").splitlines():
        if line[:6].strip() != "ATOM":
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom not in fixed:
            continue
        residue = ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)
        ordinal = observed_indices[residue] + 1
        motif.append(f"ATOM  {len(motif) + 1:5d}" + line[11:21] + f"A{ordinal:4d} " + line[27:] + "\n")
        seen.add(atom)
    if seen != fixed:
        raise DiffusionContractError("completion motif does not cover represented fixed atoms exactly")
    motif.append("TER\nEND\n")
    segments: list[str] = []
    cursor = 0
    for gap in sorted(request.gaps, key=lambda item: item.target_interval.start):
        if cursor < gap.target_interval.start:
            segments.append(f"A{cursor + 1}-{gap.target_interval.start}")
        segments.append(str(gap.target_interval.stop - gap.target_interval.start))
        cursor = gap.target_interval.stop
    if cursor < len(target.sequence):
        segments.append(f"A{cursor + 1}-{len(target.sequence)}")
    return CompletionEngineView(residues, tuple(atom_axis), tuple(indices),
                                 tuple(atom for atom in request.masks.fixed_atoms if atom in fixed), preserved,
                                 "".join(motif), "/".join(segments),
                                 request.sequence_placements[0].observed_target_indices,
                                 tuple(range(1, len(target.sequence) + 1)), (0,) * len(target.sequence))
