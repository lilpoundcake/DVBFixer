"""Research-only seeded PDBFixer proposal with retained per-seed temporary states.

No held-out coordinates, chirality reflection, H placement or production routing.
Reference is a declared numerical control, never a new production platform default.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import inspect
import math
from collections.abc import Callable
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.contracts import Addition, CompletionRequest


class BaselineUnavailable(RuntimeError):
    """Missing/broken resources are environmental, not scientific rejection."""


@dataclass(frozen=True)
class TemplateSnapshot:
    seed: int
    atoms: tuple[tuple[ExactAtomRef, str, tuple[float, float, float]], ...]
    observed_max_movement_angstrom: float


@dataclass(frozen=True)
class TemplateProposal:
    additions: tuple[Addition, ...]
    temporary_extra_atoms: tuple[ExactAtomRef, ...]
    snapshots: tuple[TemplateSnapshot, ...]
    platform: str = "Reference"
    seed_policy: str = "shared-helper-1-to-5; no reflection; retain every returned state"


def preflight() -> dict:
    from openmm import Platform
    from openmm.app import ForceField
    from pdbfixer import PDBFixer

    root = Path(inspect.getfile(PDBFixer)).parent
    required = [root/"soft.xml", *(root/"templates"/f"{name}.pdb" for name in ("VAL", "THR"))]
    # Upstream names are lower case in some installations; resolve explicitly.
    required = [p if p.is_file() else p.with_name(p.name.lower()) for p in required]
    if any(not p.is_file() for p in required):
        raise BaselineUnavailable("PDBFixer soft.xml/VAL/THR templates unavailable")
    try:
        ForceField(str(root/"soft.xml"))
        Platform.getPlatformByName("Reference")
    except Exception as exc:
        raise BaselineUnavailable(f"PDBFixer resources/Reference unusable: {exc}") from exc
    return {"pdbfixer": importlib.metadata.version("pdbfixer"), "openmm": importlib.metadata.version("openmm"),
            "platform": "Reference", "source_sha256": hashlib.sha256(Path(inspect.getfile(PDBFixer)).read_bytes()).hexdigest(),
            "resource_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [root/"soft.xml", *sorted((root/"templates").glob("*.pdb"))]},
            "helper_sha256": hashlib.sha256(Path(inspect.getfile(_helper())).read_bytes()).hexdigest()}


def _helper():
    from dvbfixer.ffutils.geometry import rebuild_missing_atoms_with_retry
    return rebuild_missing_atoms_with_retry


def template_proposal(source: bytes, request: CompletionRequest, *,
                       retain_snapshot: Callable[[TemplateSnapshot], None] | None = None) -> TemplateProposal:
    from openmm import Platform
    from openmm.unit import angstrom
    from pdbfixer import PDBFixer

    preflight()
    if hashlib.sha256(source).hexdigest() != request.source_sha256:
        raise ValueError("template source/request SHA mismatch")
    if any(o.atom.component.model != 1 or o.atom.component.occurrence != 1 or o.atom.component.alternate_location
           for o in request.observations):
        raise ValueError("template baseline requires single-frame unique observations")
    fixer = PDBFixer(pdbfile=StringIO(source.decode("ascii")))
    fixer.platform = Platform.getPlatformByName("Reference")
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # No whole-gap invention in a partial-sidechain comparator.
    fixer.findMissingAtoms()
    observed = {o.atom: o for o in request.observations}
    expected = {ExactAtomRef(request.residue, name) for name in request.missing_atoms}
    snapshots = []

    def snapshot(seed):
        positions = fixer.positions.value_in_unit(angstrom)
        atoms = []
        for a in fixer.topology.atoms():
            r = a.residue
            identity = ExactAtomRef(ComponentInstanceRef(1, r.chain.id, r.id, r.insertionCode.strip()), a.name)
            atoms.append((identity, a.element.symbol.upper(), tuple(float(x) for x in positions[a.index])))
        by_atom = {a[0]: a for a in atoms}
        if len(by_atom) != len(atoms) or not observed.keys() <= by_atom.keys():
            raise ValueError("template temporary topology lost/duplicated observed identities")
        movement = max(math.dist(o.xyz, by_atom[atom][2]) for atom, o in observed.items())
        snapshots.append(TemplateSnapshot(seed, tuple(atoms), movement))
        if retain_snapshot is not None:
            retain_snapshot(snapshots[-1])

    original = fixer.addMissingAtoms

    def traced_rebuild(seed=None):
        # Instance-local observation hook. Production helper owns all calls/seeds;
        # never invent an extra direct rebuild or retry after its result.
        original(seed=seed)
        snapshot(seed)

    fixer.addMissingAtoms = traced_rebuild
    try:
        _helper()(fixer, verbose=True, log_prefix="branched comparator: ")
    finally:
        fixer.addMissingAtoms = original
    if not snapshots:
        raise ValueError("template helper returned without a retained rebuild state")
    last = snapshots[-1]
    if last.observed_max_movement_angstrom > 1e-6:
        raise ValueError("template temporary state moved observed atoms")
    by_atom = {a[0]: a for a in last.atoms}
    if not expected <= by_atom.keys():
        raise ValueError("template proposal omitted declared missing atoms")
    extras = tuple(a for a in by_atom if a not in observed and a not in expected)
    if any(a.name != "OXT" for a in extras):
        raise ValueError("template rebuilt unrequested nonterminal heavy atoms")
    additions = tuple(Addition(a, by_atom[a][1], by_atom[a][2]) for a in sorted(expected))
    return TemplateProposal(additions, extras, tuple(snapshots))
