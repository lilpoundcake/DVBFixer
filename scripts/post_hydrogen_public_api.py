"""Private transactional public-API adapter; never registered in production."""
from __future__ import annotations

import copy
from dataclasses import dataclass

from scripts.hydrogen_initialization_guard import GuardRefusal, Input, propose, validate
from scripts.qualify_hydrogen_initialization_guard import atoms_from_rows
from scripts.research_thr2_hb_only import transport

from dvbfixer.ffutils.geometry import HydrogenGeometryError
from dvbfixer.ffutils.variants import VARIANT_TO_PARENT


class PostHGuardError(HydrogenGeometryError):
    """Geometry/ownership error, never a template-mismatch fallback request."""

    def __init__(self, reason):
        self.findings = []
        RuntimeError.__init__(self, 'Private post-H guard refused: '+str(reason))


@dataclass(frozen=True)
class Outcome:
    returned: object
    before: Input
    raw: Input
    result: object | None
    admission_error: str | None
    maximum_api_observation_rounding_nm: float


def snapshot(model):
    from openmm.unit import nanometer

    rows = []
    for atom in model.topology.atoms():
        residue = atom.residue
        rows.append(dict(chain_index=residue.chain.index, residue_index=residue.index,
            key=[residue.chain.id, residue.id, residue.insertionCode.strip(), residue.name, atom.name],
            element=atom.element.symbol if atom.element is not None else 'unknown'))
    points = tuple(tuple(float(x) for x in p) for p in model.positions.value_in_unit(nanometer))
    value = Input(atoms_from_rows(rows), points,
                  tuple((a.index, b.index) for a, b in model.topology.bonds()), frozenset())
    validate(value)
    return value


def bond_metadata(model):
    return [(a.index, b.index, str(bond.type) if bond.type is not None else None, bond.order)
            for bond in model.topology.bonds() for a, b in [bond]]


def clone_model(model):
    """Rebind every Bond endpoint; naive deepcopy(Topology) does not do that."""
    from openmm.app import Modeller, Topology

    topology, atoms = Topology(), {}
    for chain in model.topology.chains():
        copied_chain = topology.addChain(chain.id)
        for residue in chain.residues():
            copied_residue = topology.addResidue(residue.name, copied_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                atoms[atom] = topology.addAtom(atom.name, atom.element, copied_residue, atom.id,
                                               getattr(atom, 'formalCharge', None))
    for bond in model.topology.bonds():
        topology.addBond(atoms[bond[0]], atoms[bond[1]], bond.type, bond.order)
    box = model.topology.getPeriodicBoxVectors()
    if box is not None:
        topology.setPeriodicBoxVectors(copy.deepcopy(box))
    return Modeller(topology, copy.deepcopy(model.positions))


def observed_bond_metadata(before, raw, before_types, raw_types):
    old = {a.identity for a in before.atoms}

    def keyed(value, records):
        return {tuple(sorted((value.atoms[a].identity, value.atoms[b].identity))):(kind, order)
                for a, b, kind, order in records
                if value.atoms[a].identity in old and value.atoms[b].identity in old}

    if keyed(before, before_types) != keyed(raw, raw_types):
        raise GuardRefusal('API changed existing bond metadata')


def ownership(before, raw, protected_hydrogens, reference=False):
    """No System/private frame introspection: exact identity/graph difference."""
    old = {a.identity:i for i, a in enumerate(before.atoms)}
    new = {a.identity:i for i, a in enumerate(raw.atoms)}
    if not isinstance(protected_hydrogens, frozenset):
        raise GuardRefusal('explicit immutable original-H lineage required')
    if not protected_hydrogens <= old.keys() or not protected_hydrogens <= new.keys():
        raise GuardRefusal('original observed H absent before/after API; do not relabel stripped H new')
    if any(before.atoms[old[k]].element not in ('H', 'D') for k in protected_hydrogens):
        raise GuardRefusal('original-H lineage contains a heavy atom')
    if not old.keys() <= new.keys():
        raise GuardRefusal('existing atom deleted by API')
    for identity in old:
        a, b = before.atoms[old[identity]], raw.atoms[new[identity]]
        if (a.element != b.element or VARIANT_TO_PARENT.get(a.residue_name, a.residue_name)
                != VARIANT_TO_PARENT.get(b.residue_name, b.residue_name)):
            raise GuardRefusal('existing atom chemistry changed')
    created = frozenset(new.keys()-old.keys())
    if any(raw.atoms[new[k]].element != 'H' for k in created):
        raise GuardRefusal('API created a non-H atom')

    def old_graph(value):
        return {tuple(sorted((value.atoms[a].identity, value.atoms[b].identity)))
                for a, b in value.bonds
                if value.atoms[a].identity in old and value.atoms[b].identity in old}

    if old_graph(before) != old_graph(raw):
        raise GuardRefusal('API changed an existing-atom bond graph')
    error = transport([list(p) for p in before.positions_nm],
                      [list(raw.positions_nm[new[a.identity]]) for a in before.atoms],
                      'Reference' if reference else 'captured')
    value = Input(raw.atoms, raw.positions_nm, raw.bonds, created)
    validate(value)
    return value, error


def apply(model, original_call, *args, protected_hydrogens, enabled=True, journal=None, **kwargs):
    """Delegate one original call on a clone; commit only an admitted result.

    Baseline mode commits the raw API return even if observational guard
    admission fails; enabled mode refuses atomically. API errors always propagate.
    """
    from openmm import Vec3
    from openmm.unit import nanometer

    if kwargs.get('residueTemplates') or (len(args) >= 5 and args[4]):
        raise PostHGuardError('nonempty object-keyed residueTemplates are outside private scope')
    try:
        before = snapshot(model)
    except (GuardRefusal, ValueError) as error:
        raise PostHGuardError(error) from error
    # Validate caller lineage before creating H, including upstream stripped H.
    before_ids = {a.identity:a for a in before.atoms}
    if (not isinstance(protected_hydrogens, frozenset)
            or not protected_hydrogens <= before_ids.keys()
            or any(before_ids[k].element not in ('H', 'D') for k in protected_hydrogens)):
        raise PostHGuardError('missing/invalid original observed-H lineage before API')
    working = clone_model(model)
    before_types = bond_metadata(model)
    try:
        returned = original_call(working, *args, **kwargs)
    except BaseException as error:
        if journal is not None:
            try:
                failed_raw, snapshot_error = snapshot(working), None
            except (GuardRefusal, ValueError) as snapshot_failure:
                failed_raw, snapshot_error = None, str(snapshot_failure)
            journal(dict(stage='api-error', before=before, raw=failed_raw,
                         snapshot_error=snapshot_error, error_type=type(error).__name__, error=str(error)))
        raise
    try:
        raw = snapshot(working)
    except (GuardRefusal, ValueError) as error:
        if journal is not None:
            journal(dict(stage='guard-refusal', before=before, raw=None, admission_error=str(error)))
        raise PostHGuardError(error) from error
    result, reason, api_error = None, None, 0.
    raw_types = bond_metadata(working)
    platform = kwargs.get('platform', args[3] if len(args) >= 4 else None)
    reference = platform is not None and platform.getName() == 'Reference'
    try:
        owned, api_error = ownership(before, raw, protected_hydrogens, reference)
        observed_bond_metadata(before, raw, before_types, raw_types)
        result = propose(owned)
    except (GuardRefusal, ValueError) as error:
        reason = str(error)
    if journal is not None:
        journal(dict(stage='guard-refusal' if reason else 'guard-proposal', before=before, raw=raw,
                     result=result, admission_error=reason, api_observation_rounding_nm=api_error,
                     before_bond_metadata=before_types, raw_bond_metadata=raw_types,
                     reference_requested=reference, enabled=enabled, returned=returned))
    if reason and enabled:
        raise PostHGuardError(reason)
    if enabled:
        working.positions = [Vec3(*p) for p in result.positions_nm]*nanometer
    # Commit last: failures in the call, geometry or diagnostic sink leave model intact.
    model.topology, model.positions = working.topology, working.positions
    return Outcome(returned, before, raw, result, reason, api_error)
