"""Frozen all20 source/state coverage; heavy holdouts never drive construction."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import replace
from pathlib import Path

from scripts.chemical_reconstruction_ccd_controls import FIXTURES as ORIGINAL_CCD
from scripts.chemical_reconstruction_ccd_controls import ROOT, _write
from scripts.chemical_reconstruction_declared_state_controls import (
    check_pins as check_declared_pins,
)
from scripts.chemical_reconstruction_extended_h_controls import FIXTURES as GROUP_CCD
from scripts.chemical_reconstruction_hydrogen_controls import (
    FFBONDED,
    FIXED_POLICY,
    GENERATED_POLICY,
    _execute,
    backbone_types,
)
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import POLICY as HEAVY_POLICY
from scripts.chemical_reconstruction_peptide_controls import load_ccd
from scripts.chemical_reconstruction_source_controls import _bytes

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.ccd_authority import load_geometry_authority
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Connection,
    LinkKind,
    Observation,
    Operation,
    Origin,
    Ownership,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.declared_peptide import (
    PeptideStateDeclaration,
    compile_declared_hydrogens,
    derive_declared_peptide,
)
from dvbfixer.chemical_reconstruction.heldout_scoring import (
    HoldoutScoringError,
    score_heavy_holdout,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.measurement import _dot, _sub, _unit
from dvbfixer.chemical_reconstruction.observed_peptide import (
    _parse_source,
    compile_peptide_request_plan,
    select_interior_target,
)
from dvbfixer.chemical_reconstruction.peptide_state import PeptideSubstitution
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.nonprotein_reconstruction import ExternalLink
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref

PAIRS = (('ALA', 'ALA'), ('ARG', 'ARG'), ('ASN', 'ASN'), ('ASP', 'ASH'), ('CYS', 'CYS'),
         ('GLN', 'GLN'), ('GLU', 'GLH'), ('GLY', 'GLY'), ('HIS', 'HIP'), ('ILE', 'ILE'),
         ('LEU', 'LEU'), ('LYS', 'LYS'), ('MET', 'MET'), ('PHE', 'PHE'), ('PRO', 'PRO'),
         ('SER', 'SER'), ('THR', 'THR'), ('TRP', 'TRP'), ('TYR', 'TYR'), ('VAL', 'VAL'))
PROTOCOL = ROOT / 'docs/plans/chemical-reconstruction-all20-v1.md'
PROTOCOL_SHA256 = '0bafbf6d15868455c2cfe4ae822e7193c90ba09a7f6dbf18c3be8af6d53bdcb2'
FIXTURES = ROOT / 'tests/fixtures/chemical_reconstruction/ccd_all20_v1'
MANIFEST_SHA256 = 'acf270e09d298ca12c19c79013a54a85b5b845f602e863f91d6c886191e95f04'


def check_pins():
    check_declared_pins()
    for path, pin in ((PROTOCOL, PROTOCOL_SHA256), (FIXTURES/'manifest.json', MANIFEST_SHA256)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
            raise ValueError(f'frozen all20 pin mismatch: {path}')
    for directory in (ORIGINAL_CCD, GROUP_CCD, FIXTURES):
        for name, pin in json.loads((directory/'manifest.json').read_text())['components'].items():
            if hashlib.sha256((directory/f'{name}.cif').read_bytes()).hexdigest() != pin['sha256']:
                raise ValueError(f'frozen CCD source mismatch: {name}')


def load_parent(name):
    found = []
    for directory in (ORIGINAL_CCD, GROUP_CCD, FIXTURES):
        pins = json.loads((directory/'manifest.json').read_text())['components']
        if name in pins:
            found.append((directory/f'{name}.cif', pins[name]['sha256']))
    if len(found) != 1:
        raise ValueError('all20-needs-one-exact-pinned-source')
    return load_geometry_authority(found[0][0], name, expected_sha256=found[0][1])


def authority_for(source, state, identity, evidence):
    parent, alanine = load_parent(source), load_ccd('ALA')
    inventory = transport_inventory(raw_inventory(state), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
    left = ComponentInstanceRef(identity.model, '__left_chemical_witness__', '1')
    right = ComponentInstanceRef(identity.model, '__right_chemical_witness__', '1')
    leaving = 'H' if source == 'PRO' else 'H2'
    substitutions = tuple(PeptideSubstitution(ExactAtomRef(a, 'C'), ExactAtomRef(b, 'N'),
                          frozenset((ExactAtomRef(a, 'OXT'), ExactAtomRef(a, 'HXT'), ExactAtomRef(b, h))))
                          for a, b, h in ((left, identity, leaving), (identity, right, 'H2')))
    declaration = PeptideStateDeclaration(identity, state, parent.digest, alanine.digest, digest_record(inventory.canonical), substitutions)
    evidence.update(parent_authority=parent, alanine_authority=alanine, inventory_transport=inventory, declaration=declaration)
    authority = derive_declared_peptide(parent, alanine, identity, inventory.canonical, declaration)
    evidence['authority'] = authority
    return authority


def sidechain_heavy_mask(authority):
    return frozenset(a.name for a in authority.state.graph.atoms if a.element not in {'H', 'D'} and a.name not in {'N', 'CA', 'C', 'O'})


def synthetic_input(source, state, evidence):
    identity = ComponentInstanceRef(1, 'a', '82', 'B', occurrence=2)
    authority = authority_for(source, state, identity, evidence)
    graph = authority.state.graph
    points = {a.name: a.ideal_position for a in graph.atoms if a.element not in {'H', 'D'}}
    n = points['N']
    if source == 'PRO':
        neighbours = sorted(b.atom2 if b.atom1 == 'N' else b.atom1 for b in graph.bonds if 'N' in {b.atom1, b.atom2})
        a, b = (_unit(_sub(points[name], n)) for name in neighbours)
        direction = _unit(tuple(-a[i]-b[i] for i in range(3)))
    else:
        u, c = _unit(_sub(points['CA'], n)), _unit(_sub(points['C'], n))
        projection = _dot(c, u)
        v = _unit(tuple(c[i]-projection*u[i] for i in range(3)))
        angle = math.radians(121.96)
        direction = tuple(math.cos(angle)*u[i]+math.sin(angle)*v[i] for i in range(3))
    incoming = ExactAtomRef(ComponentInstanceRef(1, 'A', '81'), 'C')
    observations = tuple(Observation(ExactAtomRef(identity, a.name), a.element, points[a.name], Unit.ANGSTROM)
                         for a in graph.atoms if a.element not in {'H', 'D'}) + (
                             Observation(incoming, 'C', tuple(n[i]+1.335*direction[i] for i in range(3)), Unit.ANGSTROM),)
    pin = digest_record(('synthetic-all20-heavy-input-v1', source, state, observations, authority.declaration.digest))
    request = ReconstructionRequest(pin, Operation.HEAVY,
            (Component(identity, ChemistryKind.PROTEIN, (authority.state,)), Component(incoming.component, ChemistryKind.UNKNOWN)),
            observations, frozenset({identity}), Ownership(frozenset()),
            (Connection(ExternalLink(incoming, ExactAtomRef(identity, 'N'), 1.), LinkKind.COVALENT),),
            authority.constraints(minimum_volume=.05, planar_tolerance_angstrom=.05),
            origin=Origin('synthetic-all20-heavy-input-v1', pin))
    return request, authority, incoming


def observed_input(source, state, operation, directory, evidence):
    path = 'tests/fixtures/regressions/1FR2.pdb' if source in {'CYS', 'TRP'} else 'tests/fixtures/protein_context/1UBQ.pdb'
    original = (ROOT/path).read_bytes()
    _bytes(directory/'original.pdb', original)
    preparation = prepare_source(original)
    _bytes(directory/'prepared.pdb', preparation.prepared_text.encode('ascii'))
    evidence.update(source=path, original_sha256=preparation.original_sha256,
                    preparation={'digest': preparation.digest, 'prepared_sha256': preparation.prepared_sha256,
                                 'removed_metadata': preparation.removed_metadata})
    identity = select_interior_target(original, source, chain='A', number=None)
    authority = authority_for(source, state, identity, evidence)
    mask = sidechain_heavy_mask(authority) if operation == 'heavy' else frozenset()
    probe = prepare_source_probe(preparation, authority, withheld=mask)
    source_audit = audit_source_transport(preparation, probe)
    evidence.update(probe=probe, source_transport_audit=source_audit)
    if not source_audit.passed:
        raise ValueError('all20-source-transport-failed')
    return probe.request, authority, probe.real_peptide_partners[0]


def compile_H(base, authority, incoming, evidence):
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    types = None if authority.inventory.state_name == 'PRO' else backbone_types(authority.inventory.state_name)
    evidence.update(parameters=parameters, explicit_backbone_types=types)
    return compile_declared_hydrogens(base, authority, incoming, parameters, backbone_types=types,
                                     fixed_policy=FIXED_POLICY, generated_policy=GENERATED_POLICY)


def heldout_after_proposal(arm, authority, directory):
    if arm == 'synthetic':
        # Same-authority ideal consistency, explicitly NOT a native holdout.
        wanted = sidechain_heavy_mask(authority)
        return tuple(Observation(ExactAtomRef(authority.identity, a.name), a.element, a.ideal_position, Unit.ANGSTROM)
                     for a in authority.state.graph.atoms if a.name in wanted)
    _, atoms = _parse_source((directory/'original.pdb').read_bytes())
    wanted = sidechain_heavy_mask(authority)
    return tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position, Unit.ANGSTROM)
                 for a in atoms if _component_ref(a) == authority.identity and a.name in wanted)


def execute_heavy(request, authority, arm, directory, item, evidence):
    mask = sidechain_heavy_mask(authority)
    if not mask:
        item.update(status='not-applicable-no-sidechain-heavy-atoms', reason='GLY-empty-sidechain-mask', proposal_entry_called=False)
        return
    plan = compile_peptide_request_plan(request, HEAVY_POLICY)
    evidence['plan'] = plan
    item['proposal_entry_called'] = True
    result = propose(request, plan)
    evidence['result'] = result
    item.update(status=result.status, reason=result.reason, coordinates_produced=result.coordinates_produced,
                coordinate_constructor_executed=result.raw is not None)
    # Do not extract/reference held-out observations until AFTER proposal returns.
    reference = heldout_after_proposal(arm, authority, directory)
    evidence['heldout_reference_after_proposal'] = reference
    evidence['heavy_scores'] = {}
    for kind, candidate in (('raw', result.raw), ('rounded', result.rounded)):
        if candidate is None:
            evidence['heavy_scores'][kind] = {'status': 'not-scored-no-candidate'}
            continue
        try:
            score = score_heavy_holdout(request, plan, candidate, reference)
        except HoldoutScoringError as exc:
            evidence['heavy_scores'][kind] = {'status': 'not-scored', 'reason': str(exc)}
        else:
            evidence['heavy_scores'][kind] = {'status': 'measured-not-an-acceptance-gate', 'score': score}


def run_controls(output: Path):
    check_pins()
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for arm in ('synthetic', 'observed'):
        for operation in ('heavy', 'H'):
            for source, state in PAIRS:
                directory = output/f'{arm}-{operation}-{source}-{state}'
                directory.mkdir()
                item = {'id': directory.name, 'arm': arm, 'operation': operation, 'source_component': source,
                        'requested_state': state, 'proposal_entry_called': False}
                evidence = {'control': item}
                try:
                    request, authority, incoming = (synthetic_input(source, state, evidence) if arm == 'synthetic' else
                                                    observed_input(source, state, operation, directory, evidence))
                    if arm == 'synthetic' and operation == 'heavy':
                        generated = frozenset(ExactAtomRef(authority.identity, n) for n in sidechain_heavy_mask(authority))
                        request = replace(request, observations=tuple(o for o in request.observations if o.identity not in generated),
                                          ownership=Ownership(generated))
                    evidence['request'] = request
                    item.update(target=authority.identity.sequence_number, state_charge=authority.state.graph.net_charge,
                                declared_heavy_mask=sorted(sidechain_heavy_mask(authority)), source_observations=len(request.observations))
                    if operation == 'H':
                        view = compile_H(request, authority, incoming, evidence)
                        _execute(view, authority, item, evidence)
                    else:
                        execute_heavy(request, authority, arm, directory, item, evidence)
                except ValueError as exc:
                    item.update(status='refused', reason=str(exc))
                _write(directory/'evidence.json', evidence)
                records.append(item)
    check_pins()
    summary = {'schema': 'dvbfixer.all20-heavy-H-coverage.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'manifest_sha256': MANIFEST_SHA256, 'declared': 80, 'executed_entries': len(records),
               'controls': records, 'heavy_H_policies_differ_and_are_not_pooled': True,
               'all20_production_native_quality_or_combined_pipeline_qualified': False, 'production_enabled': False}
    _write(output/'summary.json', summary)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
