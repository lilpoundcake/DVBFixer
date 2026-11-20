"""Incoming carbonyl type is not a protonation, charge, or source-byte proof."""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest
from scripts.chemical_reconstruction_all20_controls import FFBONDED, authority_for, load_parent
from scripts.chemical_reconstruction_ccd_controls import _json

from dvbfixer.chemical_reconstruction.incoming_source import (
    audit_incoming_source,
    bind_incoming_source,
)
from dvbfixer.chemical_reconstruction.incoming_type import (
    IncomingCarbonylTypeProof,
    candidate_family,
    invariant_carbonyl_type,
    load_incoming_type_catalog,
    prove_incoming_type,
)
from dvbfixer.chemical_reconstruction.observed_peptide import select_interior_target
from dvbfixer.chemical_reconstruction.planar_port import (
    PlanarPortDeclaration,
    load_planar_port_authority,
    planar_port_constraints,
)
from dvbfixer.chemical_reconstruction.research_decode import decode_research_record
from dvbfixer.chemical_reconstruction.source_preparation import prepare_source, prepare_source_probe
from dvbfixer.chemical_reconstruction.type_proven_port import TypeProvenPlanarPortDeclaration
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.ffutils.variants import AMBER_VARIANT_TO_PARENT


@pytest.fixture(scope='module')
def case():
    data = Path('tests/fixtures/regressions/1FR2.pdb').read_bytes()
    authority = authority_for('PRO', 'PRO', select_interior_target(data, 'PRO', chain='A'), {})
    preparation = prepare_source(data)
    probe = prepare_source_probe(preparation, authority, withheld=frozenset({'CB', 'CG', 'CD'}))
    source = bind_incoming_source(preparation, probe, authority)
    catalog = load_incoming_type_catalog(Path('FF/amber19sb.ff/aminoacids.rtp'))
    proof = prove_incoming_type(probe.request, source, catalog, load_parent('HIS'))
    return preparation, probe, authority, proof


def test_HIS_complete_family_without_charge_or_state_choice(case):
    preparation, probe, authority, proof = case
    assert proof.candidate_states == ('HID', 'HIE', 'HIP')
    charges = {s.name: next(row[2] for row in s.atoms if row[0] == 'C')
               for s in proof.catalog.states if s.name in proof.candidate_states}
    assert charges == {'HID': .5973, 'HIE': .5973, 'HIP': .7341}
    assert proof.atom_type == 'C' and proof.assigned_charge is proof.selected_state is None
    donor = next(c for c in proof.request.components if c.identity == proof.source.incoming.component)
    assert not donor.states and donor.selected_state is None
    assert not audit_incoming_source(preparation, probe, authority, proof.source)


@pytest.mark.parametrize('label,expected', [
    ('HIS', ('HID', 'HIE', 'HIP')), ('HID', ('HID',)), ('HIE', ('HIE',)), ('HIP', ('HIP',)),
    ('CYS', ('CYM', 'CYS', 'CYX')), ('ASP', ('ASH', 'ASP')), ('GLU', ('GLH', 'GLU')),
    ('LYS', ('LYN', 'LYS')), ('ALA', ('ALA',)),
])
def test_family_comes_from_existing_owner(case, label, expected):
    catalog = case[3].catalog
    parent, names = candidate_family(label, catalog)
    assert names == expected and parent == AMBER_VARIANT_TO_PARENT.get(label, label)
    assert invariant_carbonyl_type(tuple(s for s in catalog.states if s.name in names)) == 'C'


@pytest.mark.parametrize('label', ['HSD', 'NHID', 'CHIE', 'his', 'UNK', ''])
def test_unknown_terminal_and_other_namespace_are_not_substituted(case, label):
    with pytest.raises(ValueError):
        candidate_family(label, case[3].catalog)


@pytest.mark.parametrize('mutation', ['missing', 'type', 'charge', 'bond', 'owner', 'pin'])
def test_catalog_cannot_omit_or_change_any_state(case, mutation):
    catalog = case[3].catalog
    if mutation == 'missing':
        changes = {'states': tuple(s for s in catalog.states if s.name != 'HIP')}
    elif mutation in {'owner', 'pin'}:
        changes = {'owner_digest' if mutation == 'owner' else 'source_sha256': '0'*64}
    else:
        state = next(s for s in catalog.states if s.name == 'HIP')
        if mutation == 'bond':
            state = replace(state, bonds=state.bonds[:-1])
        else:
            atoms = tuple((n, 'CT' if mutation == 'type' and n == 'C' else t,
                           0. if mutation == 'charge' and n == 'C' else q, g) for n, t, q, g in state.atoms)
            state = replace(state, atoms=atoms)
        changes = {'states': tuple(state if s.name == state.name else s for s in catalog.states)}
    with pytest.raises(ValueError, match='catalog'):
        replace(catalog, **changes)


def test_owner_mutation_fails_closed(case, monkeypatch):
    monkeypatch.setitem(AMBER_VARIANT_TO_PARENT, 'HIP', 'ALA')
    with pytest.raises(ValueError, match='owner'):
        case[3].validate()


@pytest.mark.parametrize('mutation', ['type', 'extra-bond', 'missing-O', 'improper', 'empty'])
def test_local_topology_and_common_type_are_required(case, mutation):
    state = next(s for s in case[3].catalog.states if s.name == 'HIP')
    if mutation == 'type':
        state = replace(state, atoms=tuple((n, 'CT' if n == 'C' else t, q, g) for n, t, q, g in state.atoms))
    elif mutation == 'extra-bond':
        state = replace(state, bonds=tuple(sorted((*state.bonds, ('C', 'CB')))))
    elif mutation == 'missing-O':
        state = replace(state, atoms=tuple(r for r in state.atoms if r[0] != 'O'))
    elif mutation == 'improper':
        state = replace(state, impropers=())
    with pytest.raises(ValueError):
        invariant_carbonyl_type(() if mutation == 'empty' else (state,))


@pytest.mark.parametrize('field,value', [
    ('candidate_states', ('HID',)), ('canonical_parent', 'ALA'), ('atom_type', 'CT'),
    ('selected_state', 'HID'), ('assigned_charge', .5973), ('scope', 'whole-residue'),
])
def test_proof_never_assigns_a_state_or_charge(case, field, value):
    with pytest.raises(ValueError):
        replace(case[3], **{field: value})


@pytest.mark.parametrize('field', ['request_digest', 'observations_digest', 'prepared_source_sha256', 'original_source_sha256', 'preparation_digest'])
def test_source_and_request_hash_binding(case, field):
    proof = case[3]
    with pytest.raises(ValueError):
        replace(proof, source=replace(proof.source, **{field: '0'*64}))


def test_forged_literal_label_is_not_source_qualification(case):
    preparation, probe, authority, proof = case
    # Exact-label semantics are pure. Only a byte audit can detect this lie.
    source = replace(proof.source, residue_name='HID')
    forged = prove_incoming_type(proof.request, source, proof.catalog, proof.ccd)
    assert forged.candidate_states == ('HID',)
    assert audit_incoming_source(preparation, probe, authority, source) == ('incoming-source-byte-slice-binding',)


def test_extra_attachment_and_missing_core_refused(case):
    proof = case[3]
    c = proof.source.incoming
    pair = tuple(sorted((c, ExactAtomRef(c.component, 'CB'))))
    source = replace(proof.source, annotation_pairs=tuple(sorted({*proof.source.annotation_pairs, pair})))
    with pytest.raises(ValueError, match='attachment'):
        replace(proof, source=source)
    with pytest.raises(ValueError, match='four'):
        replace(proof.source, core=proof.source.core[:-1])


@pytest.mark.parametrize('mutation', ['parent', 'leaving', 'order', 'charge'])
def test_CCD_local_footprint_corruption_refused(case, mutation):
    proof = case[3]
    ccd = proof.ccd
    if mutation == 'parent':
        ccd = replace(ccd, component=replace(ccd.component, component_id='ALA'))
    elif mutation == 'leaving':
        ccd = replace(ccd, component=replace(ccd.component, atoms=tuple(replace(a, leaving=False) for a in ccd.component.atoms)))
    else:
        graph = ccd.state.graph
        if mutation == 'order':
            graph = replace(graph, bonds=tuple(replace(b, order=1.) if {b.atom1, b.atom2} == {'C', 'O'} else b for b in graph.bonds))
        else:
            graph = replace(graph, atoms=tuple(replace(a, formal_charge=1) if a.name == 'C' else a for a in graph.atoms))
        with pytest.raises(ValueError):
            replace(proof, ccd=replace(ccd, state=replace(ccd.state, graph=graph)))
        return
    with pytest.raises(ValueError):
        replace(proof, ccd=ccd)


def test_v2_port_is_explicit_and_roundtrips_without_solver(case, monkeypatch):
    from scipy import optimize

    _, probe, authority, proof = case
    parameters = load_planar_port_authority(FFBONDED)
    n, ca, cd = (ExactAtomRef(authority.identity, name) for name in ('N', 'CA', 'CD'))
    old = PlanarPortDeclaration(proof.source.incoming, n, ca, cd, ('C', 'N', 'XC', 'CT'), authority.digest, parameters.digest)
    new = TypeProvenPlanarPortDeclaration(proof.source.incoming, n, ca, cd, old.atom_types,
                                        authority.digest, parameters.digest, incoming_proof=proof)
    assert new.digest != old.digest
    monkeypatch.setattr(optimize, 'least_squares', lambda *a, **k: pytest.fail('proof ran solver'))
    old_constraints = planar_port_constraints(probe.request, authority, old, parameters)
    new_constraints = planar_port_constraints(probe.request, authority, new, parameters)
    assert tuple(replace(c, label='') for c in old_constraints) == tuple(replace(c, label='') for c in new_constraints)
    assert old_constraints[0].label != new_constraints[0].label
    decoded = decode_research_record(json.dumps(_json(new)).encode(), PlanarPortDeclaration)
    assert isinstance(decoded, TypeProvenPlanarPortDeclaration) and decoded.digest == new.digest
    saved_proof = decode_research_record(json.dumps(_json(proof)).encode(), IncomingCarbonylTypeProof)
    saved_proof.validate()
    assert saved_proof.digest == proof.digest
    with pytest.raises(ValueError):
        new.validate_incoming_type(replace(probe.request, observations=tuple(reversed(probe.request.observations))))


@pytest.mark.parametrize('label', ['HID', 'HIE', 'HIP'])
def test_explicit_source_variant_is_literal_not_rewritten(case, label):
    preparation, _, authority, proof = case
    donor = proof.source.incoming.component
    lines = []
    for line in preparation.original_text.splitlines(keepends=True):
        if line.startswith('ATOM  ') and line[21:22] == donor.chain_id and line[22:26].strip() == donor.sequence_number:
            line = line[:17] + label + line[20:]
        lines.append(line)
    prepared = prepare_source(''.join(lines).encode('ascii'))
    probe = prepare_source_probe(prepared, authority, withheld=frozenset({'CB', 'CG', 'CD'}))
    source = bind_incoming_source(prepared, probe, authority)
    assert source.residue_name == label and not audit_incoming_source(prepared, probe, authority, source)
    typed = prove_incoming_type(probe.request, source, proof.catalog, proof.ccd)
    assert typed.candidate_states == (label,) and typed.selected_state is typed.assigned_charge is None


def test_probe_annotation_forgery_is_not_source_admission(case):
    preparation, probe, authority, proof = case
    c = proof.source.incoming
    forged = replace(probe, unassessed_internal_connections=(*probe.unassessed_internal_connections, (c, ExactAtomRef(c.component, 'CB'))))
    with pytest.raises(ValueError):
        bind_incoming_source(preparation, forged, authority)


def test_actual_extra_carbonyl_CONECT_is_not_ignored(case):
    preparation, _, authority, proof = case
    donor = proof.source.incoming.component
    serials = {line[12:16].strip(): int(line[6:11]) for line in preparation.original_text.splitlines()
               if line.startswith('ATOM  ') and line[21:22] == donor.chain_id and line[22:26].strip() == donor.sequence_number}
    extra = f"CONECT{serials['C']:5d}{serials['CB']:5d}\n"
    prepared = prepare_source((preparation.original_text+extra).encode('ascii'))
    probe = prepare_source_probe(prepared, authority, withheld=frozenset({'CB', 'CG', 'CD'}))
    source = bind_incoming_source(prepared, probe, authority)
    with pytest.raises(ValueError, match='attachment'):
        prove_incoming_type(probe.request, source, proof.catalog, proof.ccd)


@pytest.mark.parametrize('mutation', ['unit', 'element', 'position', 'link', 'donor-chain'])
def test_core_or_external_link_cannot_be_changed_behind_source_binding(case, mutation):
    from dvbfixer.chemical_reconstruction.contracts import Unit, digest_record

    proof = case[3]
    request = proof.request
    if mutation == 'link':
        request = replace(request, connections=tuple(c for c in request.connections if proof.source.incoming not in {c.link.component_atom, c.link.other_atom}))
    elif mutation == 'donor-chain':
        source = replace(proof.source, incoming=replace(proof.source.incoming, component=replace(proof.source.incoming.component, chain_id='a')))
        with pytest.raises(ValueError):
            replace(proof, source=source)
        return
    else:
        changes = {'unit': {'unit': Unit.NANOMETER}, 'element': {'element': 'N'}, 'position': {'position': (1., 2., 3.)}}[mutation]
        try:
            request = replace(request, observations=tuple(replace(o, **changes) if o.identity == proof.source.incoming else o for o in request.observations))
        except ValueError:
            return
    source = replace(proof.source, request_digest=request.digest, observations_digest=digest_record(request.observations))
    with pytest.raises(ValueError):
        replace(proof, request=request, source=source)
