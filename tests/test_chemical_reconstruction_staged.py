"""Staged heavy/H ownership, rejected evidence and constructor-free auditing."""
from __future__ import annotations

from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_all20_controls import (
    FFBONDED,
    HEAVY_POLICY,
    PAIRS,
    backbone_types,
    sidechain_heavy_mask,
    synthetic_input,
)

from dvbfixer.chemical_reconstruction import construction, staged_completion
from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.conformer_search import search_conformers
from dvbfixer.chemical_reconstruction.contracts import Ownership, digest_record
from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan
from dvbfixer.chemical_reconstruction.staged_audit import audit_staged_result
from dvbfixer.chemical_reconstruction.staged_completion import complete_staged


def stage_case(source='SER'):
    state = dict(PAIRS)[source]
    full, authority, incoming = synthetic_input(source, state, {})
    generated = frozenset(o.identity for o in full.observations if o.identity.component == authority.identity
                          and o.identity.name in sidechain_heavy_mask(authority))
    request = replace(full, observations=tuple(o for o in full.observations if o.identity not in generated), ownership=Ownership(generated))
    plan = compile_peptide_request_plan(request, HEAVY_POLICY)
    roots = tuple(sorted(o.identity for o in request.observations if o.identity.component in request.targets))
    parent = search_conformers(request, plan, authority, roots)
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    types = None if state == 'PRO' else backbone_types(state)
    return request, parent, authority, incoming, parameters, types


def run_case(case):
    request, parent, authority, incoming, parameters, types = case
    return complete_staged(request, parent, authority, incoming, parameters, backbone_types=types)


def readback(case, result):
    return audit_staged_result(*case, result)


@pytest.mark.parametrize('source', ['ALA', 'SER', 'THR', 'PRO', 'HIS', 'ILE'])
def test_staged_complete_masks_fixed_checkpoint_and_measured_disposition(source):
    case = stage_case(source)
    request, parent, authority, *_ = case
    parent_before, source_before = digest_record(parent), digest_record(request)
    result = run_case(case)
    assert result.status in {'staged-geometry-passed', 'H-rejected'}
    assert result.view is not None and result.hydrogen_result is not None
    assert not readback(case, result)
    assert result.raw is not None and result.rounded is not None
    checkpoint = parent.attempts[parent.selected_index].rounded
    assert result.raw.generated_heavy is checkpoint.additions
    assert result.rounded.generated_heavy is checkpoint.additions
    assert {a.identity for a in result.raw.generated_heavy} == request.ownership.generated
    assert {a.identity.name for a in result.raw.generated_hydrogens} == {a.name for a in authority.state.graph.atoms if a.element == 'H'}
    assert result.view.request.original is request
    assert digest_record(parent) == parent_before and digest_record(request) == source_before
    assert not result.md_ready and not result.native_quality_qualified


def test_independent_stage_audit_never_calls_any_coordinate_constructor(monkeypatch):
    from rdkit.Chem import rdMolTransforms

    from dvbfixer.chemical_reconstruction import (
        conformer_search,
        joint_hydrogen_geometry,
        torsion_geometry,
    )

    case = stage_case()
    result = run_case(case)
    def forbidden(*args, **kwargs):
        pytest.fail('final audit invoked coordinate construction')
    for module, name in ((staged_completion, 'propose'), (construction, 'propose'), (construction, '_internal'),
                         (conformer_search, 'search_conformers'), (conformer_search, 'rotate_fragments'),
                         (torsion_geometry, 'rotate_fragments'), (joint_hydrogen_geometry, 'joint_group_points'),
                         (joint_hydrogen_geometry, 'fit_joint_directions'), (rdMolTransforms, 'SetDihedralDeg')):
        monkeypatch.setattr(module, name, forbidden)
    assert not readback(case, result)


def test_H_failure_is_retained_and_does_not_reselect_or_move_heavy(monkeypatch):
    case = stage_case()
    parent = case[1]
    before = digest_record(parent)
    calls = []
    ordinary = staged_completion.propose
    def reject_H(request, plan):
        calls.append(request)
        proposal = ordinary(request, plan)
        return replace(proposal, status='rejected', reason='synthetic-H-negative-no-retry')
    monkeypatch.setattr(staged_completion, 'propose', reject_H)
    result = run_case(case)
    assert len(calls) == 1 and result.status == 'H-rejected' and not result.audit.passed
    assert result.hydrogen_result.reason == 'synthetic-H-negative-no-retry'
    assert digest_record(parent) == before
    assert result.raw.generated_heavy == parent.attempts[parent.selected_index].rounded.additions
    assert not readback(case, result)


@pytest.mark.parametrize('mutation', ['source-binding', 'parent-binding', 'H-moved', 'heavy-moved', 'wrong-phase',
                                     'wrong-status', 'wrong-reason', 'saved-audit', 'admission', 'rounding', 'view-plan', 'mask'])
def test_final_readback_rejects_metadata_geometry_rounding_and_role_corruption(mutation):
    case = stage_case()
    result = run_case(case)
    if mutation in {'source-binding', 'parent-binding'}:
        result = replace(result, **{'original_request_digest' if mutation == 'source-binding' else 'parent_digest': '0'*64})
    elif mutation in {'H-moved', 'heavy-moved', 'mask'}:
        field = 'generated_heavy' if mutation == 'heavy-moved' else 'generated_hydrogens'
        atoms = getattr(result.raw, field)
        altered = atoms[:-1] if mutation == 'mask' else (replace(atoms[0], position=(999., 999., 999.)), *atoms[1:])
        result = replace(result, raw=replace(result.raw, **{field: altered}))
    elif mutation == 'wrong-phase':
        result = replace(result, raw=replace(result.raw, phase='rounded-H-after-rounded-heavy'))
    elif mutation == 'wrong-status':
        result = replace(result, status='H-refused')
    elif mutation == 'wrong-reason':
        result = replace(result, reason='forged success')
    elif mutation == 'saved-audit':
        result = replace(result, audit=replace(result.audit, findings=('forged',)))
    elif mutation == 'admission':
        proposal = result.hydrogen_result
        result = replace(result, hydrogen_result=replace(proposal, admission=replace(proposal.admission, checked_constraints=())))
    elif mutation == 'rounding':
        proposal = result.hydrogen_result
        result = replace(result, hydrogen_result=replace(proposal, rounded=proposal.raw))
    else:
        view = result.view
        result = replace(result, view=replace(view, plan=replace(view.plan, policy=replace(view.plan.policy, angle_tolerance_degrees=14.))))
    assert readback(case, result)


def test_pair_coverage_includes_new_heavy_H_pairs_and_source_context():
    from itertools import combinations

    case = stage_case('ILE')
    result = run_case(case)
    coverage = result.audit.pair_coverage
    assert coverage.heavy_source > 0 and coverage.hydrogen_source > 0
    assert coverage.hydrogen_heavy > 0 and coverage.hydrogen_hydrogen > 0
    assert not readback(case, result)
    request = case[0]
    from dvbfixer.domain.structure_identity import ExactAtomRef
    edges = [{ExactAtomRef(c.identity, b.atom1), ExactAtomRef(c.identity, b.atom2)}
             for c in request.components if c.resolved is not None for b in c.resolved.graph.bonds]
    edges += [{c.link.component_atom, c.link.other_atom} for c in request.connections if c.kind.value == 'covalent']
    def neighbours(atom):
        return {other for edge in edges if atom in edge for other in edge if other != atom}
    source = {a.identity for a in request.observations}
    heavy = request.ownership.generated
    hydrogen = result.view.request.ownership.generated
    expected = {'hs': 0, 'hh': 0, 'Hs': 0, 'Hh': 0, 'HH': 0}
    for a, b in combinations(sorted(source | heavy | hydrogen), 2):
        if a in source and b in source:
            continue
        if b in neighbours(a) or neighbours(a) & neighbours(b):
            continue
        if a in hydrogen or b in hydrogen:
            key = 'HH' if a in hydrogen and b in hydrogen else 'Hh' if a in heavy or b in heavy else 'Hs'
        else:
            key = 'hh' if a in heavy and b in heavy else 'hs'
        expected[key] += 1
    assert (coverage.heavy_source, coverage.heavy_heavy, coverage.hydrogen_source, coverage.hydrogen_heavy,
            coverage.hydrogen_hydrogen) == tuple(expected[k] for k in ('hs', 'hh', 'Hs', 'Hh', 'HH'))


def test_bad_selected_parent_fails_before_H_proposal(monkeypatch):
    case = list(stage_case())
    case[1] = replace(case[1], selected_index=999)
    monkeypatch.setattr(staged_completion, 'propose', lambda *args: pytest.fail('proposed H from invalid parent'))
    result = run_case(case)
    assert result.status == 'parent-refused' and result.view is None and result.raw is None
    assert not readback(case, result)


def test_strict_known_DTO_roundtrip_then_constructor_free_readback():
    import json

    from scripts.chemical_reconstruction_ccd_controls import _json

    from dvbfixer.chemical_reconstruction.conformer_search import ConformerSearchResult
    from dvbfixer.chemical_reconstruction.contracts import ReconstructionRequest
    from dvbfixer.chemical_reconstruction.declared_peptide import DeclaredPeptideAuthority
    from dvbfixer.chemical_reconstruction.research_decode import decode_research_record
    from dvbfixer.chemical_reconstruction.staged_contracts import StagedResult

    case = stage_case('PRO')
    result = run_case(case)
    request, parent, authority, incoming, parameters, types = case
    restored = []
    for item, expected in ((request, ReconstructionRequest), (parent, ConformerSearchResult),
                           (authority, DeclaredPeptideAuthority), (result, StagedResult)):
        decoded = decode_research_record(json.dumps(_json(item)).encode(), expected)
        assert digest_record(decoded) == digest_record(item)
        restored.append(decoded)
    assert not readback((*restored[:3], incoming, parameters, types), restored[3])


@pytest.mark.parametrize('mutation', ['unknown-field', 'missing-field', 'bool-index', 'float-index', 'unknown-enum',
                                     'duplicate-set', 'wrong-container', 'wrong-literal'])
def test_decoder_rejects_unknown_or_inexact_transport(mutation):
    import json

    from scripts.chemical_reconstruction_ccd_controls import _json

    from dvbfixer.chemical_reconstruction.conformer_search import ConformerSearchResult
    from dvbfixer.chemical_reconstruction.research_decode import decode_research_record

    case = stage_case()
    saved = _json(case[1])
    if mutation == 'unknown-field':
        saved['native_hint'] = 1.
    elif mutation == 'missing-field':
        del saved['plan']
    elif mutation == 'bool-index':
        saved['selected_index'] = True
    elif mutation == 'float-index':
        saved['selected_index'] = 0.
    elif mutation == 'unknown-enum':
        saved['plan']['strategy'] = 'pickle-python-import'
    elif mutation == 'duplicate-set':
        # Use a known frozenset field in the authority-bound baseline request elsewhere.
        from dvbfixer.chemical_reconstruction.contracts import ReconstructionRequest
        request = _json(case[0])
        request['targets'].append(request['targets'][0])
        with pytest.raises(ValueError):
            decode_research_record(json.dumps(request).encode(), ReconstructionRequest)
        return
    elif mutation == 'wrong-container':
        saved['attempts'] = {}
    else:
        saved['baseline']['status'] = 'qualified'
    with pytest.raises((ValueError, TypeError)):
        decode_research_record(json.dumps(saved).encode(), ConformerSearchResult)


@pytest.mark.parametrize('data', [b'{', b'{"x":1,"x":2}', b'NaN', b'Infinity', b'1e9999', b'\xff'])
def test_private_JSON_reader_rejects_malformed_duplicate_and_nonfinite_values(data):
    from dvbfixer.chemical_reconstruction.research_decode import parse_research_json

    with pytest.raises(ValueError):
        parse_research_json(data)


def test_private_JSON_byte_depth_and_node_budgets(monkeypatch):
    from dvbfixer.chemical_reconstruction import research_decode

    monkeypatch.setattr(research_decode, 'MAX_BYTES', 10)
    with pytest.raises(ValueError, match='byte-budget'):
        research_decode.parse_research_json(b' ' * 11)
    monkeypatch.setattr(research_decode, 'MAX_BYTES', 100)
    monkeypatch.setattr(research_decode, 'MAX_DEPTH', 2)
    with pytest.raises(ValueError, match='depth-budget'):
        research_decode.parse_research_json(b'[[[[]]]]')
    monkeypatch.setattr(research_decode, 'MAX_NODES', 3)
    with pytest.raises(ValueError, match='node-or-depth'):
        research_decode.parse_research_json(b'[0,0,0]')


def test_unknown_root_class_is_not_imported_or_constructed():
    from dvbfixer.chemical_reconstruction.research_decode import decode_research_record

    class Foreign:
        def __init__(self):
            pytest.fail('foreign class constructed')
    with pytest.raises(ValueError, match='unknown-research-root'):
        decode_research_record(b'{}', Foreign)


def test_private_reader_rejects_numeric_overflow_before_DTO_construction():
    from dvbfixer.chemical_reconstruction.research_decode import parse_research_json

    with pytest.raises(ValueError, match='out-of-range'):
        parse_research_json(b'9' * 400)


def test_saved_parent_driver_refuses_missing_wrong_hash_or_oversized_input(tmp_path, monkeypatch):
    import hashlib

    from scripts import chemical_reconstruction_staged_controls as controls

    path = tmp_path/'parent.json'
    with pytest.raises(FileNotFoundError):
        controls._read_pinned(path, '0'*64)
    path.write_bytes(b'{}')
    with pytest.raises(ValueError, match='hash-mismatch'):
        controls._read_pinned(path, '0'*64)
    assert controls._read_pinned(path, hashlib.sha256(b'{}').hexdigest()) == b'{}'
    monkeypatch.setattr(controls, 'MAX_BYTES', 1)
    with pytest.raises(ValueError, match='byte-budget'):
        controls._read_pinned(path, hashlib.sha256(b'{}').hexdigest())


def test_missing_saved_parent_fails_before_output_reservation_or_science(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_staged_controls as controls

    monkeypatch.setattr(controls, 'check_pins', lambda: {'entries': [{'parent_artifact': str(tmp_path/'absent'), 'parent_sha256': '0'*64}]})
    monkeypatch.setattr(controls, 'complete_staged', lambda *args, **kwargs: pytest.fail('fabricated missing parent'))
    output = tmp_path/'output'
    with pytest.raises(FileNotFoundError):
        controls.run_controls(output)
    assert not output.exists()


def test_parent_baseline_saved_metadata_is_remeasured_not_trusted():
    case = list(stage_case())
    parent = case[1]
    baseline = replace(parent.baseline, plan_digest='0'*64)
    case[1] = replace(parent, baseline=baseline)
    result = run_case(case)
    assert result.status == 'parent-refused' and 'baseline-metadata' in result.reason
    assert not readback(case, result)


def test_partial_H_constructor_failure_keeps_raw_negative_without_false_pair_coverage(monkeypatch):
    from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate

    case = stage_case()
    ordinary = staged_completion.propose
    def partial(request, plan):
        proposal = ordinary(request, plan)
        raw = replace(proposal.raw, additions=proposal.raw.additions[:-1])
        return replace(proposal, status='construction-failed', reason='synthetic-partial-H-negative',
                       raw=raw, raw_audit=audit_candidate(request, plan, raw), rounded=None, rounded_audit=None)
    monkeypatch.setattr(staged_completion, 'propose', partial)
    result = run_case(case)
    assert result.status == 'H-rejected' and result.raw is not None and result.rounded is None
    assert result.audit.pair_coverage is None
    assert 'staged-pair-H-inventory-incomplete' in result.audit.findings
    assert not readback(case, result)
