"""Comparison coordinates never acquire a fabricated construction strategy."""
from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest
from scripts.chemical_reconstruction_all20_controls import HEAVY_POLICY, authority_for
from scripts.chemical_reconstruction_ccd_controls import _json

from dvbfixer.chemical_reconstruction import template_comparison
from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate, measure_additions
from dvbfixer.chemical_reconstruction.comparison_metrics import measure_common_heavy
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import Observation, Unit
from dvbfixer.chemical_reconstruction.heldout_scoring import (
    measure_heavy_holdout,
    score_heavy_holdout,
)
from dvbfixer.chemical_reconstruction.observed_peptide import (
    _parse_source,
    compile_peptide_request_plan,
    select_interior_target,
)
from dvbfixer.chemical_reconstruction.research_decode import decode_research_record
from dvbfixer.chemical_reconstruction.source_preparation import prepare_source, prepare_source_probe
from dvbfixer.chemical_reconstruction.template_comparison import (
    SeededHeavyBaseline,
    audit_seeded_baseline,
    inspect_temporary,
    masked_probe_source,
    seeded_heavy_baseline,
)
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref


@pytest.fixture(scope='module')
def case():
    # Historical developer source, NOT any prospective WP2P cohort source.
    data = Path('tests/fixtures/protein_context/1UBQ.pdb').read_bytes()
    identity = select_interior_target(data, 'PRO', chain='A')
    authority = authority_for('PRO', 'PRO', identity, {})
    preparation = prepare_source(data)
    probe = prepare_source_probe(preparation, authority, withheld=frozenset({'CB', 'CG', 'CD'}))
    masked, removed = masked_probe_source(preparation.prepared_text.encode('ascii'), probe.request)
    _, atoms = _parse_source(data)
    reference = tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position, Unit.ANGSTROM)
                      for a in atoms if ExactAtomRef(_component_ref(a), a.name) in probe.withheld)
    return preparation, probe.request, masked, removed, reference


def test_mask_removes_only_exact_withheld_records_and_keeps_context(case):
    preparation, request, masked, removed, reference = case
    assert len(removed) == 3 and all(line.startswith('ATOM  ') for line in removed)
    assert masked == ''.join(line for line in preparation.prepared_text.splitlines(keepends=True) if line not in removed).encode('ascii')
    template_comparison.validate_masked_source(masked, request)
    with pytest.raises(ValueError, match='observation'):
        template_comparison.validate_masked_source(preparation.prepared_text.encode('ascii'), request)
    assert len(reference) == len(request.ownership.generated) == 3


def test_external_measurement_is_same_geometry_without_strategy_forgery(case):
    _, request, _, _, reference = case
    plan = compile_peptide_request_plan(request, HEAVY_POLICY)
    proposal = propose(request, plan)
    assert measure_additions(request, plan, proposal.raw.additions) == audit_candidate(request, plan, proposal.raw)
    score = score_heavy_holdout(request, plan, proposal.raw, reference)
    assert measure_heavy_holdout(request, proposal.raw.additions, reference) == (score.atom_count, score.rmsd_angstrom, score.maximum_angstrom)
    common = measure_common_heavy(request, HEAVY_POLICY, reference, reference)
    assert common.rmsd_angstrom == common.maximum_angstrom == 0.
    assert not hasattr(common, 'strategy')
    decoded = decode_research_record(json.dumps(_json(common)).encode(), type(common))
    assert decoded == common


@pytest.mark.parametrize('mutation', ['partial', 'duplicate', 'element', 'unit'])
def test_no_partial_or_wrong_element_native_score(case, mutation):
    _, request, _, _, reference = case
    additions = reference
    if mutation == 'partial':
        additions = reference[:-1]
    elif mutation == 'duplicate':
        additions = (reference[0], reference[0], reference[2])
    else:
        additions = (replace(reference[0], **({'element': 'N'} if mutation == 'element' else {'unit': Unit.NANOMETER})), *reference[1:])
    with pytest.raises(ValueError):
        measure_common_heavy(request, HEAVY_POLICY, additions, reference)


def test_overlap_pairs_are_unique_and_not_native_acceptance(case):
    _, request, _, _, reference = case
    point = next(o.position for o in request.observations if o.identity.component not in request.targets)
    additions = (replace(reference[0], position=point), *reference[1:])
    measured = measure_common_heavy(request, HEAVY_POLICY, additions, reference)
    assert measured.severe_overlap_pairs and measured.minimum_eligible_separation_angstrom == 0.
    pairs = [(a, b) for a, b, _ in measured.severe_overlap_pairs]
    assert len(pairs) == len(set(pairs))
    assert 'generated-context-minimum-separation' in measured.geometry.findings


def mock_rebuild(case, monkeypatch, *, drift=0., extra=False, d_ca=False):
    import pdbfixer
    from openmm import Vec3
    from openmm.app import PDBFile, Topology, element
    from openmm.unit import angstrom

    import dvbfixer.ffutils.geometry as geometry
    import dvbfixer.profile_readiness as readiness

    _, request, masked, _, reference = case
    pdb = PDBFile(StringIO(masked.decode('ascii')))
    target = next(r for r in pdb.topology.residues() if r.chain.id == 'A' and r.id == '19')
    fixer = SimpleNamespace(topology=pdb.topology, positions=pdb.positions,
                            findMissingAtoms=lambda: None, missingAtoms={target: [SimpleNamespace(name=o.identity.name) for o in reference]},
                            missingTerminals={'unrequested': 'must-be-cleared'})
    calls = []
    def factory(*, pdbfile):
        assert pdbfile.read().encode('ascii') == masked
        return fixer
    def rebuild(actual, **kwargs):
        calls.append(actual)
        assert actual.missingResidues == {} and actual.missingTerminals == {}
        assert {a.name for a in actual.missingAtoms[target]} == {o.identity.name for o in reference}
        old_points = list(actual.positions.value_in_unit(angstrom))
        old_points[0] = old_points[0]+Vec3(drift, 0, 0)
        topology, points = Topology(), []
        for chain in actual.topology.chains():
            new_chain = topology.addChain(chain.id)
            for residue in chain.residues():
                new_residue = topology.addResidue(residue.name, new_chain, residue.id, residue.insertionCode)
                for atom in residue.atoms():
                    topology.addAtom(atom.name, atom.element, new_residue)
                    points.append(old_points[atom.index])
                if residue != target:
                    continue
                ca = next(a for a in residue.atoms() if a.name == 'CA')
                for o in reference:
                    topology.addAtom(o.identity.name, element.carbon, new_residue)
                    point = Vec3(*o.position)
                    if d_ca and o.identity.name == 'CB':
                        point = 2*old_points[ca.index]-point
                    points.append(point)
                if extra:
                    topology.addAtom('XX', element.carbon, new_residue)
                    points.append(Vec3(0, 0, 0))
        actual.topology = topology
        actual.positions = points*angstrom
    monkeypatch.setattr(pdbfixer, 'PDBFixer', factory)
    monkeypatch.setattr(geometry, 'rebuild_missing_atoms_with_retry', rebuild)
    monkeypatch.setattr(readiness, 'pdbfixer_resources_ready', lambda: True)
    return calls


def test_seeded_baseline_calls_shared_helper_once_and_retains_raw_context_snapshot(case, monkeypatch):
    _, request, masked, _, reference = case
    calls = mock_rebuild(case, monkeypatch)
    result = seeded_heavy_baseline(request, masked)
    assert len(calls) == 1 and result.status == 'coordinates-proposed' and result.global_l_passed
    assert not hasattr(result, 'strategy') and {o.identity for o in result.raw} == request.ownership.generated
    assert result.known_maximum_movement_angstrom < 1e-10
    decoded = decode_research_record(json.dumps(_json(result)).encode(), SeededHeavyBaseline)
    assert decoded == result
    audit_seeded_baseline(request, masked, decoded)


@pytest.mark.parametrize('mutation', ['drift', 'extra', 'D'])
def test_seeded_failures_are_retained_not_repaired(case, monkeypatch, mutation):
    _, request, masked, _, _ = case
    calls = mock_rebuild(case, monkeypatch, drift=.1 if mutation == 'drift' else 0., extra=mutation == 'extra', d_ca=mutation == 'D')
    result = seeded_heavy_baseline(request, masked)
    assert len(calls) == 1 and result.status == 'refused'
    if mutation != 'extra':
        assert result.raw and result.temporary_atoms
        audit_seeded_baseline(request, masked, result)


def test_missing_assets_are_unavailable_and_never_call_rebuild(case, monkeypatch):
    import dvbfixer.ffutils.geometry as geometry
    import dvbfixer.profile_readiness as readiness

    _, request, masked, _, _ = case
    monkeypatch.setattr(readiness, 'pdbfixer_resources_ready', lambda: False)
    monkeypatch.setattr(geometry, 'rebuild_missing_atoms_with_retry', lambda *a, **k: pytest.fail('missing assets rebuilt'))
    result = seeded_heavy_baseline(request, masked)
    assert result.status == 'unavailable' and not result.raw and result.reason == 'pdbfixer-resources-not-ready'
    audit_seeded_baseline(request, masked, result)


@pytest.mark.parametrize('mutation', ['raw', 'rounding', 'fixed', 'global-L', 'environment', 'source'])
def test_saved_template_corruption_fails_without_rebuilding(case, monkeypatch, mutation):
    _, request, masked, _, _ = case
    calls = mock_rebuild(case, monkeypatch)
    result = seeded_heavy_baseline(request, masked)
    if mutation == 'raw':
        result = replace(result, raw=result.raw[:-1])
    elif mutation == 'rounding':
        result = replace(result, rounded=result.rounded[:-1])
    elif mutation == 'fixed':
        result = replace(result, known_maximum_movement_angstrom=1.)
    elif mutation == 'global-L':
        result = replace(result, global_l_passed=False)
    elif mutation == 'environment':
        result = replace(result, environment=())
    else:
        result = replace(result, masked_source_sha256='0'*64)
    with pytest.raises(ValueError):
        audit_seeded_baseline(request, masked, result)
    assert len(calls) == 1


def test_unknown_duplicate_or_missing_temporary_atoms_are_not_spliced_away(case):
    _, request, _, _, reference = case
    temporary = request.observations+reference
    assert inspect_temporary(request, reference, temporary) == 0.
    for changed in (temporary[:-1], (*temporary, temporary[0])):
        with pytest.raises(ValueError):
            inspect_temporary(request, reference, changed)


def test_missing_protocol_refuses_before_output_or_science(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_comparison_controls as controls

    monkeypatch.setattr(controls, 'PROTOCOL', tmp_path/'missing')
    with pytest.raises(FileNotFoundError):
        controls.run_controls(tmp_path/'out')
    assert not (tmp_path/'out').exists()


def test_cohort_is_three_systems_nine_dependent_cells():
    from scripts.chemical_reconstruction_comparison_controls import check_pins

    cohort, _ = check_pins()
    assert len(cohort['sources']) == len(cohort['masks']) == 3 and cohort['declared_cells'] == 9
    assert {r['id'] for r in cohort['sources']} == {'1TIM', '1CRN', '2CI2'}
    assert next(r for r in cohort['sources'] if r['id'] == '2CI2')['chain'] == 'I'
    assert all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest() == r['sha256'] for r in cohort['sources'])


def test_source_refusals_keep_all_nine_cells_without_constructors(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_comparison_controls as controls

    def refuse(*args):
        raise ValueError('synthetic-source-refusal')
    monkeypatch.setattr(controls, 'prepare_cell', refuse)
    monkeypatch.setattr(controls, 'propose', lambda *a: pytest.fail('source refusal constructed'))
    monkeypatch.setattr(controls, 'seeded_heavy_baseline', lambda *a: pytest.fail('source refusal rebuilt'))
    summary = controls.run_controls(tmp_path/'run')
    assert len(summary['controls']) == 9
    assert all(r['status'] == 'source-refused' and not r['template_called'] for r in summary['controls'])
    assert controls.readback_controls(tmp_path/'run') == summary


def test_reference_extracted_after_all_available_arms_and_template_runs_when_ours_rejects(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_comparison_controls as controls

    _, catalog = controls.check_pins()
    row = {'id': 'historical-test', 'path': 'tests/fixtures/protein_context/1UBQ.pdb',
           'chain': 'A', 'residue_number': '19', 'insertion_code': ''}
    monkeypatch.setattr(controls, 'check_pins', lambda: ({'sources': [row], 'masks': [['CB', 'CG', 'CD']]}, catalog))
    calls = []
    ordinary = controls.propose
    def rejected(request, plan):
        calls.append('original')
        return replace(ordinary(request, plan), status='rejected', reason='synthetic-negative')
    def template(request, source):
        calls.append('template')
        return SeededHeavyBaseline(request.digest, hashlib.sha256(source).hexdigest(), 'unavailable', 'synthetic-unavailable')
    reference = controls.heldout_after_proposal
    def after(*args):
        assert calls == ['original', 'template']
        calls.append('reference')
        return reference(*args)
    monkeypatch.setattr(controls, 'propose', rejected)
    monkeypatch.setattr(controls, 'seeded_heavy_baseline', template)
    monkeypatch.setattr(controls, 'heldout_after_proposal', after)
    monkeypatch.setattr(controls, 'propose_projection', lambda *a: pytest.fail('rejected original projected'))
    summary = controls.run_controls(tmp_path/'run')
    assert calls == ['original', 'template', 'reference']
    assert not summary['controls'][0]['projection_called'] and summary['controls'][0]['template_called']


def test_unexpected_failure_keeps_evidence_and_is_not_scientific_refusal(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_comparison_controls as controls

    def broken(*args):
        raise RuntimeError('synthetic-programming-error')
    monkeypatch.setattr(controls, 'prepare_cell', broken)
    with pytest.raises(RuntimeError, match='programming'):
        controls.run_controls(tmp_path/'run')
    saved = json.loads((tmp_path/'run/1TIM-CD/evidence.json').read_text())
    assert 'status' not in saved['control'] and not saved['control']['template_called']
    assert not (tmp_path/'run/summary.json').exists()
