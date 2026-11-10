"""Prospective enrollment does not inspect reconstruction/scoring outcomes."""
from __future__ import annotations

import json
from io import StringIO
from pathlib import Path

import pytest
from scripts import prospective_branched_cohort as module
from scripts import prospective_scoring_followup as followup


def atom(serial=1, chain='A', x=1., name='N', residue='VAL', number=1, occupancy=1., altloc=' '):
    return (f'ATOM  {serial:5d} {name:^4}{altloc}{residue:3} {chain}{number:4d}    '
            f'{x:8.3f}{2.:8.3f}{3.:8.3f}{occupancy:6.2f}{20.:6.2f}          {name[0]:>2}  \n')


def test_projection_keeps_first_model_not_later_model_or_lowercase_chain():
    first, second, lowercase = atom(), atom(x=9.), atom(serial=2, chain='a')
    original = ('HEADER test\nMODEL       10\n'+first+lowercase+'TER\nENDMDL\nMODEL       11\n'+second+'ENDMDL\nEND\n').encode()
    projected = module.project(original)
    assert first.encode() in projected and second.encode() not in projected and lowercase.encode() not in projected
    assert b'MODEL ' not in projected and projected.endswith(b'TER\nEND\n')


def test_relevant_cross_chain_metadata_is_retained_not_silently_truncated():
    selected = atom()
    ss = list(' '*80+'\n')
    ss[:6] = 'SSBOND'
    ss[15], ss[29] = 'A', 'B'
    metadata = ''.join(ss)
    conect = 'CONECT    1  100\n'
    original = (metadata+selected+atom(serial=100, chain='B')+'TER\n'+conect+'END\n').encode()
    projected = module.project(original)
    assert metadata.encode() in projected and conect.encode() in projected
    assert atom(serial=100, chain='B').encode() not in projected
    assert projected.index(conect.encode()) > projected.index(selected.encode())


def test_altloc_partial_occupancy_and_explicit_ter_are_not_repaired():
    line = atom(occupancy=.5, altloc='B')
    original = (line+'TER\n'+atom(serial=2, number=2)+'TER\nEND\n').encode()
    assert module.project(original) == original


def test_absent_chain_refuses_without_fallback():
    with pytest.raises(ValueError, match='do not substitute'):
        module.project((atom(chain='a')+'END\n').encode())


def test_enrollment_only_reads_chosen_policy_and_source_hashes(tmp_path, monkeypatch):
    policy = tmp_path/'policy.json'
    policy.write_text(json.dumps(dict(source_hashes={}, selection_policy='fixed rule', hydrogen_protocol='fixed H',
                                      init_flags='fixed flags', expected_rosetta_revision='fixed revision')))
    monkeypatch.setattr(module, 'POLICY', policy)
    monkeypatch.setattr(module, 'implementation_hashes', lambda: {'example': 'frozen'})
    module.enroll(tmp_path/'enrollment.json')
    result = json.loads((tmp_path/'enrollment.json').read_bytes())
    assert [s['id'] for s in result['systems']] == ['1CRN', '1AKI', '1BTA']
    assert result['expected_masks'] == 48 and not result['outcome_filtering']
    assert not result['replacement_structures_allowed']


def test_changed_enrolled_implementation_refuses_before_reading_inputs(tmp_path, monkeypatch):
    record = tmp_path/'enrollment.json'
    record.write_text(json.dumps({'source_hashes': {'example': 'old'}}))
    monkeypatch.setattr(module, 'implementation_hashes', lambda: {'example': 'new'})
    with pytest.raises(ValueError, match='implementation changed'):
        module.freeze(record, tmp_path/'no-inputs', tmp_path/'out')
    assert not (tmp_path/'out').exists()


def test_annotation_order_amendment_is_explicit_and_cannot_change_score_code(tmp_path, monkeypatch):
    record = tmp_path/'enrollment.json'
    data = dict(source_hashes={'scripts/prospective_branched_cohort.py': 'old', 'score.py': 'same'},
                systems=[{'id': 'unchanged'}], score_selection_policy='fixed rule')
    record.write_text(json.dumps(data))
    original = record.read_bytes()
    monkeypatch.setattr(module, 'implementation_hashes', lambda: {'scripts/prospective_branched_cohort.py': 'new', 'score.py': 'same'})
    module.amend(record, tmp_path/'amended.json')
    amended = json.loads((tmp_path/'amended.json').read_bytes())
    assert record.read_bytes() == original and amended['systems'] == data['systems']
    assert amended['score_selection_policy'] == 'fixed rule'
    assert amended['protocol_amendment']['reconstruction_results_already_observed']
    monkeypatch.setattr(module, 'implementation_hashes', lambda: {'scripts/prospective_branched_cohort.py': 'new', 'score.py': 'changed'})
    with pytest.raises(ValueError, match='only projection serialization'):
        module.amend(record, tmp_path/'invalid.json')


@pytest.mark.parametrize('changed', ['score_function', 'seed', 'extra_heavy_policy', 'selection_policy'])
def test_followup_rejects_policy_retuning_including_nonselection_fields(tmp_path, monkeypatch, changed):
    keys = ('score_function', 'seed', 'radii', 'init_flags', 'hydrogen_protocol', 'extra_heavy_policy',
            'virtual_atom_policy', 'selection_policy', 'native_policy', 'topology_scope', 'expected_rosetta_revision')
    policy = {k: 'fixed-'+k for k in keys}
    policy['source_hashes'] = {}
    path = tmp_path/'policy.json'
    path.write_text(json.dumps(policy))
    monkeypatch.setattr(followup, 'POLICY', path)
    enrollment = dict(score_policy_plan_sha256=followup.sha(path.read_bytes()), source_hashes={},
                      score_selection_policy=policy['selection_policy'], hydrogen_protocol=policy['hydrogen_protocol'],
                      init_flags=policy['init_flags'], expected_rosetta_revision=policy['expected_rosetta_revision'])
    followup.check_policy(enrollment, policy)
    with pytest.raises(ValueError, match='scoring policy changed'):
        followup.check_policy(enrollment, dict(policy, **{changed: 'retuned'}))


def test_all_48_masks_freeze_without_admission_or_proposal_success_filtering(tmp_path, monkeypatch):
    inputs = tmp_path/'inputs'
    inputs.mkdir()
    # Deliberately missing backbone/distal atoms and containing a measured H.
    # Admission would fail; enrollment must still retain every declared mask.
    source = ''.join(atom(serial=i, residue=name, number=i) for i, name in enumerate(('VAL', 'THR', 'VAL', 'THR'), 1))
    source += atom(serial=5, residue='VAL', number=1, name='H')+'TER\nEND\n'
    for label in module.SYSTEMS:
        (inputs/(label+'.pdb')).write_text(source)
    hashes = {'example': 'unchanged'}
    monkeypatch.setattr(module, 'implementation_hashes', lambda: hashes)
    enrollment = tmp_path/'enrollment.json'
    enrollment.write_text(json.dumps(dict(source_hashes=hashes, systems=[{'id': n} for n in module.SYSTEMS], expected_masks=48)))
    module.freeze(enrollment, inputs, tmp_path/'plan')
    plan = json.loads((tmp_path/'plan/plan.json').read_bytes())
    assert len(plan['cases']) == 48 and plan['no_own_success_filter']
    assert ' H ' in (tmp_path/'plan/1CRN-reference.pdb').read_text()


def test_projected_conect_is_readable_by_openmm_without_coordinate_changes():
    from openmm.app import PDBFile

    a, b = atom(name='N'), atom(serial=2, name='CA')
    projected = module.project((a+b+'TER\nCONECT    1    2\nEND\n').encode())
    pose = PDBFile(StringIO(projected.decode()))
    assert pose.topology.getNumAtoms() == 2


@pytest.mark.parametrize(('label', 'atoms', 'disulfides'), [('1CRN', 327, 3), ('1AKI', 1001, 4)])
def test_deposited_projection_preserves_disulfides_and_physical_atom_count(label, atoms, disulfides):
    from openmm.app import PDBFile

    source = Path('tests/fixtures/protein_context/scoring_followup', label+'.pdb').read_bytes()
    projected = module.project(source)
    pose = PDBFile(StringIO(projected.decode()))
    assert pose.topology.getNumAtoms() == atoms
    assert sum(a.name == b.name == 'SG' for a, b in pose.topology.bonds()) == disulfides


def test_barstar_observed_hydrogens_are_not_removed_to_make_admission_pass():
    source = Path('tests/fixtures/protein_context/scoring_followup/1BTA.pdb').read_bytes()
    projection = module.project(source)
    hydrogen_lines = [l for l in source.splitlines(keepends=True) if l.startswith(b'ATOM  ') and l[76:78].strip() == b'H']
    assert len(hydrogen_lines) == 714
    assert all(line in projection for line in hydrogen_lines)
