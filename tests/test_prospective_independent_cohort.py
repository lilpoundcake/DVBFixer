"""Identity-only pre-download selection and exact-byte first-chain projection."""
import copy
import json
from dataclasses import asdict

import pytest
from scripts import prospective_independent_cohort as module


def atom(serial, name='CA', residue='VAL', chain='a', resid=1, altloc=' ', occupancy=.5):
    return (f'ATOM  {serial:5d} {name:>4}{altloc}{residue:>3} {chain}{resid:4d}    '
            f'{1.:8.3f}{2.:8.3f}{3.:8.3f}{occupancy:6.2f}{10.:6.2f}          {name[0]:>2}\n')


def test_first_chain_is_case_sensitive_and_preserves_h_altloc_and_occupancy():
    first = atom(1, altloc='A')
    hydrogen = atom(2, name='H', altloc='A')
    unrelated = atom(3, chain='A')
    source = ('HEADER    CONTROL\n'+first+hydrogen+'TER\n'+unrelated+'TER\nEND\n').encode()
    output, chain = module.project_first(source)
    assert chain == 'a' and output == (first+hydrogen+'TER\nEND\n').encode()


def test_first_coordinate_model_not_best_or_most_complete_model():
    first = atom(1, chain='B')
    second = atom(2, chain='A', occupancy=1.)
    output, chain = module.project_first(('MODEL        1\n'+first+'ENDMDL\nMODEL        2\n'+second+'ENDMDL\nEND\n').encode())
    assert chain == 'B' and output == (first+'TER\nEND\n').encode()
    with pytest.raises(ValueError, match='no ATOM chain'):
        module.project_first(('MODEL        1\nENDMDL\nMODEL        2\n'+second+'ENDMDL\n').encode())


def test_conect_is_after_coordinates_and_excluded_endpoints_are_not_silently_deleted():
    original = ('CONECT    1    9\n'+atom(1)+atom(9, chain='A')+'END\n').encode()
    output, chain = module.project_first(original)
    assert chain == 'a'
    assert output == (atom(1)+'TER\nCONECT    1    9\nEND\n').encode()


def test_selector_counts_identity_without_geometry_or_occupancy_filtering_and_retains_shortfall():
    source = ''.join([atom(1, resid=20, altloc='A'), atom(2, name='CB', resid=20, altloc='B'),
                      atom(3, residue='THR', resid=21), atom(4, resid=22), atom(5, resid=23)]).encode()
    targets, counts = module.selected_targets(source)
    assert counts == {'VAL': 2, 'THR': 1}
    assert [(n, r.sequence_number) for n, r in targets] == [('VAL', '20'), ('THR', '21'), ('VAL', '22')]


def test_enrollment_reads_no_coordinates_and_calls_no_proposal_or_preflight(tmp_path, monkeypatch):
    def never(*args, **kwargs):
        raise AssertionError('pre-download enrollment must not inspect scientific outcomes')

    for name in ('project_first', 'selected_targets', 'complete_sidechain', 'template_proposal', 'preflight'):
        monkeypatch.setattr(module, name, never)
    policy = tmp_path/'policy.json'
    policy.write_text(json.dumps(dict(source_hashes={}, selection_policy='unchanged', hydrogen_protocol='unchanged',
                                      init_flags='unchanged', expected_rosetta_revision='unchanged')))
    monkeypatch.setattr(module, 'POLICY', policy)
    monkeypatch.setattr(module, 'implementation_hashes', lambda: {})
    output = tmp_path/'enrollment.json'
    module.enroll(output)
    record = json.loads(output.read_bytes())
    assert [s['id'] for s in record['systems']] == ['2CI2', '1R69', '1CTF']
    assert record['max_cases'] == 48 and not record['replacement_structures_allowed']
    with pytest.raises(FileExistsError):
        module.enroll(output)


def test_template_readback_is_exact_and_rejects_moved_observations_and_changed_additions(tmp_path):
    from scripts.audit_independent_cohort import verify_template

    from dvbfixer.domain.structure_identity import ExactAtomRef
    from dvbfixer.protein_completion.contracts import Addition, CompletionRequest, Observation

    target = module.ComponentInstanceRef(1, 'I', '22', 'A')
    observed = tuple(Observation(ExactAtomRef(target, name), name[0], (float(i), 0., 0.)) for i, name in enumerate(('N', 'CA', 'C', 'O')))
    generated = tuple(Addition(ExactAtomRef(target, name), 'C', (float(i), 1., 0.)) for i, name in enumerate(('CB', 'CG1', 'CG2')))
    request = CompletionRequest('a'*64, target, 'VAL', observed, ('CB', 'CG1', 'CG2'))
    snapshot = dict(seed=1, atoms=[(asdict(a.atom), a.element, a.xyz) for a in (*observed, *generated)], observed_max_movement_angstrom=0.)
    proposal = json.loads(json.dumps(dict(additions=[asdict(a) for a in generated], temporary_extra_atoms=[], snapshots=[snapshot])))
    path = tmp_path/'template-seed-1.json'
    path.write_text(json.dumps(proposal['snapshots'][0]))
    assert verify_template(request, proposal, tmp_path) == generated
    changed = copy.deepcopy(proposal)
    changed['additions'][0]['xyz'][0] += 1.
    with pytest.raises(ValueError, match='last snapshot'):
        verify_template(request, changed, tmp_path)
    changed = copy.deepcopy(proposal)
    changed['snapshots'][0]['atoms'][0][2][0] += 1.
    path.write_text(json.dumps(changed['snapshots'][0]))
    with pytest.raises(ValueError, match='moved observations'):
        verify_template(request, changed, tmp_path)


@pytest.mark.parametrize('identifier,chain,counts', [('2CI2', 'I', {'VAL': 2, 'THR': 2}),
                                                   ('1R69', 'A', {'VAL': 2, 'THR': 2}),
                                                   ('1CTF', 'A', {'VAL': 2, 'THR': 1})])
def test_tracked_prospective_sources_retain_actual_chain_and_selection_shortfall(identifier, chain, counts):
    original = module.Path(f'tests/fixtures/protein_context/independent_followup/{identifier}.pdb').read_bytes()
    projected, chosen_chain = module.project_first(original)
    targets, actual = module.selected_targets(projected)
    assert chosen_chain == chain and actual == counts
    assert len(targets) == sum(counts.values())
