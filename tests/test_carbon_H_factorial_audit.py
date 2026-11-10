"""Independent fixed-zero slots and matched denominator retention."""
import copy
import gzip
import json

import pytest
from scripts.audit_carbon_H_factorial import check_factorial_initial, known_replay, summarize
from scripts.research_carbon_H_factorial import initialize
from scripts.research_dunbrack_prior import sha


def fixture():
    rows = [dict(key=['a', '82', 'A', 'ALA', n], element=e) for n, e in
        [('CB', 'C'), ('CA', 'C'), ('C', 'C'), ('HB1', 'H'), ('HB2', 'H'), ('HB3', 'H')]]
    points = [[0., 0., 0.], [.15, 0., 0.], [.15, .15, 0.],
        [0., 0., .1], [0., .01, .11], [0., .02, .09]]
    case = dict(methyl=dict(groups=[dict(parent=0, parent_key=rows[0]['key'], axis=1, frame=2, hydrogens=[3, 4, 5])]),
        methylene=dict(groups=[]))
    return rows, points, case


def test_independent_fixed_zero_geometry_accepts_original_constructor():
    rows, points, case = fixture()
    check_factorial_initial(rows, points, initialize(points, case, 'all_methyl_zero'), case, 'all_methyl_zero')
    check_factorial_initial(rows, points, points, case, 'unchanged')


@pytest.mark.parametrize('fault', ['heavy', 'radius', 'phase', 'slots', 'incomplete'])
def test_unowned_coordinate_radius_phase_or_slot_drift_is_refused(fault):
    rows, points, case = fixture()
    after = initialize(points, case, 'all_methyl_zero')
    if fault == 'heavy':
        after[1][0] += 1e-10
    elif fault == 'radius':
        after[3] = [2*x for x in after[3]]
    elif fault == 'phase':
        after[3][1], after[3][2] = after[3][2], -after[3][1]
    elif fault == 'slots':
        after[4], after[5] = after[5], after[4]
    else:
        after.pop()
    with pytest.raises(ValueError):
        check_factorial_initial(rows, points, after, case, 'all_methyl_zero')


def test_inherited_refusals_and_partial_execution_do_not_shrink_denominator():
    common = dict(id='source', seed=1, conditioning='prepared-view', context_kind='prepared-protein-view')
    completed = dict(common, arm='unchanged', status='completed', exact_original_nm_endpoint=True,
        nodes=[dict(HA=dict(all_acute=[])), dict(HA=dict(all_acute=[]))])
    rows = [completed, *[dict(common, arm=arm, status='refused_source')
        for arm in ('all_methyl_zero', 'all_methylene', 'combined_zero')]]
    before = copy.deepcopy(rows)
    result = summarize(rows)
    assert result['enrolled'] == 4 and result['statuses'] == {'completed': 1, 'refused_source': 3}
    assert result['complete_matched_contexts'] == 0
    assert result['strata']['prepared-protein-view']['arms']['unchanged']['matched_contexts'] == 0
    assert rows == before


def test_known_replay_requires_original_published_run_binding(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    root = tmp_path/'tmp/pdbfixer-research/graph-methyl-control-run-01'
    folder = root/'synthetic/seed-1/graph_tetrahedral_zero'
    folder.mkdir(parents=True)
    raw = gzip.compress(json.dumps(dict(initial_nm=[[0., 0., 0.]], final_nm=[[1., 0., 0.]])).encode())
    (folder/'path.json.gz').write_bytes(raw)
    manifest = dict(status='completed', artifacts={'path.json.gz': sha(raw)})
    (folder/'result.json').write_text(json.dumps(manifest))
    run_raw = json.dumps(dict(cases=[manifest])).encode()
    (root/'evidence.json').write_bytes(run_raw)
    public = tmp_path/'docs/research/graph-methyl-control-evidence.json'
    public.parent.mkdir(parents=True)
    public.write_text(json.dumps(dict(run_sha256=sha(run_raw))))
    case = dict(id='synthetic', seed=1, context_kind='prepared-protein-view')
    result = known_replay(case, 'all_methyl_zero')
    assert result['path_sha256'] == sha(raw)
    assert result['data']['final_nm'] == [[1., 0., 0.]]
    (root/'evidence.json').write_text(json.dumps(dict(cases=[manifest], changed=True)))
    with pytest.raises(ValueError, match='published run binding'):
        known_replay(case, 'all_methyl_zero')
