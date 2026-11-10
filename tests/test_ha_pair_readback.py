"""Independent private-pair inventory/sum failures, no OpenMM execution."""
import copy
import gzip
import json
import math
from pathlib import Path

import pytest
from scripts.audit_ha_pair_readback import audit
from scripts.research_dunbrack_prior import sha


def fixture(tmp_path, fault=None):
    radius = math.sqrt(8/9)
    points = [[0., 0., 0.], [radius, 0., 1/3], [-radius/2, radius*math.sqrt(3)/2, 1/3],
        [-radius/2, -radius*math.sqrt(3)/2, 1/3], [0., 0., -.1]]
    rows = [dict(key=['A', '16', '', 'GLN', n], element=n[0], xyz=p)
            for n, p in zip(('CA', 'N', 'C', 'CB', 'HA'), points, strict=True)]
    source = tmp_path/'source'
    source.mkdir()
    (source/'h-internal-before.json.gz').write_bytes(gzip.compress(json.dumps(dict(atoms=rows)).encode(), mtime=0))
    case = dict(id='synthetic', seed=1, source_arm='raw', source=str(source))
    parent = tmp_path/'parent.json'
    parent.write_text(json.dumps(dict(cases=[case])))
    prior = tmp_path/'prior.json'
    prior.write_text(json.dumps(dict(cases=[dict(case, nodes=[dict(label='initial', ordinal=None,
        force_components=[dict(type='NonbondedForce', HA_force_kj_mol_nm=[0., 0., 0.])])])])))
    script = Path('scripts/audit_ha_nonbonded_pairs.py').resolve()
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={str(script): sha(script.read_bytes())},
        path_plan=str(parent), prior=str(prior), openmm_version='synthetic', policy={})))
    terms = [dict(partner_index=i, partner=r['key'], element=r['element'], force_kj_mol_nm=[0., 0., 0.],
        outward_projection_kj_mol_nm=0., LJ_outward_projection_kj_mol_nm=0., Coulomb_outward_projection_kj_mol_nm=0.) for i, r in enumerate(rows[:4])]
    categories = {n: dict(count=c, outward_projection_kj_mol_nm=0.) for n, c in
        [('H_partner', 0), ('heavy_partner', 4), ('same_residue', 4), ('other_residue', 0)]}
    summary = dict(label='initial', ordinal=None, pair_count=4, categories=categories,
        total_force_kj_mol_nm=[0., 0., 0.], largest_inward_axis_pairs=[])
    private = copy.deepcopy(summary)
    private['terms'] = terms
    if fault == 'duplicate':
        terms[-1] = terms[0]
    elif fault == 'force':
        terms[0]['force_kj_mol_nm'][0] = 1.
    elif fault == 'selected':
        private['largest_inward_axis_pairs'] = terms[:1]
        summary['largest_inward_axis_pairs'] = terms[:1]
    elif fault == 'category':
        private['categories']['heavy_partner']['count'] = 3
        summary['categories']['heavy_partner']['count'] = 3
    root = tmp_path/'run'
    root.mkdir()
    (root/'executed-plan.json').write_bytes(plan.read_bytes())
    raw_pairs = root/'pairs.jsonl.gz'
    raw_pairs.write_bytes(gzip.compress((json.dumps(private)+'\n').encode(), mtime=0))
    cell = dict(case, status='evaluated', nodes=[summary], raw_artifact=raw_pairs.name, raw_sha256=sha(raw_pairs.read_bytes()))
    (root/'evidence.json').write_text(json.dumps(dict(plan_sha256=sha(plan.read_bytes()), cases=[cell], summary=dict(HA_pair_observations=4))))
    return plan, root, tmp_path/'audit.json'


def test_complete_saved_private_inventory_can_be_audited(tmp_path):
    plan, root, output = fixture(tmp_path)
    audit(plan, root, output)
    assert json.loads(output.read_bytes())['summary']['independently_checked_pairs'] == 4


@pytest.mark.parametrize('fault', ['duplicate', 'force', 'selected', 'category'])
def test_independent_pair_audit_does_not_accept_partial_or_inconsistent_evidence(tmp_path, fault):
    plan, root, output = fixture(tmp_path, fault)
    with pytest.raises(ValueError):
        audit(plan, root, output)
    assert not output.exists()
