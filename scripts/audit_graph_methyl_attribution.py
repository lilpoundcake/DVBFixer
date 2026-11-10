"""Saved-only independent term inventories, force sums and callback lineage."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def check_vector(actual, expected):
    if len(actual) != 3 or len(expected) != 3 or any(not math.isfinite(x) for x in (*actual, *expected)) or any(not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-4) for a, b in zip(actual, expected)):
        raise ValueError('attribution native vector closure differs')
    return max(abs(a-b) for a, b in zip(actual, expected))


def audit(plan_path, root, output):
    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw or [c['case'] for c in run['cases']] != plan['cases']:
        raise ValueError('complete saved attribution enrollment differs')
    counts, native_errors, gradient_errors = 0, [], []
    reports = []
    for cell in run['cases']:
        case = cell['case']
        source = Path(case['source'])
        rows = json.loads((source/'H-ownership.json').read_bytes())['rows']
        ha = next(i for i, r in enumerate(rows) if r['key'] == plan['target'])
        path = json.loads(gzip.decompress(Path(case['path']).read_bytes()))
        expected_stages = ['initial', *[f"callback-{i}" for i in range(len(path['callbacks']))], 'final']
        if [n['stage'] for n in cell['nodes']] != expected_stages:
            raise ValueError('attribution node/path lineage differs')
        pair_path = root/cell['raw_artifact']
        if sha(pair_path.read_bytes()) != cell['raw_sha256']:
            raise ValueError('retained raw pair artifact changed')
        n = 0
        with gzip.open(pair_path, 'rt') as stream:
            for n, line in enumerate(stream, 1):
                data = json.loads(line)
                summary = cell['nodes'][n-1]
                if data['summary'] != summary:
                    raise ValueError('raw/public pair summary differs')
                terms = data['pairs']
                if [t['index'] for t in terms] != [i for i in range(len(rows)) if i != ha]:
                    raise ValueError('complete HA partner inventory differs')
                if any(t['partner'] != rows[t['index']]['key'] or t['element'] != rows[t['index']]['element'] for t in terms):
                    raise ValueError('pair partner identity differs')
                summed = [sum(t['force_kj_mol_nm'][j] for t in terms) for j in range(3)]
                check_vector(summed, summary['total_nonbonded_force_kj_mol_nm'])
                nb = [c for c in summary['force_components'] if c['type'] == 'NonbondedForce']
                if len(nb) != 1:
                    raise ValueError('nonunique native nonbonded force')
                native_errors.append(check_vector(summed, nb[0]['force_kj_mol_nm']))
                top = sorted((t for t in terms if t['outward_projection_kj_mol_nm'] < 0), key=lambda t: (t['outward_projection_kj_mol_nm'], t['index']))[:3]
                if top != summary['largest_inward_pairs']:
                    raise ValueError('declared exact top-three inward reporting differs')
                if summary['stage'].startswith('callback-'):
                    index = int(summary['stage'].split('-')[1])
                    gradient = path['callbacks'][index]['gradient_kj_mol_nm'][3*ha:3*ha+3]
                    total = [sum(c['force_kj_mol_nm'][j] for c in summary['force_components']) for j in range(3)]
                    gradient_errors.append(check_vector(total, [-x for x in gradient]))
                counts += len(terms)
        if n != len(cell['nodes']):
            raise ValueError('raw pair nodes omitted')
        reports.append(cell)
    if counts != run['summary']['HA_pair_observations']:
        raise ValueError('raw all-pair count differs')
    write(output, dict(schema='dvbfixer.graph-methyl-ASP49-attribution-audit.v1', cases=reports,
        summary=dict(contexts=len(reports), nodes=sum(len(c['nodes']) for c in reports), HA_pair_observations=counts,
            maximum_native_pair_closure_error_kj_mol_nm=max(native_errors, default=None),
            maximum_original_callback_gradient_error_kj_mol_nm=max(gradient_errors, default=None)),
        plan_sha256=sha(raw), run_sha256=sha(run_raw), policy=plan['policy'],
        audit_script_sha256=sha(Path(__file__).read_bytes()), production_changed=False,
        minimization_executed=False, scientific_qualification=False, original_thr2_causality_resolved=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
