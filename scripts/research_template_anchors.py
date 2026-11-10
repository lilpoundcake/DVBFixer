"""Separate frozen anchor-set ablation reusing the unchanged stage tracer.

The tracer still records all anchors supplied by upstream. A separate transport
record explicitly identifies the subset actually used in the proper fit. Only a
process-local research-function hook changes; installed PDBFixer is never edited.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pdbfixer.pdbfixer as pf
from scripts import research_template_causality as tracer
from scripts.audit_template_causality import read_pdb
from scripts.branched_template_baseline import preflight
from scripts.research_rosetta_joint import sha
from scripts.research_rosetta_repack import write_json

PARENT_PLAN = Path('tmp/pdbfixer-research/causality-plan-01.json')
PARENT_RUN = Path('tmp/pdbfixer-research/causality-run-01')


def make_anchor_overlay(names, expected_reference, expected_template, policy, proper_fit, calls):
    """Names and coordinates are independently source/resource bound, not guessed."""
    if len(names) != 4 or set(names) != {'N', 'CA', 'C', 'O'} or policy not in {'all', 'ncac'}:
        raise ValueError('exact four supplied backbone anchors and declared policy required')
    indices = list(range(4)) if policy == 'all' else [i for i, name in enumerate(names) if name != 'O']

    def fit(points1, points2):
        p, q = np.asarray(points1, dtype=float), np.asarray(points2, dtype=float)
        if (p.shape != (4, 3) or q.shape != (4, 3) or not np.allclose(p, expected_reference, rtol=0., atol=1e-12)
                or not np.allclose(q, expected_template, rtol=0., atol=1e-12)):
            raise ValueError('unexpected overlay anchors/order; no silent fallback')
        used_p, used_q = p[indices], q[indices]
        result = proper_fit(used_p, used_q)
        calls.append(dict(supplied_names=list(names), used_names=[names[i] for i in indices], used_indices=indices,
                          supplied_reference_points_nm=p.tolist(), supplied_template_points_nm=q.tolist(),
                          used_reference_points_nm=used_p.tolist(), used_template_points_nm=used_q.tolist(),
                          rotation=np.asarray(result[1]).tolist(), translate_template_nm=np.asarray(result[0]).tolist(),
                          translate_reference_nm=np.asarray(result[2]).tolist()))
        return result

    return fit


def anchor_data(case, environment):
    target = case['target']
    key = target['chain_id'], target['sequence_number'], target['insertion_code']
    observed = {tuple(a['key'])[3]: a['xyz'] for a in read_pdb(Path(case['input_path']).read_bytes()) if tuple(a['key'])[:3] == key}
    filenames = [name for name in environment['resource_sha256'] if Path(name).stem.upper() == case['residue_name']]
    if len(filenames) != 1 or set(observed) != {'N', 'CA', 'C', 'O'}:
        raise ValueError('control requires exactly observed N/CA/C/O at target')
    path = Path(pf.__file__).parent/filenames[0]
    raw = path.read_bytes()
    if sha(raw) != environment['resource_sha256'][filenames[0]]:
        raise ValueError('canonical template resource changed')
    atoms = [a for a in read_pdb(raw) if a['key'][3] in observed]
    names = [a['key'][3] for a in atoms]
    return dict(names=names, reference_points_nm=[[x/10. for x in observed[n]] for n in names],
                template_points_nm=[[x/10. for x in a['xyz']] for a in atoms], resource=filenames[0], resource_sha256=sha(raw))


def plan(output):
    parent_raw = PARENT_PLAN.read_bytes()
    parent = json.loads(parent_raw)
    for filename, digest in parent['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('previous tracer implementation changed')
    if preflight() != parent['environment']:
        raise ValueError('scientific environment differs from parent')
    parent_run_raw = (PARENT_RUN/'evidence.json').read_bytes()
    if json.loads(parent_run_raw)['plan_sha256'] != sha(parent_raw):
        raise ValueError('parent run binding differs')
    cases = []
    for original in parent['cases']:
        case = dict(original)
        if sha(Path(case['input_path']).read_bytes()) != case['input_sha256']:
            raise ValueError('original masked input changed')
        case['anchor_data'] = anchor_data(case, parent['environment'])
        cases.append(case)
    hashes = dict(parent['source_hashes'])
    for filename in ('scripts/research_template_anchors.py', 'scripts/audit_template_causality.py',
                     'scripts/compare_branched_completion.py', 'src/dvbfixer/pdb_component_reconstruction.py',
                     'src/dvbfixer/domain/structure_identity.py', 'src/dvbfixer/diagnose/geometry.py'):
        hashes[filename] = sha(Path(filename).read_bytes())
    for path in sorted(Path('src/dvbfixer/protein_completion').glob('*.py')):
        hashes[str(path)] = sha(path.read_bytes())
    record = dict(schema='dvbfixer.template-anchor-plan.v1', cases=cases, environment=parent['environment'],
                  parent_plan_path=str(PARENT_PLAN), parent_plan_sha256=sha(parent_raw),
                  parent_run_path=str(PARENT_RUN), parent_run_sha256=sha(parent_run_raw),
                  arms=[dict(anchor_set=anchors, relaxation=relaxation) for anchors in ('all', 'ncac') for relaxation in ('normal', 'none')],
                  source_hashes=hashes, installed_pdbfixer_sha256=parent['installed_pdbfixer_sha256'],
                  intervention='proper SVD fit with original four supplied anchors or exactly N/CA/C subset, preserving upstream order; observed O stays fixed but does not steer subset fit',
                  trace_transport='unchanged tracer records four SUPPLIED anchors; separate anchor-policy.json records actual USED anchors and transform for every construction',
                  controls='all-anchor arms must reproduce corresponding saved parent proper arms exactly, including complete trace/iterations/coordinates',
                  geometry='unchanged comparator raw/rounded/global/target stereo and steric thresholds; no relaxation suppression interpreted as valid output',
                  scope='retrospective four-case anchor mechanism, not independent-family accuracy or a public repair',
                  public_capability=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('xb') as stream:
        stream.write((json.dumps(record, sort_keys=True, indent=2, allow_nan=False)+'\n').encode())


def run(plan_path, output):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    for filename, digest in record['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('frozen anchor experiment implementation changed')
    if preflight() != record['environment'] or sha(Path(pf.__file__).read_bytes()) != record['installed_pdbfixer_sha256']:
        raise ValueError('installed source/resources changed')
    if sha(Path(record['parent_plan_path']).read_bytes()) != record['parent_plan_sha256'] or sha((Path(record['parent_run_path'])/'evidence.json').read_bytes()) != record['parent_run_sha256']:
        raise ValueError('parent experiment changed')
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results = {}
    for case in record['cases']:
        if sha(Path(case['input_path']).read_bytes()) != case['input_sha256'] or anchor_data(case, record['environment']) != case['anchor_data']:
            raise ValueError('source-bound anchor mapping changed')
        data = case['anchor_data']
        for arm in record['arms']:
            name = case['id']+'-'+arm['anchor_set']+'-'+arm['relaxation']
            folder = output/name
            original_proper = tracer.proper_overlay
            calls = []
            hook = make_anchor_overlay(data['names'], data['reference_points_nm'], data['template_points_nm'],
                                       arm['anchor_set'], original_proper, calls)
            try:
                with patch.object(tracer, 'proper_overlay', hook):
                    traced = tracer.run_arm(case, dict(overlay='proper', relaxation=arm['relaxation']), folder)
                if tracer.proper_overlay is not original_proper:
                    raise ValueError('research proper-fit hook was not restored')
                if len(calls) != len(traced['overlays']):
                    raise ValueError('used-anchor transport omitted a fit')
                reproduced = None
                if arm['anchor_set'] == 'all':
                    previous = Path(record['parent_run_path'])/(case['id']+'-proper-'+arm['relaxation'])
                    if (folder/'evidence.json').read_bytes() != (previous/'evidence.json').read_bytes():
                        raise ValueError('all-anchor complete parent trace does not reproduce')
                    for event in traced['events']:
                        if (folder/event['coordinates']).read_bytes() != (previous/event['coordinates']).read_bytes():
                            raise ValueError('all-anchor full stage payload does not reproduce')
                    reproduced = True
                transport = dict(schema='dvbfixer.template-anchor-transport.v1', case_id=case['id'], arm=arm,
                                 calls=calls, research_hook_restored=True, all_anchor_parent_trace_exact=reproduced)
                write_json(folder/'anchor-policy.json', transport)
            except Exception as error:
                write_json(output/'execution-error.json', dict(case=name, type=type(error).__name__, error=str(error)))
                raise
            results[name] = dict(case_id=case['id'], arm=arm, trace_sha256=sha((folder/'evidence.json').read_bytes()),
                                 transport_sha256=sha((folder/'anchor-policy.json').read_bytes()),
                                 seed_count=sum(e['stage'] == 'returned-seed' for e in traced['events']))
    write_json(output/'evidence.json', dict(schema='dvbfixer.template-anchor-run.v1', plan_sha256=sha(raw), results=results,
                                           installed_source_unchanged=sha(Path(pf.__file__).read_bytes()) == record['installed_pdbfixer_sha256']))
    print(json.dumps(dict(cases=len(record['cases']), arms=len(results), returned_seeds=sum(r['seed_count'] for r in results.values())), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.output)
    else:
        run(args.plan, args.output)
