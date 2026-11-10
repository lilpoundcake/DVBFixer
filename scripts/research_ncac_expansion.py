"""Frozen all-whole-mask NCAC evaluation, not a production atom-completion route.

Keep four seen mechanism cases separate from seven cases unexecuted under NCAC.
No accuracy/reference data enters reconstruction, anchor choice or seed selection.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import pdbfixer.pdbfixer as pf
from scripts import research_template_causality as tracer
from scripts.branched_template_baseline import preflight
from scripts.research_rosetta_joint import sha, verify
from scripts.research_rosetta_repack import write_json
from scripts.research_template_anchors import anchor_data, make_anchor_overlay

MASK_ROOT = Path('tmp/pdbfixer-research/independent-masks-01')
COMPARISON_ROOT = Path('tmp/pdbfixer-research/independent-comparison-01')
ANCHOR_PLAN = Path('tmp/pdbfixer-research/anchor-plan-01.json')
ANCHOR_RUN = Path('tmp/pdbfixer-research/anchor-run-01')
PRIOR_AUDIT = Path('docs/research/rosetta-independent-controls-evidence.json')
PRIOR_SCORE_PLAN = Path('tmp/pdbfixer-research/independent-score-plan-01/plan.json')
PRIOR_SCORE_RUN = Path('tmp/pdbfixer-research/independent-score-run-01/evidence.json')


def whole_cases(cases, seen_ids):
    selected = [c for c in cases if 'CB' in c['mask'] and len(c['mask']) == 3]
    ids = [c['id'] for c in selected]
    if len(ids) != len(set(ids)) or not set(seen_ids) <= set(ids):
        raise ValueError('whole-mask inventory duplicates or omits previously tested case')
    return [dict(c, ncac_previously_executed=c['id'] in seen_ids) for c in selected]


def plan(output):
    mask_raw = (MASK_ROOT/'plan.json').read_bytes()
    mask_plan = json.loads(mask_raw)
    verify(mask_plan, MASK_ROOT)
    anchor_raw = ANCHOR_PLAN.read_bytes()
    previous = json.loads(anchor_raw)
    environment = preflight()
    if environment != previous['environment']:
        raise ValueError('scientific environment differs from previous NCAC controls')
    for filename, digest in previous['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('previous NCAC implementation changed')
    selected = whole_cases(mask_plan['cases'], [c['id'] for c in previous['cases']])
    if len(selected) != 11 or sum(c['ncac_previously_executed'] for c in selected) != 4:
        raise ValueError('complete eleven-mask/four-seen/seven-unseen-variant cohort required')
    cases = []
    for item in selected:
        case = dict(item, input_path=str(COMPARISON_ROOT/item['id']/'input.pdb'),
                    input_sha256=mask_plan['artifacts'][item['source']],
                    prior_proposal_path=str(COMPARISON_ROOT/item['id']/'template-proposal.json'))
        if sha(Path(case['input_path']).read_bytes()) != case['input_sha256']:
            raise ValueError('saved source differs from frozen mask')
        case['prior_proposal_sha256'] = sha(Path(case['prior_proposal_path']).read_bytes())
        case['anchor_data'] = anchor_data(case, environment)
        cases.append(case)
    inputs = {str(MASK_ROOT/'plan.json'): sha(mask_raw), str(ANCHOR_PLAN): sha(anchor_raw),
              str(ANCHOR_RUN/'evidence.json'): sha((ANCHOR_RUN/'evidence.json').read_bytes()),
              str(COMPARISON_ROOT/'evidence.json'): sha((COMPARISON_ROOT/'evidence.json').read_bytes()),
              str(PRIOR_AUDIT): sha(PRIOR_AUDIT.read_bytes()), str(PRIOR_SCORE_PLAN): sha(PRIOR_SCORE_PLAN.read_bytes()),
              str(PRIOR_SCORE_RUN): sha(PRIOR_SCORE_RUN.read_bytes())}
    hashes = dict(previous['source_hashes'])
    for filename in ('scripts/research_ncac_expansion.py', 'scripts/audit_template_anchors.py',
                     'scripts/unfiltered_branched_comparison.py', 'scripts/summarize_rosetta_fixed_scoring.py'):
        hashes[filename] = sha(Path(filename).read_bytes())
    record = dict(schema='dvbfixer.ncac-expansion-plan.v1', cases=cases, methods=['stock', 'ncac'],
                  environment=environment, installed_pdbfixer_sha256=sha(Path(pf.__file__).read_bytes()),
                  inputs=inputs, source_hashes=hashes, mask_root=str(MASK_ROOT), comparison_root=str(COMPARISON_ROOT),
                  seen_anchor_run=str(ANCHOR_RUN), prior_audit=str(PRIOR_AUDIT), prior_score_plan=str(PRIOR_SCORE_PLAN),
                  prior_score_run=str(PRIOR_SCORE_RUN), source_case_count=len(mask_plan['cases']),
                  selection='ALL whole-sidechain masks in original enrolled order, never filtered by native accuracy/geometry or earlier template outcome',
                  variant='unchanged proper N/CA/C-only fit, normal upstream relaxation and original helper 1..5; no new guard/reflection/scoring',
                  controls='every stock state must reproduce original untraced template seeds; four seen NCAC complete traces must reproduce previous NCAC-normal arm',
                  reporting='commit all proposals and geometry-only audits before reading native coordinates; raw and rounded exact/symmetry reporting; all-raw diagnostic and mutual-feasible accuracy denominators remain distinct',
                  scope='four retrospective mechanism cases plus seven cases new only to the NCAC variant; original native/baseline outcomes already known; not blinded/new-family qualification',
                  public_capability=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('xb') as stream:
        stream.write((json.dumps(record, sort_keys=True, indent=2, allow_nan=False)+'\n').encode())


def run(plan_path, output):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    for filename, digest in {**record['inputs'], **record['source_hashes']}.items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('frozen expansion input/implementation changed')
    if preflight() != record['environment']:
        raise ValueError('scientific environment changed')
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results = {}
    for case in record['cases']:
        if sha(Path(case['input_path']).read_bytes()) != case['input_sha256'] or sha(Path(case['prior_proposal_path']).read_bytes()) != case['prior_proposal_sha256']:
            raise ValueError('original source/template proposal changed')
        for method in record['methods']:
            name = case['id']+'-'+method
            root = output/name
            calls, reproduced = [], None
            try:
                if method == 'stock':
                    trace = tracer.run_arm(case, dict(overlay='stock', relaxation='normal'), root)
                else:
                    data = case['anchor_data']
                    if anchor_data(case, record['environment']) != data:
                        raise ValueError('source-bound NCAC anchor mapping changed')
                    original = tracer.proper_overlay
                    hook = make_anchor_overlay(data['names'], data['reference_points_nm'], data['template_points_nm'], 'ncac', original, calls)
                    with patch.object(tracer, 'proper_overlay', hook):
                        trace = tracer.run_arm(case, dict(overlay='proper', relaxation='normal'), root)
                    if tracer.proper_overlay is not original or len(calls) != len(trace['overlays']):
                        raise ValueError('NCAC hook restoration/used-anchor inventory differs')
                    if case['ncac_previously_executed']:
                        parent = Path(record['seen_anchor_run'])/(case['id']+'-ncac-normal')
                        if (root/'evidence.json').read_bytes() != (parent/'evidence.json').read_bytes():
                            raise ValueError('seen NCAC complete trace does not reproduce')
                        for event in trace['events']:
                            if (root/event['coordinates']).read_bytes() != (parent/event['coordinates']).read_bytes():
                                raise ValueError('seen NCAC stage bytes do not reproduce')
                        reproduced = True
                    write_json(root/'anchor-policy.json', dict(method=method, calls=calls, research_hook_restored=True,
                                                               seen_ncac_trace_exact=reproduced))
            except Exception as error:
                write_json(output/'execution-error.json', dict(case=name, type=type(error).__name__, error=str(error)))
                raise
            results[name] = dict(case_id=case['id'], method=method, ncac_previously_executed=case['ncac_previously_executed'],
                                 trace_sha256=sha((root/'evidence.json').read_bytes()),
                                 transport_sha256=sha((root/'anchor-policy.json').read_bytes()) if method == 'ncac' else None,
                                 seed_count=sum(e['stage'] == 'returned-seed' for e in trace['events']), seen_ncac_trace_exact=reproduced)
    write_json(output/'evidence.json', dict(schema='dvbfixer.ncac-expansion-run.v1', plan_sha256=sha(raw), results=results,
                                           installed_source_unchanged=sha(Path(pf.__file__).read_bytes()) == record['installed_pdbfixer_sha256']))
    print(json.dumps(dict(cases=len(record['cases']), method_arms=len(results), seed_states=sum(r['seed_count'] for r in results.values())), indent=2))


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
