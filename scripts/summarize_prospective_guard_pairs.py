"""Predeclared ALL-arm-pair saved accuracy summaries, never a profile selector."""
from __future__ import annotations

import argparse
import itertools
import json
import math
import statistics
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write

METRICS=('exact_name_rmsd_angstrom','generated_only_val_symmetry_rmsd_angstrom')


def summarize(arms,cases):
    rows=[]
    for first,second in itertools.combinations(arms,2):
        names=[n for n,c in cases.items() if first in c['accuracies'] and second in c['accuracies']]
        result={}
        for arm in (first,second):
            result[arm]={}
            for metric in METRICS:
                values=[cases[n]['accuracies'][arm]['rounded'][metric] for n in names]
                if any(not math.isfinite(v) for v in values):
                    raise ValueError('finite matched accuracy values required')
                result[arm][metric]=statistics.median(values) if values else None
        rows.append(dict(first=first,second=second,matched_count=len(names),matched_mask_ids=names,medians=result))
    return rows


def plan(ledger_path,report_path,output):
    if report_path.exists():
        raise ValueError('pair-summary protocol must be frozen before native report exists')
    raw=ledger_path.read_bytes()
    ledger=json.loads(raw)
    if ledger['native_coordinates_read'] or ledger['profile_selected']:
        raise ValueError('native-free all-arm choice ledger required')
    write(output,dict(schema='dvbfixer.prospective-guard-pairs-plan.v1',ledger_path=str(ledger_path),
        ledger_sha256=sha(raw),expected_report_path=str(report_path),arm_order=ledger['arm_order'],
        source_sha256=sha(Path(__file__).read_bytes()),metrics=METRICS,
        enrollment='ALL unordered arm pairs, ALL shared selected outputs per pair; no accuracy/native-defect/margin filtering',
        interpretation='matched descriptive medians only, no preferred profile/significance/qualification claim'))


def audit(plan_path,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    if record['source_sha256']!=sha(Path(__file__).read_bytes()):
        raise ValueError('pre-report pair summary code changed')
    ledger_raw=Path(record['ledger_path']).read_bytes()
    if sha(ledger_raw)!=record['ledger_sha256']:
        raise ValueError('all-arm selections changed')
    report_raw=Path(record['expected_report_path']).read_bytes()
    report=json.loads(report_raw)
    if report['selection_ledger_sha256']!=sha(ledger_raw) or not report['native_reporting_after_all_choices']:
        raise ValueError('native report not bound to prior choices')
    rows=summarize(record['arm_order'],report['cases'])
    write(output,dict(schema='dvbfixer.prospective-guard-pairs-audit.v1',pairs=rows,
        pair_count=len(rows),plan_sha256=sha(raw),report_sha256=sha(report_raw),
        native_coordinates_read=False,native_accuracy_read=True,profile_selected=False,
        scientific_qualification=False,production_changed=False))
    print(json.dumps(dict(pairs=len(rows),matched_counts=[r['matched_count'] for r in rows]),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--ledger',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('audit')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.ledger,args.report,args.output)
    else:
        audit(args.plan,args.output)
