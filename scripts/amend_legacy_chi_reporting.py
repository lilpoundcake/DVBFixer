"""Explicit saved-report correction: chi must use actual relaxed backbone.

The shared fixed-observation reporter's RMSD is valid for legacy outputs, but its
chi helper combines new additions with original anchors. Correct only that field;
retain v1 report/source, no resampling, scientific selection or RMSD modification.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.audit_legacy_comparator import atoms
from scripts.research_dunbrack_prior import sha, write

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef


def actual_chi(points,residue_name):
    tip='OG1' if residue_name=='THR' else 'CG1'
    chi=signed_dihedral_degrees(*(points[n] for n in ('N','CA','CB',tip)))
    if chi is None:
        raise ValueError('actual-frame chi undefined')
    return chi


def plan(legacy_plan,run_root,old_report,output):
    record=json.loads(legacy_plan.read_bytes())
    paths=[legacy_plan,old_report,Path('scripts/audit_legacy_comparator.py'),
           Path('scripts/amend_legacy_chi_reporting.py')]
    for case in record['cases']:
        paths.append(Path(record['mask_root'])/case['reference'])
        for name in ('prepared.pdb','minimized.pdb'):
            path=run_root/case['id']/name
            if path.exists():
                paths.append(path)
    write(output,dict(schema='dvbfixer.legacy-chi-amendment-plan.v1',legacy_plan=str(legacy_plan),
        old_report=str(old_report),run_root=str(run_root),inputs={str(p):sha(p.read_bytes()) for p in paths},
        correction='only named_chi1_error_degrees uses actual full output N/CA/CB/tip; old fixed-observation helper used hybrid anchors; exact/symmetry RMSD and all choices unchanged',
        native_accuracy_already_seen=True,resampling=False,profile_selected=False))


def run(plan_path,output):
    amendment=json.loads(plan_path.read_bytes())
    for p,d in amendment['inputs'].items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('bound chi-report amendment input changed')
    legacy=json.loads(Path(amendment['legacy_plan']).read_bytes())
    old_raw=Path(amendment['old_report']).read_bytes()
    record=json.loads(old_raw)
    before=json.loads(old_raw)
    for case,result in zip(legacy['cases'],record['cases'],strict=True):
        target=ComponentInstanceRef(**case['target'])
        native=Path(legacy['mask_root'])/case['reference']
        truth={n:a.position for (component,n),a in atoms(native.read_bytes()).items() if component==target}
        for stage,filename in (('prepare','prepared.pdb'),('minimize','minimized.pdb')):
            if stage not in result['stages']:
                continue
            actual=Path(amendment['run_root'])/case['id']/filename
            points={n:a.position for (component,n),a in atoms(actual.read_bytes()).items() if component==target}
            residue_name=next(a.residue_name for (component,_),a in atoms(actual.read_bytes()).items() if component==target)
            chi,ref=actual_chi(points,residue_name),actual_chi(truth,residue_name)
            result['stages'][stage]['original_frame_accuracy']['named_chi1_error_degrees']=abs((chi-ref+180)%360-180)
    # The old native RMSDs/coverage/geometry/movement are not adapted to the result.
    check=json.loads(json.dumps(record))
    for old_case,new_case in zip(before['cases'],check['cases'],strict=True):
        for stage in old_case['stages']:
            new_case['stages'][stage]['original_frame_accuracy']['named_chi1_error_degrees']=old_case['stages'][stage]['original_frame_accuracy']['named_chi1_error_degrees']
    if check!=before:
        raise ValueError('chi amendment modified another scientific result')
    record.update(schema='dvbfixer.legacy-comparator-report.v2',
        chi_reporting_amendment_plan_sha256=sha(plan_path.read_bytes()),superseded_report_sha256=sha(old_raw),
        chi_scope='actual full-output backbone/tip; no hybrid original anchors',rmsd_unchanged=True)
    write(output,record)
    print('Corrected actual-frame chi only; all RMSD, coverage, geometry and motion unchanged.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--legacy-plan',type=Path,required=True)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--old-report',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('run')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.legacy_plan,args.run,args.old_report,args.output)
    else:
        run(args.plan,args.output)
