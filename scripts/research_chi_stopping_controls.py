"""Predeclared numerical stopping-policy controls, never native-tuned profiles.

Reuse the immutable target-chi runner and the same complete starting states.
Three sequential controls isolate energy criterion, tolerance and gradient cutoff.
No cross-profile selector, changed objective/weights, new starts or production.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

from scripts import research_generated_chi_relax as parent
from scripts.research_dunbrack_prior import sha, write

BASE=Path('tmp/pdbfixer-research/generated-chi-relax-plan-01.json')


def profiles():
    return [dict(id='absolute-0p01-gmax1',type='lbfgs_armijo_nonmonotone_atol',tolerance=.01,gradient_cutoff=1.),
            dict(id='absolute-1e-6-gmax1',type='lbfgs_armijo_nonmonotone_atol',tolerance=1e-6,gradient_cutoff=1.),
            dict(id='absolute-1e-6-gmax1e-4',type='lbfgs_armijo_nonmonotone_atol',tolerance=1e-6,gradient_cutoff=1e-4)]


def check_settings(actual,expected):
    for key in ('type','tolerance','gradient_cutoff','max_iter','cartesian','neighbor_list'):
        if actual[key]!=expected[key]:
            raise ValueError('runtime minimizer setting differs: '+key)


def plan(output):
    raw=BASE.read_bytes()
    base=json.loads(raw)
    for p,d in {**base['inputs'],**base['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('original relaxation input changed')
    source=Path('tmp/pdbfixer-research/minimizer-source-1a56185/Minimizer.cc')
    provenance=source.with_name('provenance.json')
    origin=json.loads(provenance.read_bytes())
    if sha(source.read_bytes())!=origin['sha256'] or origin['source_revision']!=base['expected_revision']:
        raise ValueError('upstream minimizer source binding differs')
    output.mkdir(parents=True,exist_ok=False)
    plans=[]
    for profile in profiles():
        record=json.loads(raw)
        record['schema']='dvbfixer.generated-chi-relax-plan.v2-stopping-control'
        for p in (BASE,source,provenance,Path('tmp/pdbfixer-research/generated-chi-relax-run-01/evidence.json')):
            record['inputs'][str(p)]=sha(p.read_bytes())
        record['source_hashes']['scripts/research_chi_stopping_controls.py']=sha(Path(__file__).read_bytes())
        record['minimizer'].update({k:v for k,v in profile.items() if k!='id'})
        record['control_id']=profile['id']
        record['control_enrollment']='ALL 42 original feasible mean starts in all 15 pools; seven inherited nonexecutions retained; no native-error, movement or score-margin filter'
        record['control_interpretation']='sequential stopping-policy sensitivity only; no native reporting or choice of best control; max200 fixed; not a convergence certificate'
        path=output/(profile['id']+'.json')
        write(path,record)
        plans.append(dict(id=profile['id'],filename=path.name,sha256=sha(path.read_bytes())))
    write(output/'suite.json',dict(schema='dvbfixer.chi-stopping-suite.v1',profiles=plans,
        baseline_plan_sha256=sha(raw),native_accuracy_used=False,cross_profile_selection=False))


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    if record['schema']!='dvbfixer.generated-chi-relax-plan.v2-stopping-control':
        raise ValueError('separate stopping-control plan required')
    import pyrosetta.rosetta.protocols.minimization_packing as packing

    original=packing.MinMover
    created=[]

    def make():
        mover=original()
        if mover.min_options().gmax_cutoff_for_convergence()!=1.:
            raise ValueError('original default gradient cutoff differs')
        mover.min_options().gmax_cutoff_for_convergence(record['minimizer']['gradient_cutoff'])
        created.append(mover)
        return mover

    with patch.object(packing,'MinMover',make):
        parent.run(plan_path,output)
    if packing.MinMover is not original:
        raise ValueError('process-local stopping adapter not restored')
    if len(created)!=sum(len(c['cells']) for c in record['cases']):
        raise ValueError('expected starting-state minimizer inventory differs')
    actual=[]
    for mover in created:
        options=mover.min_options()
        row=dict(type=mover.min_type(),tolerance=mover.tolerance(),gradient_cutoff=options.gmax_cutoff_for_convergence(),
                 max_iter=options.max_iter(),cartesian=mover.cartesian(),neighbor_list=mover.nb_list())
        check_settings(row,record['minimizer'])
        actual.append(row)
    write(output/'stopping-control-settings.json',dict(control_id=record['control_id'],hook_restored=True,
        actual_minimizers=actual,native_coordinates_read=False,weights_modified=False,cross_profile_selection=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('run')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.output)
    else:
        run(args.plan,args.output)
