"""Logging-only replay of all mode-C starts; numerical settings stay frozen.

Capture native Minimizer diagnostics per apply, require exact saved endpoints.
No changed optimizer settings, native accuracy, profile choice or production use.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

from scripts import research_generated_chi_relax as parent
from scripts.research_chi_stopping_controls import check_settings
from scripts.research_dunbrack_prior import sha, write

BASE=Path('tmp/pdbfixer-research/chi-stopping-plans-01/absolute-1e-6-gmax1e-4.json')
BASE_RUN=Path('tmp/pdbfixer-research/chi-stopping-absolute-1e-6-gmax1e-4-run-01')


@contextmanager
def capture_native(path):
    """Local child-process fd capture, restored even when apply fails."""
    sys.stdout.flush()
    sys.stderr.flush()
    stdout,stderr=os.dup(1),os.dup(2)
    try:
        with path.open('xb') as capture:
            os.dup2(capture.fileno(),1)
            os.dup2(capture.fileno(),2)
            try:
                yield
            finally:
                sys.stdout.flush()
                sys.stderr.flush()
    finally:
        os.dup2(stdout,1)
        os.dup2(stderr,2)
        os.close(stdout)
        os.close(stderr)


def classify(text):
    maximum='LBFGS MAX CYCLES 200 EXCEEDED, BUT FUNC NOT CONVERGED!'
    failed='Line search failed even after resetting Hessian; aborting at iter#'
    return dict(explicit_max200_warning=maximum in text,explicit_line_search_abort=failed in text,
                termination_certified=False)


def plan(output):
    raw=BASE.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('mode-C state changed')
    for p in (BASE,BASE_RUN/'executed-plan.json',BASE_RUN/'evidence.json'):
        record['inputs'][str(p)]=sha(p.read_bytes())
    run=json.loads((BASE_RUN/'evidence.json').read_bytes())
    for case in record['cases']:
        for cell in case['cells']:
            folder=BASE_RUN/case['id']/cell['declaration']['id']
            for filename in ('before-atoms.json','after-atoms.json','result.json'):
                record['inputs'][str(folder/filename)]=sha((folder/filename).read_bytes())
    record['schema']='dvbfixer.chi-termination-plan.v1'
    record['source_hashes']['scripts/research_chi_termination.py']=sha(Path(__file__).read_bytes())
    record['init_flags']+=' -unmute core.optimization.Minimizer'
    record['termination_replay']=dict(parent_run=str(BASE_RUN),parent_plan_sha256=sha(raw),
        logging_only=True,exact_endpoint_traces_and_results_required=True,
        enrollment='ALL42 original mode-C starts, same15 pools/seven inherited nonexecutions, no finite-descent or geometry filtering',
        interpretation='explicit max-cycle or line-search messages only; missing diagnostics do not certify normal convergence')
    if run['plan_sha256']!=sha(raw):
        raise ValueError('mode-C execution binding differs')
    write(output,record)


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    if record['schema']!='dvbfixer.chi-termination-plan.v1':
        raise ValueError('logging-only frozen plan required')
    import pyrosetta.rosetta.protocols.minimization_packing as packing

    original=packing.MinMover
    inventory=[(c,cell) for c in record['cases'] for cell in c['cells']]
    traces=[]

    class LoggedMover:
        def __init__(self):
            self.inner=original()
            if self.inner.min_options().gmax_cutoff_for_convergence()!=1.:
                raise ValueError('unexpected default gradient cutoff')
            self.inner.min_options().gmax_cutoff_for_convergence(record['minimizer']['gradient_cutoff'])

        def __getattr__(self,name):
            return getattr(self.inner,name)

        def apply(self,pose):
            index=len(traces)
            case,cell=inventory[index]
            actual=dict(type=self.inner.min_type(),tolerance=self.inner.tolerance(),
                gradient_cutoff=self.inner.min_options().gmax_cutoff_for_convergence(),
                max_iter=self.inner.min_options().max_iter(),cartesian=self.inner.cartesian(),neighbor_list=self.inner.nb_list())
            check_settings(actual,record['minimizer'])
            if self.inner.min_options().silent():
                raise ValueError('minimizer silent flag hides required diagnostics')
            path=output/case['id']/cell['declaration']['id']/'minimizer.log'
            with capture_native(path):
                self.inner.apply(pose)
            text=path.read_text()
            traces.append(dict(id=case['id'],index=cell['declaration']['index'],path=str(path),sha256=sha(path.read_bytes()),
                actual_settings=actual,**classify(text)))

    with patch.object(packing,'MinMover',LoggedMover):
        parent.run(plan_path,output)
    if packing.MinMover is not original or len(traces)!=len(inventory):
        raise ValueError('logging adapter not restored or complete inventory differs')
    comparisons=[]
    root=Path(record['termination_replay']['parent_run'])
    for (case,cell),trace in zip(inventory,traces,strict=True):
        old=root/case['id']/cell['declaration']['id']
        new=output/case['id']/cell['declaration']['id']
        exact={name:(old/name).read_bytes()==(new/name).read_bytes() for name in ('before-atoms.json','after-atoms.json','result.json','pose-snapshot.json.gz')}
        if not all(exact.values()):
            write(output/'replay-error.json',dict(id=case['id'],index=cell['declaration']['index'],exact=exact))
            raise ValueError('logging-only replay changed a saved endpoint/result')
        comparisons.append(dict(id=case['id'],index=cell['declaration']['index'],exact=exact))
    write(output/'termination-evidence.json',dict(schema='dvbfixer.chi-termination-run.v1',traces=traces,
        comparisons=comparisons,hook_restored=True,native_coordinates_read=False,weights_modified=False,
        plan_sha256=sha(plan_path.read_bytes())))
    print(json.dumps(dict(endpoints=len(traces),exact_endpoint_repeats=len(comparisons),
        explicit_max200=sum(t['explicit_max200_warning'] for t in traces),
        explicit_line_search_abort=sum(t['explicit_line_search_abort'] for t in traces)),indent=2))


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
