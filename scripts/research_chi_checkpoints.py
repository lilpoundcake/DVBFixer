"""Native derivative-check checkpoints on ALL original mode-C starts.

Use supported AtomTreeMinimizer diagnostics, not Python C++ virtual callbacks.
Keep numerical settings and all evidence; any changed endpoint is NOT an exact
mode-C return trace. No native accuracy, selected outputs or production changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from unittest.mock import patch

from scripts import research_generated_chi_relax as parent
from scripts.research_chi_stopping_controls import check_settings
from scripts.research_dunbrack_prior import sha, write

BASE=Path('tmp/pdbfixer-research/chi-stopping-plans-01/absolute-1e-6-gmax1e-4.json')
BASE_RUN=Path('tmp/pdbfixer-research/chi-stopping-absolute-1e-6-gmax1e-4-run-01')


def norm_bounds(norm,dimensions):
    """Only bounds on max-gradient follow from an L2 norm, not a full vector."""
    if dimensions<1 or not math.isfinite(norm) or norm<0:
        raise ValueError('finite nonnegative norm and positive dimension required')
    return dict(maximum_component_lower_bound=norm/math.sqrt(dimensions),maximum_component_upper_bound=norm)


def plan(output):
    raw=BASE.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('bound mode-C source changed')
    for p in (BASE,BASE_RUN/'evidence.json',BASE_RUN/'executed-plan.json',
              Path('docs/research/chi-derivatives-evidence.json'),
              Path('tmp/pdbfixer-research/chi-calltrace-synthetic-01.log'),
              Path('tmp/pdbfixer-research/chi-calltrace-synthetic-02.log')):
        record['inputs'][str(p)]=sha(p.read_bytes())
    for root in ('atomtree-source-1a56185','atomtree-derivative-source-1a56185','minmover-source-1a56185'):
        for p in sorted((Path('tmp/pdbfixer-research')/root).iterdir()):
            record['inputs'][str(p)]=sha(p.read_bytes())
    for case in record['cases']:
        for cell in case['cells']:
            folder=BASE_RUN/case['id']/cell['declaration']['id']
            for name in ('before-atoms.json','after-atoms.json','result.json'):
                record['inputs'][str(folder/name)]=sha((folder/name).read_bytes())
    record['source_hashes']['scripts/research_chi_checkpoints.py']=sha(Path(__file__).read_bytes())
    record['schema']='dvbfixer.chi-checkpoints-plan.v1'
    record['checkpoint_protocol']=dict(parent_run=str(BASE_RUN),
        enrollment='ALL42 original mode-C starts/all15 pools/seven inherited nonexecutions; no final-gradient/native/geometry filtering',
        api='same target-only MoveMap and ScoreFunction; native AtomTreeMinimizer.run with same MinimizerOptions plus deriv_check=True, deriv_check_to_stdout=False; same pre-score as MinMover.inner_run_minimizer',
        native_steps_degrees=[.001,.002,.004,.008,.016],
        interpretation='native derivative-evaluation checkpoint norms, NOT full iteration/function/return-branch trace; L2 norm bounds max-component only',
        sparse_data_policy='DOF_DataPoint rows are populated by native checker only above its fixed0.001 analytic-or-numeric threshold; invalid default rows are ABSENT, never zero-gradient evidence',
        comparison='compare all before/after/result/initial-snapshot bytes to modeC; retain changed states and block exact-replay interpretation, never restore coordinates or relax comparison',
        report='all checkpoint norms and sparse raw diagnostics retained; no native metrics, profile selector or new scientific output')
    write(output,record)


def serialize_checkpoints(result):
    if result is None or result.n_deriv_check_results()<1:
        raise ValueError('native derivative checkpoints missing')
    rows=[]
    for i in range(1,result.n_deriv_check_results()+1):
        checkpoint=result.deriv_check_result(i)
        n=checkpoint.nangles()
        if n<1 or checkpoint.nsteps()!=5:
            raise ValueError('native diagnostic vector/step inventory differs')
        norm=float(checkpoint.best_norm_analytic())
        if not math.isfinite(norm) or norm<0:
            raise ValueError('nonfinite native analytic norm')
        points=[]
        for dof in range(1,n+1):
            for step in range(1,6):
                point=checkpoint.dof_step_data(dof,step)
                identifier=point.dof_id()
                if not identifier.valid():
                    continue
                points.append(dict(vector_index=dof,step_index=step,
                    dof_id=dict(residue=identifier.rsd(),atom=identifier.atomno(),type=str(identifier.type())),
                    analytic_per_degree=float(point.ana_deriv()),numeric_per_degree=float(point.num_deriv()),
                    dof_value=float(point.dof_val()),f00=float(point.f00()),f11=float(point.f11()),f22=float(point.f22())))
        rows.append(dict(checkpoint=i,dimensions=n,steps=5,analytic_l2_norm_per_degree=norm,
            max_component_bounds_per_degree=norm_bounds(norm,n),sparse_populated_points=points,
            minimum_absolute_derivative_deviation_per_dof=[float(checkpoint.abs_deriv_dev(j)) for j in range(1,n+1)]))
    return rows


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    if record['schema']!='dvbfixer.chi-checkpoints-plan.v1':
        raise ValueError('separate native checkpoint plan required')
    import pyrosetta.rosetta.protocols.minimization_packing as packing
    from pyrosetta.rosetta.core.optimization import AtomTreeMinimizer

    original=packing.MinMover
    inventory=[(case,cell) for case in record['cases'] for cell in case['cells']]
    created=[]

    class CheckpointMover:
        def __init__(self):
            self.inner=original()
            if self.inner.min_options().gmax_cutoff_for_convergence()!=1.:
                raise ValueError('default gradient cutoff differs')
            self.inner.min_options().gmax_cutoff_for_convergence(record['minimizer']['gradient_cutoff'])
            self.case,self.cell=inventory[len(created)]
            created.append(self)

        def __getattr__(self,name):
            return getattr(self.inner,name)

        def apply(self,pose):
            settings=dict(type=self.inner.min_type(),tolerance=self.inner.tolerance(),
                gradient_cutoff=self.inner.min_options().gmax_cutoff_for_convergence(),
                max_iter=self.inner.min_options().max_iter(),cartesian=self.inner.cartesian(),neighbor_list=self.inner.nb_list())
            check_settings(settings,record['minimizer'])
            if self.inner.abs_score_convergence_threshold()!=0.:
                raise ValueError('unexpected repeat-minimization policy')
            options=self.inner.min_options()
            options.deriv_check(True)
            options.deriv_check_to_stdout(False)
            movemap=self.inner.movemap(pose).clone()
            sf=self.inner.score_function()
            sf(pose)  # matches MinMover::inner_run_minimizer's pre-score
            native=AtomTreeMinimizer()
            native.run(pose,movemap,sf,options)
            checkpoints=serialize_checkpoints(native.deriv_check_result())
            folder=output/self.case['id']/self.cell['declaration']['id']
            write(folder/'checkpoints.json',dict(checkpoints=checkpoints,settings=settings,
                deriv_check=True,deriv_check_to_stdout=False,python_virtual_callback=False))

    with patch.object(packing,'MinMover',CheckpointMover):
        parent.run(plan_path,output)
    if packing.MinMover is not original or len(created)!=len(inventory):
        raise ValueError('checkpoint adapter not restored/inventory differs')
    root=Path(record['checkpoint_protocol']['parent_run'])
    comparisons=[]
    for case,cell in inventory:
        old=root/case['id']/cell['declaration']['id']
        new=output/case['id']/cell['declaration']['id']
        exact={name:(old/name).exists() and (new/name).exists() and (old/name).read_bytes()==(new/name).read_bytes()
               for name in ('before-atoms.json','after-atoms.json','result.json','pose-snapshot.json.gz')}
        checkpoint_path=new/'checkpoints.json'
        comparisons.append(dict(id=case['id'],index=cell['declaration']['index'],exact=exact,
            endpoint_exact=all(exact.values()),checkpoint_sha256=sha(checkpoint_path.read_bytes()) if checkpoint_path.exists() else None))
    write(output/'checkpoint-evidence.json',dict(schema='dvbfixer.chi-checkpoints-run.v1',comparisons=comparisons,
        plan_sha256=sha(plan_path.read_bytes()),hook_restored=True,native_coordinates_read=False,weights_modified=False))
    print(json.dumps(dict(starts=len(comparisons),exact_endpoint_repeats=sum(c['endpoint_exact'] for c in comparisons),
        changed_or_missing_endpoints=sum(not c['endpoint_exact'] for c in comparisons)),indent=2))


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
