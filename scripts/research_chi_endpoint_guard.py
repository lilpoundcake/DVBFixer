"""Private bounded actual-endpoint-gradient guard on ALL original mode-C starts.

Same inner optimizer; at most three total cycles, checked at the returned pose.
No native accuracy, optimizer grid, proposal selection or production capability.
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
from scripts.research_generated_chi_relax import atom_state, movement

BASE=Path('tmp/pdbfixer-research/chi-stopping-plans-01/absolute-1e-6-gmax1e-4.json')
BASE_RUN=Path('tmp/pdbfixer-research/chi-stopping-absolute-1e-6-gmax1e-4-run-01')


def guard_decision(maximum,cycle,limit,cutoff):
    if not math.isfinite(maximum) or maximum<0 or not 1<=cycle<=limit:
        raise ValueError('finite gradient and bounded cycle required')
    return 'endpoint_gradient_met' if maximum<=cutoff else 'budget_exhausted' if cycle==limit else 'continue'


def plan(output):
    raw=BASE.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('original mode-C input changed')
    for p in (BASE,BASE_RUN/'evidence.json',BASE_RUN/'executed-plan.json',
              Path('docs/research/chi-checkpoint-components-evidence.json'),
              Path('tmp/pdbfixer-research/chi-derivatives-plan-01.json')):
        record['inputs'][str(p)]=sha(p.read_bytes())
    for case in record['cases']:
        for cell in case['cells']:
            folder=BASE_RUN/case['id']/cell['declaration']['id']
            for name in ('before-atoms.json','after-atoms.json','result.json'):
                record['inputs'][str(folder/name)]=sha((folder/name).read_bytes())
    for name in ('scripts/research_chi_endpoint_guard.py','scripts/research_chi_derivatives.py'):
        record['source_hashes'][name]=sha(Path(name).read_bytes())
    record['schema']='dvbfixer.chi-endpoint-guard-plan.v1'
    record['guard_protocol']=dict(parent_run=str(BASE_RUN),total_cycle_limit=3,
        cutoff_per_degree=record['minimizer']['gradient_cutoff'],
        center_score_absolute_bound=1e-8,mapped_endpoint_numerical_bound_angstrom=1e-9,
        enrollment='ALL42 original starts/all15 pools/seven inherited nonexecutions, no native/gradient/geometry filtering',
        algorithm='one unchanged mode-C inner MinMover cycle, actual returned-pose gradient on independent clone, stop if <=gmax else reenter identical inner optimizer; maximum three total cycles; no new candidate starts',
        first_cycle='full returned atom trace and raw score must reproduce original mode C EXACTLY; no fallback/coordinate restore',
        gradient='same target-only MinimizerMap and AtomTreeMultifunc.dfunc, degree-scaled variables, independent ScoreFunction clone; record every phase and all measured gradients',
        geometry='independently audit every retained cycle, not used to choose continuation; no candidate ranking/native reporting or geometry threshold changes',
        interpretation='bounded guarded-termination experiment, not guaranteed convergence/global minimum/accuracy or production qualification')
    write(output,record)


def endpoint_gradient(pose,movemap,sf,case,bound,score_bound):
    from pyrosetta.rosetta.core.id import CHI, PHI
    from pyrosetta.rosetta.core.optimization import AtomTreeMultifunc, MinimizerMap
    from pyrosetta.rosetta.utility import vector1_double

    observations=atom_state(pose,case)
    work=pose.clone()
    scorer=sf.clone()
    before_score=float(scorer(work))
    minmap=MinimizerMap()
    minmap.setup(work,movemap)
    scorer.setup_for_minimizing(work,minmap)
    func=AtomTreeMultifunc(work,minmap,scorer,False,False)
    x=vector1_double()
    gradient=vector1_double()
    for _ in range(minmap.nangles()):
        x.append(0.)
        gradient.append(0.)
    minmap.copy_dofs_from_pose(work,x)
    center=float(func(x))
    func.dfunc(x,gradient)
    nodes=list(minmap.dof_nodes())
    rows=[]
    info=work.pdb_info()
    for node,value in zip(nodes,gradient,strict=True):
        torsion=node.torsion_id()
        index=torsion.rsd()
        key=[info.chain(index),info.number(index),info.icode(index),work.residue(index).name3()]
        if node.type()!=PHI or torsion.type()!=CHI or key!=case['target_key'] or not math.isclose(minmap.torsion_scale_factor(node),180/math.pi,rel_tol=0.,abs_tol=1e-12) or not math.isfinite(value):
            raise ValueError('actual endpoint derivative identity/scale/finite check failed')
        rows.append(dict(chi=torsion.torsion(),analytic_per_degree=float(value)))
    if not rows or sorted(r['chi'] for r in rows)!=list(range(1,len(rows)+1)):
        raise ValueError('actual endpoint gradient vector incomplete')
    moved=movement(observations,atom_state(work,case))
    if max(moved.values())>bound or atom_state(pose,case)!=observations:
        raise ValueError('gradient mapping moved endpoint or protected observations')
    if not math.isfinite(center) or abs(center-before_score)>score_bound:
        raise ValueError('endpoint derivative objective differs from raw score')
    scorer.finalize_after_minimizing(work)
    return dict(components=rows,maximum_per_degree=max(abs(r['analytic_per_degree']) for r in rows),
                raw_score=before_score,minimizer_center_score=center,mapped_motion=moved,original_pose_unchanged=True)


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    if record['schema']!='dvbfixer.chi-endpoint-guard-plan.v1':
        raise ValueError('separate bounded endpoint-guard plan required')
    import pyrosetta.rosetta.protocols.minimization_packing as packing

    original=packing.MinMover
    inventory=[(case,cell) for case in record['cases'] for cell in case['cells']]
    created=[]
    traces=[]
    settings=record['guard_protocol']

    class GuardedMover:
        def __init__(self):
            self.inner=original()
            if self.inner.min_options().gmax_cutoff_for_convergence()!=1.:
                raise ValueError('unexpected default gradient cutoff')
            self.inner.min_options().gmax_cutoff_for_convergence(record['minimizer']['gradient_cutoff'])
            self.case,self.cell=inventory[len(created)]
            created.append(self)

        def __getattr__(self,name):
            return getattr(self.inner,name)

        def apply(self,pose):
            actual=dict(type=self.inner.min_type(),tolerance=self.inner.tolerance(),
                gradient_cutoff=self.inner.min_options().gmax_cutoff_for_convergence(),
                max_iter=self.inner.min_options().max_iter(),cartesian=self.inner.cartesian(),neighbor_list=self.inner.nb_list())
            check_settings(actual,record['minimizer'])
            if self.inner.deriv_check() or self.inner.abs_score_convergence_threshold()!=0.:
                raise ValueError('unexpected native diagnostic/repeat options')
            baseline=Path(settings['parent_run'])/self.case['id']/self.cell['declaration']['id']
            folder=output/self.case['id']/self.cell['declaration']['id']
            initial=atom_state(pose,self.case['case'])
            phases=[]
            for cycle in range(1,settings['total_cycle_limit']+1):
                self.inner.apply(pose)
                atoms=atom_state(pose,self.case['case'])
                score=float(self.inner.score_function()(pose))
                if not math.isfinite(score):
                    raise ValueError('nonfinite guarded raw score')
                if cycle==1:
                    old=json.loads((baseline/'after-atoms.json').read_bytes())['atoms']
                    old_score=json.loads((baseline/'result.json').read_bytes())['final_raw_score']
                    if atoms!=old or score!=old_score:
                        raise ValueError('first inner cycle not an exact mode-C replay')
                moves=movement(initial,atoms)
                if moves['protected_max_motion_angstrom']>record['protected_numerical_bound_angstrom']:
                    raise ValueError('guarded cycle moved protected atoms; no restore')
                gradient=endpoint_gradient(pose,self.inner.movemap(pose),self.inner.score_function(),self.case['case'],
                    settings['mapped_endpoint_numerical_bound_angstrom'],settings['center_score_absolute_bound'])
                if abs(gradient['raw_score']-score)>settings['center_score_absolute_bound']:
                    raise ValueError('cloned gradient score differs from guarded raw pose')
                decision=guard_decision(gradient['maximum_per_degree'],cycle,settings['total_cycle_limit'],settings['cutoff_per_degree'])
                phase=dict(cycle=cycle,atoms=atoms,raw_score=score,gradient=gradient,decision=decision,movement=moves)
                path=folder/f'cycle-{cycle}.json'
                write(path,phase)
                phases.append(dict(cycle=cycle,phase_sha256=sha(path.read_bytes()),decision=decision))
                if decision!='continue':
                    break
            traces.append(dict(id=self.case['id'],index=self.cell['declaration']['index'],phases=phases,
                first_cycle_exact=True,final_decision=phases[-1]['decision'],actual_settings=actual))

    with patch.object(packing,'MinMover',GuardedMover):
        parent.run(plan_path,output)
    if packing.MinMover is not original:
        raise ValueError('endpoint guard adapter not restored')
    write(output/'guard-evidence.json',dict(schema='dvbfixer.chi-endpoint-guard-run.v1',traces=traces,
        plan_sha256=sha(plan_path.read_bytes()),hook_restored=True,native_coordinates_read=False,weights_modified=False,
        enrolled_starts=len(inventory),created_minimizers=len(created)))
    print(json.dumps(dict(starts=len(inventory),completed=len(traces),cycles=sum(len(t['phases']) for t in traces),
        endpoint_gradient_met=sum(t['final_decision']=='endpoint_gradient_met' for t in traces),
        budget_exhausted=sum(t['final_decision']=='budget_exhausted' for t in traces)),indent=2))


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
