"""Private analytic/finite derivative diagnosis on ALL saved mode-C endpoints.

MinimizerMap uses degree-scaled torsions; no optimizer run, native inputs, new
candidate selection, changed objective, parameters or production capability.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from unittest.mock import patch

from scripts.research_chi_stationarity import save_atoms
from scripts.research_dunbrack_prior import sha, write
from scripts.research_generated_chi_relax import atom_state, movement
from scripts.research_rosetta_fixed_scoring import pose_edges, score_cell

BASE=Path('tmp/pdbfixer-research/chi-stationarity-plan-01.json')


def derivative_comparison(analytic,numeric,absolute=1e-6,relative=1e-3):
    if not all(math.isfinite(v) for v in (analytic,numeric)):
        raise ValueError('nonfinite derivative')
    error=abs(analytic-numeric)
    bound=absolute+relative*max(abs(analytic),abs(numeric))
    return dict(absolute_error=error,reporting_bound=bound,within_reporting_bound=error<=bound)


def plan(output):
    raw=BASE.read_bytes()
    record=json.loads(raw)
    record['inputs'][str(BASE)]=sha(raw)
    for name in ('chi-stationarity-evidence.json','chi-termination-evidence.json'):
        p=Path('docs/research')/name
        record['inputs'][str(p)]=sha(p.read_bytes())
    for root in (Path('tmp/pdbfixer-research/atomtree-source-1a56185'),
                 Path('tmp/pdbfixer-research/atomtree-derivative-source-1a56185')):
        for p in sorted(root.iterdir()):
            record['inputs'][str(p)]=sha(p.read_bytes())
    record['source_hashes']['scripts/research_chi_derivatives.py']=sha(Path(__file__).read_bytes())
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('bound source or endpoint changed')
    record['schema']='dvbfixer.chi-derivatives-plan.v1'
    record['derivative_protocol']=dict(enrollment='ALL42 mode-C endpoints/all15 pools/seven inherited nonexecutions, including geometry failures; no gradient/native filtering',
        api='MinimizerMap setup with only target chi; ScoreFunction.setup_for_minimizing; AtomTreeMultifunc(x) and dfunc(x,g); no Minimizer.run',
        units='torsion_scale_factor=180/pi, internal atom-tree PHI radians mapped to minimizer variables in degrees; report gradients per degree and converted per radian',
        numeric_steps_degrees=record['steps_degrees'],comparison_absolute=1e-6,comparison_relative=1e-3,
        comparison_policy='diagnostic reporting bound only, not a scientific geometry or production qualification gate; preserve all mismatches',
        cutoff_per_degree=1e-4,endpoint_atoms_bound_angstrom=1e-9,
        probe_policy='independent signed vectors from endpoint x, complete vector order/TorsionID inventory and atom traces retained; no accepted output candidates',
        interpretation='unit-aware derivative consistency, not exact termination attribution, convergence certification or native accuracy')
    write(output,record)


def materialize(pyrosetta,sf,record,row,folder):
    from pyrosetta.rosetta.core.id import AtomID
    from pyrosetta.rosetta.numeric import xyzVector_double_t

    original=pyrosetta.pose_from_pdb
    captured=[]

    def capture(path):
        pose=original(path)
        captured.append(pose)
        return pose

    with patch.object(pyrosetta,'pose_from_pdb',capture):
        initial=score_cell(pyrosetta,sf,Path(record['score_plan_root']),folder,row['case'],row['cell']['declaration'])
    initial=json.loads(json.dumps(initial,allow_nan=False))
    if len(captured)!=1 or pyrosetta.pose_from_pdb is not original:
        raise ValueError('initial capture/hook inventory differs')
    if any(row['cell']['saved'].get(k)!=v for k,v in initial.items()) or (folder/'pose-snapshot.json.gz').read_bytes()!=Path(row['cell']['snapshot']).read_bytes():
        raise ValueError('initial complete score/pose not exactly replayed')
    pose=captured[0]
    before=json.loads((Path(row['folder'])/'before-atoms.json').read_bytes())['atoms']
    after=json.loads((Path(row['folder'])/'after-atoms.json').read_bytes())['atoms']
    if atom_state(pose,row['case'])!=before:
        raise ValueError('initial atom trace differs')
    info=pose.pdb_info()
    ids={}
    for i in range(1,pose.total_residue()+1):
        residue=pose.residue(i)
        for a in range(1,residue.natoms()+1):
            key=(info.chain(i),info.number(i),info.icode(i),residue.name3(),residue.atom_name(a).strip())
            if key in ids:
                raise ValueError('duplicate atom identity')
            ids[key]=AtomID(a,i)
    graph=pose_edges(pose,tuple(row['case']['target_key']))
    for a,b in zip(before,after,strict=True):
        if a['key']!=b['key']:
            raise ValueError('endpoint atom ordering differs')
        if a['mobile']:
            pose.set_xyz(ids[tuple(a['key'])],xyzVector_double_t(*b['xyz']))
    score=float(sf(pose))
    if atom_state(pose,row['case'])!=after or pose_edges(pose,tuple(row['case']['target_key']))!=graph:
        raise ValueError('endpoint materialization not exact')
    if not math.isfinite(score) or abs(score-row['result']['final_raw_score'])>record['materialized_score_absolute_bound']:
        raise ValueError('endpoint score recomputation outside declared bound')
    return pose,ids[tuple(row['case']['target_key']+['CA'])].rsd(),after,score


def run(plan_path,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen derivative input changed')
    import pyrosetta
    from pyrosetta.rosetta.core.id import CHI, PHI
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.optimization import AtomTreeMultifunc, MinimizerMap
    from pyrosetta.rosetta.utility import vector1_double

    if record['expected_revision'] not in pyrosetta.version():
        raise ValueError('unexpected PyRosetta revision')
    pyrosetta.init(record['init_flags'])
    sf=pyrosetta.create_score_function(record['score_function'])
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    settings=record['derivative_protocol']
    results=[]

    def vector(values):
        result=vector1_double()
        for value in values:
            result.append(value)
        return result

    for row in record['endpoints']:
        folder=output/row['id']/row['cell']['declaration']['id']
        folder.mkdir(parents=True)
        try:
            pose,target,endpoint,score=materialize(pyrosetta,sf,record,row,folder)
            work=pose.clone()
            mm=MoveMap()
            mm.set_bb(False)
            mm.set_jump(False)
            mm.set_chi(False)
            mm.set_chi(target,True)
            minmap=MinimizerMap()
            minmap.setup(work,mm)
            sf.setup_for_minimizing(work,minmap)
            func=AtomTreeMultifunc(work,minmap,sf,False,False)
            x=vector([0.]*minmap.nangles())
            minmap.copy_dofs_from_pose(work,x)
            values=list(x)
            nodes=list(minmap.dof_nodes())
            if len(nodes)!=work.residue(target).nchi() or len(nodes)!=len(values):
                raise ValueError('target chi/minimizer vector inventory differs')
            mapped=[]
            for index,node in enumerate(nodes,1):
                torsion=node.torsion_id()
                if node.type()!=PHI or torsion.type()!=CHI or torsion.rsd()!=target or not math.isclose(minmap.torsion_scale_factor(node),180/math.pi,rel_tol=0.,abs_tol=1e-12):
                    raise ValueError('unexpected DOF type, target identity or torsion scaling')
                mapped.append(dict(vector_index=index,chi=torsion.torsion(),scale_factor=minmap.torsion_scale_factor(node)))
            if sorted(m['chi'] for m in mapped)!=list(range(1,work.residue(target).nchi()+1)):
                raise ValueError('chi mapping is not complete/unique')
            center=float(func(x))
            gradient=vector([0.]*len(values))
            func.dfunc(x,gradient)
            derivatives=list(gradient)
            if not all(math.isfinite(v) for v in derivatives) or abs(center-score)>record['materialized_score_absolute_bound']:
                raise ValueError('nonfinite gradient or minimizer objective differs')
            base_atoms=atom_state(work,row['case'])
            base_motion=movement(endpoint,base_atoms)
            if max(base_motion.values())>settings['endpoint_atoms_bound_angstrom']:
                raise ValueError('base DOF mapping/derivatives moved endpoint beyond bound')
            base_digest=save_atoms(folder/'mapped-endpoint-atoms.json.gz',base_atoms)
            probes=[]
            for mapped_row in mapped:
                i=mapped_row['vector_index']
                for h in settings['numeric_steps_degrees']:
                    for sign in record['signed_steps']:
                        v=list(values)
                        v[i-1]+=sign*h
                        energy=float(func(vector(v)))
                        atoms=atom_state(work,row['case'])
                        motion=movement(endpoint,atoms)
                        if not math.isfinite(energy) or motion['protected_max_motion_angstrom']>record['protected_bound_angstrom']:
                            raise ValueError('finite-DOF probe nonfinite or protected motion')
                        filename=f'var{i}-h{h}-sign{sign}.json.gz'
                        probes.append(dict(vector_index=i,chi=mapped_row['chi'],step_degrees=h,sign=sign,score=energy,
                            movement=motion,atoms_file=filename,atoms_sha256=save_atoms(folder/filename,atoms)))
            # The scientific saved endpoint is a separate pose and is never reset.
            if atom_state(pose,row['case'])!=endpoint:
                raise ValueError('derivative clones mutated saved endpoint')
            sf.finalize_after_minimizing(work)
            result=dict(status='evaluated',id=row['id'],index=row['result']['index'],initial_exact_replay=True,
                endpoint_atoms_exact=True,endpoint_unchanged=True,endpoint_score=score,minimizer_center_score=center,
                mapped_endpoint_motion=base_motion,mapped_endpoint_atoms_sha256=base_digest,
                variables=values,mapping=mapped,analytic_derivatives_per_degree=derivatives,probes=probes)
        except Exception as error:
            result=dict(status='execution_blocker',id=row['id'],index=row['result']['index'],error_type=type(error).__name__,error=str(error))
        write(folder/'result.json',result)
        results.append(result)
    write(output/'evidence.json',dict(schema='dvbfixer.chi-derivatives-run.v1',results=results,
        plan_sha256=sha(raw),native_coordinates_read=False,weights_modified=False,minimizer_executed=False))
    print(json.dumps(dict(endpoints=len(results),evaluated=sum(r['status']=='evaluated' for r in results),
        blockers=sum(r['status']!='evaluated' for r in results),analytic_dofs=sum(len(r.get('mapping',[])) for r in results)),indent=2))


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
