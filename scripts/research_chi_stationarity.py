"""Frozen native-free finite local probes of ALL saved mode-C chi endpoints.

No minimization, new proposal selection, native inputs or optimizer adaptation.
Materialize only saved mobile atoms; never move/restore protected observations.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path
from unittest.mock import patch

from scripts.research_dunbrack_prior import sha, write
from scripts.research_generated_chi_relax import atom_state, movement
from scripts.research_rosetta_fixed_scoring import pose_edges, score_cell

PARENT=Path('tmp/pdbfixer-research/chi-stopping-plans-01/absolute-1e-6-gmax1e-4.json')
RUN=Path('tmp/pdbfixer-research/chi-stopping-absolute-1e-6-gmax1e-4-run-01')


def central_slopes(center,minus,plus,step_degrees):
    if step_degrees<=0 or not all(math.isfinite(v) for v in (center,minus,plus,step_degrees)):
        raise ValueError('finite energies and positive finite step required')
    h=math.radians(step_degrees)
    return dict(slope_per_radian=(plus-minus)/(2*h),
                curvature_per_radian_squared=((plus-center)+(minus-center))/(h*h),
                best_signed_energy_change=min(minus-center,plus-center))


def plan(output):
    raw=PARENT.read_bytes()
    parent=json.loads(raw)
    inputs=dict(parent['inputs'])
    for p in (PARENT,RUN/'executed-plan.json',RUN/'evidence.json',RUN/'stopping-control-settings.json'):
        inputs[str(p)]=sha(p.read_bytes())
    if (RUN/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('parent executed plan differs')
    evidence=json.loads((RUN/'evidence.json').read_bytes())
    rows=[]
    for case in parent['cases']:
        for cell,result in zip(case['cells'],evidence['results'][case['id']]['cells'],strict=True):
            if result['status']!='refined' or cell['declaration']['index']!=result['index']:
                raise ValueError('all parent endpoints required, including geometry failures')
            folder=RUN/case['id']/cell['declaration']['id']
            for p in (folder/'before-atoms.json',folder/'after-atoms.json',folder/'result.json'):
                inputs[str(p)]=sha(p.read_bytes())
            rows.append(dict(id=case['id'],case=case['case'],cell=cell,folder=str(folder),result=result))
    sources=dict(parent['source_hashes'])
    sources['scripts/research_chi_stationarity.py']=sha(Path(__file__).read_bytes())
    for p,d in {**inputs,**sources}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('bound parent state changed')
    write(output,dict(schema='dvbfixer.chi-stationarity-plan.v1',endpoints=rows,inputs=inputs,source_hashes=sources,
        init_flags=parent['init_flags'],expected_revision=parent['expected_revision'],score_function=parent['score_function'],
        score_plan_root=parent['score_plan_root'],inherited_nonexecutions=parent['inherited_nonexecutions'],
        steps_degrees=[.01,.005],signed_steps=[-1,1],all_target_chi=True,
        protected_bound_angstrom=1e-9,materialized_score_absolute_bound=1e-8,descent_reporting_bound=1e-7,
        protocol='ALL42 mode-C endpoints, every target chi including proton chi, each signed step from independent clone of saved endpoint; no geometry/native filtering; no candidate/profile selection',
        materialization='initial pose/scores exactly replayed; set_xyz ONLY saved mobile atoms, require entire endpoint atom trace exact and topology unchanged; energy recomputation bound1e-8 is not exact score replay',
        interpretation='finite directional energy sensitivity, not analytic derivative validation, convergence proof, optimizer choice or native accuracy',
        native_coordinates_read=False,production_enabled=False))


def save_atoms(path,atoms):
    path.write_bytes(gzip.compress(json.dumps(atoms,sort_keys=True,allow_nan=False).encode(),mtime=0))
    return sha(path.read_bytes())


def run(plan_path,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen stationarity state changed')
    import pyrosetta
    from pyrosetta.rosetta.core.id import AtomID
    from pyrosetta.rosetta.numeric import xyzVector_double_t

    if record['expected_revision'] not in pyrosetta.version():
        raise ValueError('unexpected PyRosetta revision')
    pyrosetta.init(record['init_flags'])
    sf=pyrosetta.create_score_function(record['score_function'])
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results=[]
    for row in record['endpoints']:
        folder=output/row['id']/row['cell']['declaration']['id']
        folder.mkdir(parents=True)
        try:
            original_import=pyrosetta.pose_from_pdb
            captured=[]

            def capture(path,captured=captured,original_import=original_import):
                pose=original_import(path)
                captured.append(pose)
                return pose

            with patch.object(pyrosetta,'pose_from_pdb',capture):
                initial=score_cell(pyrosetta,sf,Path(record['score_plan_root']),folder,row['case'],row['cell']['declaration'])
            initial=json.loads(json.dumps(initial,allow_nan=False))
            if len(captured)!=1 or pyrosetta.pose_from_pdb is not original_import:
                raise ValueError('initial pose capture differs or hook not restored')
            if any(row['cell']['saved'].get(k)!=v for k,v in initial.items()) or (folder/'pose-snapshot.json.gz').read_bytes()!=Path(row['cell']['snapshot']).read_bytes():
                raise ValueError('initial complete pose/score not exactly replayed')
            pose=captured[0]
            before=json.loads((Path(row['folder'])/'before-atoms.json').read_bytes())['atoms']
            after=json.loads((Path(row['folder'])/'after-atoms.json').read_bytes())['atoms']
            if atom_state(pose,row['case'])!=before:
                raise ValueError('initial all-atom trace differs')
            ids={}
            info=pose.pdb_info()
            for i in range(1,pose.total_residue()+1):
                residue=pose.residue(i)
                for a in range(1,residue.natoms()+1):
                    key=(info.chain(i),info.number(i),info.icode(i),residue.name3(),residue.atom_name(a).strip())
                    if key in ids:
                        raise ValueError('duplicate atom identity')
                    ids[key]=AtomID(a,i)
            original_graph=pose_edges(pose,tuple(row['case']['target_key']))
            for a,b in zip(before,after,strict=True):
                if a['key']!=b['key']:
                    raise ValueError('saved endpoint atom ordering differs')
                if a['mobile']:
                    pose.set_xyz(ids[tuple(a['key'])],xyzVector_double_t(*b['xyz']))
            endpoint=atom_state(pose,row['case'])
            if endpoint!=after or pose_edges(pose,tuple(row['case']['target_key']))!=original_graph:
                raise ValueError('saved endpoint materialization not exact or graph differs')
            score=float(sf(pose))
            if not math.isfinite(score) or abs(score-row['result']['final_raw_score'])>record['materialized_score_absolute_bound']:
                raise ValueError('materialized endpoint score outside declared recomputation bound')
            if atom_state(pose,row['case'])!=after:
                raise ValueError('endpoint scoring moved atoms')
            target=ids[tuple(row['case']['target_key']+['CA'])].rsd()
            probes=[]
            for chi in range(1,pose.residue(target).nchi()+1):
                for h in record['steps_degrees']:
                    for sign in record['signed_steps']:
                        probe=pose.clone()
                        probe.set_chi(chi,target,pose.chi(chi,target)+sign*h)
                        energy=float(sf(probe))
                        atoms=atom_state(probe,row['case'])
                        moves=movement(endpoint,atoms)
                        if not math.isfinite(energy) or moves['protected_max_motion_angstrom']>record['protected_bound_angstrom'] or pose_edges(probe,tuple(row['case']['target_key']))!=original_graph:
                            raise ValueError('probe nonfinite, protected motion or changed graph')
                        filename=f'chi{chi}-h{h}-sign{sign}.json.gz'
                        probes.append(dict(chi=chi,step_degrees=h,sign=sign,score=energy,movement=moves,
                            atoms_file=filename,atoms_sha256=save_atoms(folder/filename,atoms)))
            if atom_state(pose,row['case'])!=endpoint:
                raise ValueError('probing mutated the saved endpoint')
            result=dict(status='probed',id=row['id'],index=row['result']['index'],endpoint_score=score,
                saved_endpoint_score=row['result']['final_raw_score'],target_chi_count=pose.residue(target).nchi(),
                endpoint_unchanged=True,initial_exact_replay=True,endpoint_atoms_exact=True,probes=probes)
        except Exception as error:
            result=dict(status='execution_blocker',id=row['id'],index=row['result']['index'],error_type=type(error).__name__,error=str(error))
        write(folder/'result.json',result)
        results.append(result)
    write(output/'evidence.json',dict(schema='dvbfixer.chi-stationarity-run.v1',results=results,
        plan_sha256=sha(raw),native_coordinates_read=False,weights_modified=False))
    print(json.dumps(dict(endpoints=len(results),probed=sum(r['status']=='probed' for r in results),
        blockers=sum(r['status']!='probed' for r in results),probes=sum(len(r.get('probes',[])) for r in results)),indent=2))


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
