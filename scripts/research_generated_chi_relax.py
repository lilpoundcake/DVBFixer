"""Private target-chi relaxation from ALL feasible saved library-mean starts.

Only generated distal VAL/THR heavy atoms and target sidechain H may move. CB,
observations, backbone, context H and virtual atoms are protected. No packing,
weight changes, native inputs, production artifacts or coordinate restoration.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from unittest.mock import patch

from scripts.research_dunbrack_prior import sha, write

SCORE_PLAN=Path('tmp/pdbfixer-research/dunbrack-means-generation-01/score-plan/plan.json')
SCORE_RUN=Path('tmp/pdbfixer-research/dunbrack-means-score-run-01')


def is_mobile(key,generated,target,is_hydrogen,is_backbone,is_virtual):
    if is_virtual:
        return False
    if is_hydrogen:
        return tuple(key[:4])==tuple(target) and not is_backbone
    return tuple(key) in {tuple(k) for k in generated} and key[-1]!='CB'


def movement(before,after):
    old={tuple(a['key']):a for a in before}
    new={tuple(a['key']):a for a in after}
    if len(old)!=len(before) or len(new)!=len(after) or old.keys()!=new.keys():
        raise ValueError('atom inventory changed during target-chi relaxation')
    maximum=0.
    mobile_maximum=0.
    for key,a in old.items():
        b=new[key]
        if any(a[k]!=b[k] for k in ('mobile','is_hydrogen','is_virtual','is_backbone')):
            raise ValueError('atom mobility/type classification changed')
        if any(not math.isfinite(v) for v in b['xyz']):
            raise ValueError('nonfinite relaxed coordinate')
        distance=math.dist(a['xyz'],b['xyz'])
        if a['mobile']:
            mobile_maximum=max(mobile_maximum,distance)
        else:
            maximum=max(maximum,distance)
    return dict(protected_max_motion_angstrom=maximum,mobile_max_motion_angstrom=mobile_maximum)


def plan(output):
    raw=SCORE_PLAN.read_bytes()
    prior=json.loads(raw)
    original=json.loads((SCORE_RUN/'evidence.json').read_bytes())
    if original['plan_sha256']!=sha(raw) or (SCORE_RUN/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('original mean score plan/run binding differs')
    means_protocol=Path('tmp/pdbfixer-research/dunbrack-means-plan-01.json')
    means=json.loads(means_protocol.read_bytes())
    parent_cases={c['id']:c for c in means['cases']}
    inputs={str(p):sha(p.read_bytes()) for p in (SCORE_PLAN,SCORE_RUN/'evidence.json',SCORE_RUN/'executed-plan.json',SCORE_RUN/'environment.json',means_protocol,
        Path('tmp/pdbfixer-research/dunbrack-prior-run-02/environment.json'))}
    cases=[]
    for case in prior['cases']:
        name=case['id'].removesuffix('-means')
        source=SCORE_PLAN.parent/case['id']/'source.pdb'
        inputs[str(source)]=sha(source.read_bytes())
        cells=[]
        saved=original['results'][case['id']]
        for cell in case['cells']:
            score=next(c for c in saved['cells'] if c['id']==cell['id'])
            if score['score_status']!='scored' or not score['import_audit']['accepted']:
                raise ValueError('all original feasible starts must have valid scores')
            snapshot=SCORE_RUN/case['id']/cell['id']/'pose-snapshot.json.gz'
            for p in (SCORE_PLAN.parent/cell['pdb'],snapshot):
                inputs[str(p)]=sha(p.read_bytes())
            cells.append(dict(declaration=cell,saved=score,snapshot=str(snapshot)))
        cases.append(dict(id=name,group=parent_cases[name]['group'],case=case,source=str(source),cells=cells))
    hashes=dict(prior['source_hashes'])
    for p in ('scripts/research_generated_chi_relax.py','scripts/research_rosetta_fixed_scoring.py'):
        hashes[p]=sha(Path(p).read_bytes())
    for p,d in hashes.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('original scientific implementation changed')
    write(output,dict(schema='dvbfixer.generated-chi-relax-plan.v1',cases=cases,
        inherited_nonexecutions=[dict(id=c['id'],group=c['group']) for c in means['cases'] if c['status']!='generation_ready'],
        inputs=inputs,source_hashes=hashes,score_plan_root=str(SCORE_PLAN.parent),
        init_flags=prior['init_flags'],expected_revision=prior['expected_rosetta_revision'],score_function=prior['score_function'],
        minimizer=dict(type='lbfgs_armijo_nonmonotone',tolerance=.01,max_iter=200,cartesian=False,neighbor_list=False),
        mobility='all bb/jump/chi disabled except target residue chi; only declared generated distal heavy atoms and target sidechain H allowed; CB protected even though generated',
        protected_numerical_bound_angstrom=1e-9,
        selection='after independent raw/rounded geometry and mobility audit: minimum final raw refined-state total score among all accepted refined candidates; stable original index tie; failed refinements/geometry remain retained',
        score_semantics='raw relaxed state with target sidechain H; NOT fresh rounded-heavy import scoring; no claim equivalent to old fixed-heavy score protocol',
        replay='canonical JSON metadata and complete compressed original pose snapshots must repeat exactly before minimizing',
        scope='retrospective bounded intervention on all 42 feasible starts/15 pools, seven inherited nonexecutions; no new independent qualification',
        native_coordinates_read=False,weights_modified=False,production_enabled=False))


def atom_state(pose,case):
    info=pose.pdb_info()
    rows=[]
    for i in range(1,pose.total_residue()+1):
        residue=pose.residue(i)
        for a in range(1,residue.natoms()+1):
            atom_type=residue.atom_type(a)
            key=[info.chain(i),info.number(i),info.icode(i),residue.name3(),residue.atom_name(a).strip()]
            h,v=atom_type.is_hydrogen(),atom_type.is_virtual()
            bb=residue.atom_is_backbone(a)
            xyz=residue.xyz(a)
            rows.append(dict(key=key,xyz=[xyz.x,xyz.y,xyz.z],is_hydrogen=h,is_virtual=v,is_backbone=bb,
                mobile=is_mobile(key,case['generated'],case['target_key'],h,bb,v)))
    return rows


def run(plan_path,output):
    from scripts.research_rosetta_fixed_scoring import pose_edges, score_cell

    raw=plan_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen chi-relax input changed')
    import pyrosetta

    if record['expected_revision'] not in pyrosetta.version():
        raise ValueError('unexpected PyRosetta version')
    pyrosetta.init(record['init_flags'])
    sf=pyrosetta.create_score_function(record['score_function'])
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.pack.dunbrack import RotamerLibrary
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover

    library=json.loads(Path('tmp/pdbfixer-research/dunbrack-prior-run-02/environment.json').read_bytes())
    binary=Path(RotamerLibrary.get_instance().get_binary_name())
    if sha(binary.read_bytes())!=library['binary_sha256']:
        raise ValueError('installed library differs from original mean queries')
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results={}
    for row in record['cases']:
        cells=[]
        for cell in row['cells']:
            folder=output/row['id']/cell['declaration']['id']
            folder.mkdir(parents=True)
            original_import=pyrosetta.pose_from_pdb
            captured=[]

            def capture(path,original_import=original_import,captured=captured):
                pose=original_import(path)
                captured.append(pose)
                return pose

            try:
                with patch.object(pyrosetta,'pose_from_pdb',capture):
                    repeated=score_cell(pyrosetta,sf,Path(record['score_plan_root']),folder,row['case'],cell['declaration'])
                repeated=json.loads(json.dumps(repeated,allow_nan=False))
                if pyrosetta.pose_from_pdb is not original_import or len(captured)!=1:
                    raise ValueError('pose capture not restored/inventory differs')
                if any(cell['saved'].get(k)!=v for k,v in repeated.items()) or (folder/'pose-snapshot.json.gz').read_bytes()!=Path(cell['snapshot']).read_bytes():
                    raise ValueError('original complete score/state not reproduced exactly')
                pose=captured[0]
                before=atom_state(pose,row['case'])
                write(folder/'before-atoms.json',dict(atoms=before))
                info=pose.pdb_info()
                indices=[i for i in range(1,pose.total_residue()+1) if [info.chain(i),info.number(i),info.icode(i),pose.residue(i).name3()]==row['case']['target_key']]
                if len(indices)!=1:
                    raise ValueError('target residue not uniquely resolved')
                target=indices[0]
                graph_before=pose_edges(pose,tuple(row['case']['target_key']))
                mm=MoveMap()
                mm.set_bb(False)
                mm.set_jump(False)
                mm.set_chi(False)
                mm.set_chi(target,True)
                settings=record['minimizer']
                minimizer=MinMover()
                minimizer.movemap(mm)
                minimizer.score_function(sf)
                minimizer.min_type(settings['type'])
                minimizer.tolerance(settings['tolerance'])
                minimizer.max_iter(settings['max_iter'])
                minimizer.cartesian(settings['cartesian'])
                minimizer.nb_list(settings['neighbor_list'])
                minimizer.apply(pose)
                final_score=float(sf(pose))
                after=atom_state(pose,row['case'])
                write(folder/'after-atoms.json',dict(atoms=after))
                moves=movement(before,after)
                if moves['protected_max_motion_angstrom']>record['protected_numerical_bound_angstrom']:
                    raise ValueError('minimizer moved protected atoms; no restore/fallback')
                if pose_edges(pose,tuple(row['case']['target_key']))!=graph_before or not math.isfinite(final_score):
                    raise ValueError('relaxed graph changed or score is nonfinite')
                terms={str(t):dict(raw=float(pose.energies().total_energies()[t]),weight=float(sf.get_weight(t)))
                       for t in sf.get_nonzero_weighted_scoretypes()}
                result=dict(status='refined',index=cell['declaration']['index'],initial_score=cell['saved']['total_score'],
                    final_raw_score=final_score,score_terms=terms,movement=moves,
                    before_sha256=sha((folder/'before-atoms.json').read_bytes()),after_sha256=sha((folder/'after-atoms.json').read_bytes()),
                    original_pose_exact=True,graph_unchanged=True,only_target_chi_enabled=True,
                    convergence_claim=False,coordinate_restore_performed=False)
            except Exception as error:
                result=dict(status='execution_or_mobility_blocker',index=cell['declaration']['index'],
                            error_type=type(error).__name__,error=str(error),coordinate_restore_performed=False)
            write(folder/'result.json',result)
            cells.append(result)
        results[row['id']]=dict(group=row['group'],cells=cells)
    if sha(binary.read_bytes())!=library['binary_sha256']:
        raise ValueError('installed library changed during refinement')
    write(output/'evidence.json',dict(schema='dvbfixer.generated-chi-relax-run.v1',results=results,
        plan_sha256=sha(raw),native_coordinates_read=False,weights_modified=False))
    print(json.dumps(dict(pools=len(results),starts=sum(len(r['cells']) for r in results.values()),
        refined=sum(c['status']=='refined' for r in results.values() for c in r['cells']),
        blockers=sum(c['status']!='refined' for r in results.values() for c in r['cells'])),indent=2))


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
