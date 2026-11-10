"""Private installed-PyRosetta prior queries, never a table export or production API.

Rank saved geometry-feasible NCAC proposals by their library rotamer-family mass.
No new coordinates, poses, packing, H optimization, score fitting or native input.
Every cohort here has already been evaluated: this new ranking is retrospective.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

FLAGS='-constant_seed -jran 7 -mute all -multithreading:total_threads 1 -packing:pack_missing_sidechains false -packing:no_optH true'
REVISION='1a56185c2592611dec4c9c75ddc9468cd2227c1f'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def write(path,record):
    with path.open('x') as stream:
        json.dump(record,stream,sort_keys=True,indent=2,allow_nan=False)
        stream.write('\n')


def select_prior(candidates,well_probabilities):
    """Prior-only family ranking; geometry first, stable index ties, no RMSD."""
    if set(well_probabilities)!={1,2,3} or any(not math.isfinite(p) or not 0<=p<=1 for p in well_probabilities.values()):
        raise ValueError('three finite family probabilities required')
    if not math.isclose(sum(well_probabilities.values()),1.,abs_tol=5e-6):
        raise ValueError('family probabilities not normalized')
    if [c['index'] for c in candidates]!=[0,1,2] or {c['well'] for c in candidates}!={1,2,3}:
        raise ValueError('exactly one saved candidate per family required')
    feasible=[c for c in candidates if c['common_accepted']]
    return max(feasible,key=lambda c:(well_probabilities[c['well']],-c['index']))['index'] if feasible else None


def plan(output):
    from scripts.compare_branched_completion import measure
    from scripts.research_rosetta_joint import sha as scientific_sha

    from dvbfixer.diagnose.geometry import signed_dihedral_degrees
    from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
    from dvbfixer.protein_completion.backbone import backbone_angles
    from dvbfixer.protein_completion.branched_source import plan_branched_source
    from dvbfixer.protein_completion.contracts import Addition

    sources=[('old_eleven',Path('tmp/pdbfixer-research/ncac-expansion-plan-01.json'),
              Path('tmp/pdbfixer-research/ncac-rotamers-score-plan-01/geometry.json')),
             ('new_inputs_eleven',Path('tmp/pdbfixer-research/ncac-new-system-masks-01/plan.json'),
              Path('tmp/pdbfixer-research/ncac-new-system-run-01/evidence.json'))]
    inputs,hashes,cases={},{},[]
    for group,source_plan_path,geometry_path in sources:
        inputs[str(source_plan_path)]=sha(source_plan_path.read_bytes())
        inputs[str(geometry_path)]=sha(geometry_path.read_bytes())
        source_plan=json.loads(source_plan_path.read_bytes())
        hashes.update(source_plan['source_hashes'])
        geometry=json.loads(geometry_path.read_bytes())
        for original in source_plan['cases']:
            identifier=original['id']
            if group=='old_eleven':
                source_path=Path(original['input_path'])
                candidates=geometry['cases'][identifier]['candidates']
                state='evaluated'
            else:
                folder=geometry_path.parent/identifier
                source_path=folder/'input.pdb'
                row=geometry['results'][identifier]
                state=row.get('ncac_status',row['status'])
                if state!='evaluated':
                    cases.append(dict(id=identifier,group=group,status='parent_not_executed',parent_status=state,
                                      reason=row.get('reason',row.get('error')),source_path=str(source_path)))
                    inputs[str(source_path)]=sha(source_path.read_bytes())
                    continue
                candidate_path=folder/'ncac-proposals.json'
                inputs[str(candidate_path)]=sha(candidate_path.read_bytes())
                candidates=json.loads(candidate_path.read_bytes())['candidates']
            source=source_path.read_bytes()
            inputs[str(source_path)]=sha(source)
            target=ComponentInstanceRef(**original['target'])
            request=plan_branched_source(source,target).request
            try:
                phi,psi=backbone_angles(request)
            except ValueError as error:
                cases.append(dict(id=identifier,group=group,status='backbone_refused',reason=str(error),source_path=str(source_path)))
                continue
            measured=[]
            for i,c in enumerate(candidates):
                points=c['points']
                added=tuple(Addition(ExactAtomRef(target,n),n[0],tuple(points[n])) for n in request.missing_atoms)
                common,_=measure(source,request,added,anchored=False)
                if common.accepted!=c['common_audit']['accepted']:
                    raise ValueError('saved candidate feasibility changed')
                branch='CG1' if request.residue_name=='VAL' else 'OG1'
                chi=signed_dihedral_degrees(*(points[n] for n in ('N','CA','CB',branch)))
                if chi is None:
                    raise ValueError('saved chi is degenerate')
                measured.append(dict(index=i,chi_degrees=chi,common_accepted=common.accepted))
            if len(measured)!=3:
                raise ValueError('all three original relative alternatives required')
            cases.append(dict(id=identifier,group=group,status='query_ready',residue_name=request.residue_name,
                              phi_degrees=phi,psi_degrees=psi,candidates=measured,source_path=str(source_path),
                              target=original['target']))
    for p,d in hashes.items():
        if scientific_sha(Path(p).read_bytes())!=d:
            raise ValueError('prior frozen implementation changed')
    hashes[str(Path(__file__).relative_to(Path.cwd()) if Path(__file__).is_absolute() else Path(__file__))]=sha(Path(__file__).read_bytes())
    output.parent.mkdir(parents=True,exist_ok=True)
    write(output,dict(schema='dvbfixer.pyrosetta-prior-plan.v1',cases=cases,inputs=inputs,source_hashes=hashes,
        init_flags=FLAGS,expected_revision=REVISION,
        library_policy='record active RotamerLibrary binary name/SHA256 and canonical chi definition; do not claim original upstream equivalence',
        query='get_all_rotamer_samples at input-only phi/psi; bin_rotameric_chi for each unchanged candidate; all three wells, no cutoff',
        selection='largest library rotamer-family mass among common-feasible candidates; stable original index tie; no Gaussian offset or ref2015 sum',
        controls='query phi+360 and psi-360, require same well/probability/mean/sigma within 1e-6; separately repeat complete API output',
        scope='retrospective ranking on both already evaluated cohorts; parent refusals/blockers retained, not new qualification',
        native_coordinates_read=False,table_redistribution=False,production_enabled=False))


def run(plan_path,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen query input/implementation changed')
    import pyrosetta

    if record['expected_revision'] not in pyrosetta.version():
        raise ValueError('unexpected PyRosetta revision')
    pyrosetta.init(record['init_flags'])
    from pyrosetta.rosetta.core.chemical import ChemicalManager, aa_from_name
    from pyrosetta.rosetta.core.pack.dunbrack import RotamerLibrary
    from pyrosetta.rosetta.utility import fixedsizearray1_double_5UL_t

    manager=RotamerLibrary.get_instance()
    binary=Path(manager.get_binary_name())
    if not manager.decide_read_from_binary() or not binary.is_file():
        raise ValueError('active installed binary must be present and identified')
    database=Path(pyrosetta.__file__).parent/'database'
    binary_raw=binary.read_bytes()
    environment=dict(version=pyrosetta.version(),binary_relative_path=str(binary.relative_to(database)),
        binary_sha256=sha(binary_raw),binary_bytes=len(binary_raw),read_binary=True,init_flags=record['init_flags'])
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    write(output/'environment.json',environment)
    results={}
    for case in record['cases']:
        name=case['id']
        if case['status']!='query_ready':
            results[name]=dict(status=case['status'],group=case['group'],reason=case.get('reason'),
                               parent_status=case.get('parent_status'))
            continue
        lib=manager.get_library_by_aa(aa_from_name(case['residue_name']))
        residue_type=ChemicalManager.get_instance().residue_type_set('fa_standard').name_map(case['residue_name'])
        atoms=[residue_type.atom_name(i).strip() for i in residue_type.chi_atoms(1)]
        expected=['N','CA','CB','CG1' if case['residue_name']=='VAL' else 'OG1']
        if atoms!=expected or lib.nchi()!=1:
            raise ValueError('library chi convention does not match measured candidate chi')

        def query(phi,psi,lib=lib):
            bb=fixedsizearray1_double_5UL_t()
            bb[1],bb[2]=phi,psi
            values=[]
            for sample in lib.get_all_rotamer_samples(bb):
                if sample.nchi()!=1 or sample.nrchi_sample():
                    raise ValueError('single-chi rotameric sample required')
                values.append(dict(well=int(sample.rot_well()[1]),probability=float(sample.probability()),
                    mean_degrees=float(sample.chi_mean()[1]),sigma_degrees=float(sample.chi_sd()[1])))
            if len(values)!=3 or {v['well'] for v in values}!={1,2,3}:
                raise ValueError('complete three-well API result required')
            if any(not all(math.isfinite(v[k]) for k in ('probability','mean_degrees','sigma_degrees')) or v['sigma_degrees']<=0 for v in values):
                raise ValueError('invalid sample parameters')
            return sorted(values,key=lambda v:v['well'])

        samples=query(case['phi_degrees'],case['psi_degrees'])
        periodic=[query(case['phi_degrees']+360.,case['psi_degrees']),query(case['phi_degrees'],case['psi_degrees']-360.)]
        for check in periodic:
            for a,b in zip(samples,check,strict=True):
                if a['well']!=b['well'] or any(not math.isclose(a[k],b[k],rel_tol=0.,abs_tol=1e-6) for k in ('probability','mean_degrees','sigma_degrees')):
                    raise ValueError('periodic API control differs')
        candidates=[dict(c,well=int(lib.bin_rotameric_chi(c['chi_degrees'],1))) for c in case['candidates']]
        selected=select_prior(candidates,{s['well']:s['probability'] for s in samples})
        results[name]=dict(status='queried',group=case['group'],library_type=type(lib).__name__,chi_atoms=atoms,
            phi_degrees=case['phi_degrees'],psi_degrees=case['psi_degrees'],samples=samples,periodic_controls=periodic,
            candidates=candidates,selected_index=selected)
    if sha(binary.read_bytes())!=environment['binary_sha256']:
        raise ValueError('installed library changed during queries')
    write(output/'selection-before-native.json',dict(results=results,native_coordinates_read=False,
        plan_sha256=sha(raw),environment_sha256=sha((output/'environment.json').read_bytes())))
    print(json.dumps(dict(planned=len(results),queried=sum(r['status']=='queried' for r in results.values()),
                         selected=sum(r.get('selected_index') is not None for r in results.values()),
                         active_library=environment['binary_relative_path']),indent=2))


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
