"""Private source-geometry factorial, retaining measured chi and original atoms.

Raw vs distal-only tetrahedral rebuild vs full generated-sidechain rebuild.
Identical enforced H rebuild, declared seeds and Reference numerical controls.
No production mutation or native-accuracy-dependent construction/selection.
"""
from __future__ import annotations

import argparse
import inspect
import json
import os
import random
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from scripts.audit_legacy_comparator import atoms
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify
from scripts.research_legacy_start_controls import compressed, platform_properties, snapshot_atoms

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.geometry import extend, l_cb

PARENT=Path('tmp/pdbfixer-research/legacy-comparator-plan-01.json')
RUN=Path('tmp/pdbfixer-research/legacy-comparator-run-01')
ARMS=('raw','beta_only','full_generated')
SEEDS=(1,2)


def regularize(points,residue,arm):
    if arm not in ARMS or residue not in BRANCHES:
        raise ValueError('declared geometry arm and canonical residue required')
    branch=BRANCHES[residue]
    chi=signed_dihedral_degrees(*(points[n] for n in ('N','CA','CB',branch.first)))
    if chi is None:
        raise ValueError('source chi undefined; no alternate angle selected')
    output=dict(points)
    if arm!='raw' and residue=='THR':
        if arm=='full_generated':
            output['CB']=l_cb(points['N'],points['CA'],points['C'])
        output[branch.first]=extend(output['N'],output['CA'],output['CB'],branch.first_length,TETRAHEDRAL_ANGLE,chi)
        output[branch.second]=extend(output['N'],output['CA'],output['CB'],branch.second_length,TETRAHEDRAL_ANGLE,chi+branch.second_chi_offset)
    rebuilt=signed_dihedral_degrees(*(output[n] for n in ('N','CA','CB',branch.first)))
    if rebuilt is None or abs((rebuilt-chi+180)%360-180)>1e-8:
        raise ValueError('geometry regularization did not preserve measured chi')
    return output,chi


def rewrite(source,target,points,allowed):
    lines=[]
    changed=[]
    for line in source.decode().splitlines(keepends=True):
        if line.startswith(('ATOM  ','HETATM')):
            key=(line[21],line[22:26].strip(),line[26].strip(),line[12:16].strip())
            if key[:3]==(target.chain_id,target.sequence_number,target.insertion_code) and key[-1] in allowed:
                xyz=''.join(f'{v:8.3f}' for v in points[key[-1]])
                if len(xyz)!=24:
                    raise ValueError('rebuilt coordinate exceeds PDB width')
                candidate=line[:30]+xyz+line[54:]
                if candidate!=line:
                    changed.append(key[-1])
                line=candidate
        lines.append(line)
    return ''.join(lines).encode(),changed


def plan(output):
    parent=json.loads(PARENT.read_bytes())
    verify(parent)
    inputs={str(PARENT):sha(PARENT.read_bytes())}
    cases=[]
    for case in parent['cases']:
        folder=RUN/case['id']
        prepared,dat=folder/'prepared.pdb',folder/'prepared.dat'
        source=prepared.read_bytes()
        target=ComponentInstanceRef(**case['target'])
        owned=DatRecord.load(dat).added_keys()
        for name in case['mask']:
            if (target.chain_id,target.sequence_number,target.insertion_code,name) not in owned:
                raise ValueError('planned rebuild would change an originally observed target atom')
        full=atoms(source)
        residue=full[(target,'CA')].residue_name
        points={n:a.position for (r,n),a in full.items() if r==target}
        for arm in ARMS:
            regularize(points,residue,arm)
        for p in (prepared,dat):
            inputs[str(p)]=sha(p.read_bytes())
        cases.append(dict(id=case['id'],target=case['target'],mask=case['mask'],residue=residue,
                          prepared=str(prepared),dat=str(dat)))
    hashes=dict(parent['source_hashes'])
    for name in ('scripts/research_start_geometry.py','scripts/research_legacy_start_controls.py',
                 'scripts/audit_legacy_start_controls.py','src/dvbfixer/protein_completion/geometry.py',
                 'src/dvbfixer/protein_completion/branched.py','src/dvbfixer/protein_completion/validation.py'):
        hashes[name]=sha(Path(name).read_bytes())
    write(output,dict(schema='dvbfixer.private-start-geometry-plan.v1',cases=cases,inputs=inputs,source_hashes=hashes,
        python=sys.executable,protocol=dict(arms=ARMS,seeds=SEEDS,
            source='ALL six original legacy prepared/dat pairs; all target rebuilt heavy names must be dat-owned',
            interventions='THR beta_only rebuilds OG1/CG2 at original CB; full_generated rebuilds CB from observed N/CA/C then branches; both retain original measured chi1 before rounding; VAL arms are identity controls',
            geometry='unchanged illustrative l_cb and branch lengths/angles; not a calibrated geometry library or fitted native conformer',
            h_control='--rebuild-h BOTH arms/all36 cells; fixed Python seeds only around addHydrogens, restore state; explicitly Reference in internal H context AND native minimization',
            cli='unchanged production minimize --no-solvent --atom-naming standard --rebuild-h --platform Reference, same dat/default restraints/iterations',
            capture='source rewrite/ownership/chi evidence, H-after, all native phase coordinates/system XML; no source coordinate restoration',
            admission='all planned geometry attempts execute, even if initial shape/clashes fail; no success-dependent pool selection',
            limitations='source-geometry intervention precedes H reconstruction: H differences are mediators, not held-identical coordinates; no exact causal isolation of beta vs all H interactions',
            timeout_seconds_per_cell=1800,no_adaptation='no new seeds/angle sampling/profile grid/native metrics or production integration'),
        native_coordinates_read=False,production_changed=False))


def materialize(case,arm,folder):
    source=Path(case['prepared']).read_bytes()
    target=ComponentInstanceRef(**case['target'])
    full=atoms(source)
    points={n:a.position for (r,n),a in full.items() if r==target}
    rebuilt,chi=regularize(points,case['residue'],arm)
    allowed=() if case['residue']=='VAL' or arm=='raw' else ('OG1','CG2') if arm=='beta_only' else tuple(case['mask'])
    owned=DatRecord.load(case['dat']).added_keys()
    if any((target.chain_id,target.sequence_number,target.insertion_code,n) not in owned for n in allowed):
        raise ValueError('geometry rebuild changed an unowned heavy atom')
    modified,changed=rewrite(source,target,rebuilt,allowed)
    (folder/'prepared.pdb').write_bytes(modified)
    (folder/'prepared.dat').write_bytes(Path(case['dat']).read_bytes())
    write(folder/'input-change.json',dict(arm=arm,residue=case['residue'],allowed_atoms=allowed,changed_atom_names=changed,
        original_chi_degrees=chi,unrounded_points=rebuilt,source_sha256=sha(source),input_sha256=sha(modified),
        val_identity_control=case['residue']=='VAL',coordinates_restored=False))


def worker(plan_path,case_id,seed,arm,folder):
    import openmm
    from openmm.app import Modeller, Simulation

    from dvbfixer.cli import main

    record=json.loads(plan_path.read_bytes())
    verify(record)
    if seed not in SEEDS or arm not in ARMS:
        raise ValueError('undeclared geometry cell')
    native_h,native_m=Modeller.addHydrogens,Simulation.minimizeEnergy
    h_calls=[]
    phases=[]
    reference=openmm.Platform.getPlatformByName('Reference')

    def hydrogen(model,*args,**kwargs):
        bound=inspect.signature(native_h).bind(model,*args,**kwargs)
        bound.apply_defaults()
        if bound.arguments['platform'] is not None:
            raise ValueError('unexpected preselected H platform')
        kwargs=dict(kwargs,platform=reference)
        saved=random.getstate()
        random.seed(seed)
        try:
            result=native_h(model,*args,**kwargs)
        finally:
            random.setstate(saved)
        i=len(h_calls)+1
        compressed(folder/f'h-{i}-after.json.gz',dict(atoms=snapshot_atoms(model.topology,model.positions)))
        h_calls.append(dict(call=i,seed=seed,platform='Reference',python_random_state_restored=random.getstate()==saved))
        return result

    def capture(simulation,phase,stage):
        state=simulation.context.getState(positions=True,velocities=True,parameters=True,energy=True)
        rows=snapshot_atoms(simulation.topology,state.getPositions())
        compressed(folder/f'phase-{phase}-{stage}.json.gz',dict(atoms=rows))
        (folder/f'phase-{phase}-{stage}-system.xml').write_text(openmm.XmlSerializer.serialize(simulation.system))
        if phase==1 and stage=='before':
            write(folder/'topology.json',dict(bonds=[[a.index,b.index] for a,b in simulation.topology.bonds()]))
            (folder/'initial-integrator.xml').write_text(openmm.XmlSerializer.serialize(simulation.integrator))
            (folder/'initial-state.xml').write_text(openmm.XmlSerializer.serialize(state))
        return dict(stage=stage,platform=platform_properties(simulation.context.getPlatform(),simulation.context))

    def minimize(simulation,*args,**kwargs):
        phase=len(phases)+1
        before=capture(simulation,phase,'before')
        result=native_m(simulation,*args,**kwargs)
        after=capture(simulation,phase,'after')
        phases.append(dict(phase=phase,before=before,after=after))
        return result

    sys.argv=['dvbfixer','minimize',str(folder/'prepared.pdb'),'--dat',str(folder/'prepared.dat'),
        '--no-solvent','--atom-naming','standard','--rebuild-h','--platform','Reference','-o',str(folder/'minimized.pdb')]
    with patch.object(Modeller,'addHydrogens',hydrogen),patch.object(Simulation,'minimizeEnergy',minimize):
        main()
    if Modeller.addHydrogens is not native_h or Simulation.minimizeEnergy is not native_m:
        raise ValueError('geometry-control hooks not restored')
    verify(record)
    write(folder/'trace.json',dict(id=case_id,seed=seed,arm=arm,phases=phases,h_calls=h_calls,
        hooks_restored=True,plan_sha256=sha(plan_path.read_bytes()),native_coordinates_read=False))


def run(plan_path,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    verify(record)
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    env=dict(os.environ,PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH',''),
        PYTHONPATH=os.pathsep.join((str(Path('src').resolve()),str(Path('.').resolve()))))
    results=[]
    for case in record['cases']:
        for seed in record['protocol']['seeds']:
            for arm in record['protocol']['arms']:
                folder=output/case['id']/f'seed-{seed}'/arm
                folder.mkdir(parents=True,exist_ok=False)
                row=dict(id=case['id'],seed=seed,arm=arm)
                try:
                    materialize(case,arm,folder)
                    argv=[record['python'],str(Path(__file__).resolve()),'worker','--plan',str(plan_path.resolve()),
                        '--case',case['id'],'--seed',str(seed),'--arm',arm,'--output',str(folder.resolve())]
                    with (folder/'minimize.log').open('xb') as log:
                        process=subprocess.run(argv,env=env,stdout=log,stderr=subprocess.STDOUT,
                            timeout=record['protocol']['timeout_seconds_per_cell'],check=False)
                    row.update(status='completed' if process.returncode==0 else 'execution_blocker',returncode=process.returncode)
                except Exception as error:
                    row.update(status='execution_blocker',error_type=type(error).__name__,error=str(error))
                row['artifacts']={p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
                write(folder/'result.json',row)
                results.append(row)
    verify(record)
    write(output/'evidence.json',dict(schema='dvbfixer.private-start-geometry-run.v1',cases=results,plan_sha256=sha(raw),
        native_coordinates_read=False,production_changed=False,coordinates_restored=False))
    print(json.dumps(dict(enrolled=len(results),completed=sum(r['status']=='completed' for r in results),
        execution_blockers=sum(r['status']!='completed' for r in results)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--output',type=Path,required=True)
    for name in ('run','worker'):
        p=sub.add_parser(name)
        p.add_argument('--plan',type=Path,required=True)
        p.add_argument('--output',type=Path,required=True)
        if name=='worker':
            p.add_argument('--case',required=True)
            p.add_argument('--seed',type=int,required=True)
            p.add_argument('--arm',choices=ARMS,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.output)
    elif args.command=='run':
        run(args.plan,args.output)
    else:
        worker(args.plan,args.case,args.seed,args.arm,args.output)
