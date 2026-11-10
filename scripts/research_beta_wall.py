"""Private THR-beta wall mechanism/stress experiment on frozen complete starts.

No production imports are modified. Reference-only two native phases, one wall
strength and one declared artificial inversion drive; no native accuracy/grid.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.audit_legacy_start_controls import load_atoms
from scripts.audit_legacy_stereo_trace import state_report
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify
from scripts.research_legacy_start_controls import compressed

CAPTURE=Path('tmp/pdbfixer-research/legacy-start-controls-capture-01')
PARENT=Path('tmp/pdbfixer-research/legacy-start-controls-plan-01.json')
ARMS=('unchanged','wall','drive','drive_wall')
SCALE_NM3=.0025  # 2.5 A^3, explicit fixed numerical scale, not a per-case fit.
FLOOR_A3=1.0
WALL_K=10000.
DRIVE_K=1000.
Q_DEF=('bq=(ax*(by*cz-bz*cy)+ay*(bz*cx-bx*cz)+az*(bx*cy-by*cx))/bscale;'
       'ax=x1-x2;ay=y1-y2;az=z1-z2;bx=x3-x2;by=y3-y2;bz=z3-z2;'
       'cx=x4-x2;cy=y4-y2;cz=z4-z2')


def beta_force(indices,kind):
    from openmm import CustomCompoundBondForce

    if kind=='wall':
        expression='0.5*bwall_k*step(bfloor-bq)*(bfloor-bq)^2;'+Q_DEF
    elif kind=='drive':
        expression='0.5*bdrive_k*(bq+1)^2;'+Q_DEF
    else:
        raise ValueError('undeclared beta force kind')
    force=CustomCompoundBondForce(4,expression)
    force.addGlobalParameter('bscale',SCALE_NM3)
    if kind=='wall':
        force.addGlobalParameter('bwall_k',WALL_K)
        force.addGlobalParameter('bfloor',FLOOR_A3/1000/SCALE_NM3)
        force.setForceGroup(30)
    else:
        force.addGlobalParameter('bdrive_k',DRIVE_K)
        force.setForceGroup(31)
    force.setName('PRIVATE_BETA_'+kind.upper())
    force.addBond(indices,[])
    return force


def admission(stereo):
    target=stereo['target']
    if target['residue'][-1]!='THR':
        return 'not_applicable_non_THR'
    if target['target_cip'].get('CB')!='R' or target['beta_volume_angstrom3']<=0:
        return 'refused_initial_non_R_beta'
    if target['target_cip'].get('CA')!='S' or stereo['negative_ca']:
        return 'refused_initial_ca_stereo'
    return 'admitted_positive_R_beta'


def restraint(system):
    import openmm

    candidates=[f for f in system.getForces() if isinstance(f,openmm.CustomExternalForce)
        and f.getEnergyFunction()=='0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        and [f.getPerParticleParameterName(i) for i in range(f.getNumPerParticleParameters())]==['k','x0','y0','z0']]
    if len(candidates)!=1:
        raise ValueError('exact original positional restraint not uniquely identified')
    return candidates[0]


def plan(output):
    parent=json.loads(PARENT.read_bytes())
    verify(parent)
    run_raw=(CAPTURE/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(PARENT.read_bytes()):
        raise ValueError('original complete-start capture differs')
    inputs={str(PARENT):sha(PARENT.read_bytes()),str(CAPTURE/'evidence.json'):sha(run_raw)}
    declarations={c['id']:c for c in parent['cases']}
    cases=[]
    for result in run['cases']:
        folder=CAPTURE/result['id']/f"seed-{result['seed']}"
        if result['status']!='completed':
            raise ValueError('original complete-start pool has execution blocker')
        for name,d in result['artifacts'].items():
            if sha((folder/name).read_bytes())!=d:
                raise ValueError('original captured start changed')
        capture=json.loads((folder/'capture.json').read_bytes())
        phase=capture['phases'][0]
        if phase['phase']!=1 or phase['settings']!=dict(max_iterations=1000,tolerance_kj_mol_nm=10.):
            raise ValueError('unexpected captured native phase settings')
        prefix=phase['before']['prefix']
        for suffix in ('system.xml','integrator.xml','state.xml','atoms.json.gz','topology.json'):
            path=folder/f'{prefix}-{suffix}'
            inputs[str(path)]=sha(path.read_bytes())
        case=dict(id=result['id'],seed=result['seed'],source=str(folder),prefix=prefix,
            target=declarations[result['id']]['target'])
        rows=load_atoms(folder/f'{prefix}-atoms.json.gz')
        case['initial_stereo']=state_report(rows,case['target'])
        case['admission']=admission(case['initial_stereo'])
        cases.append(case)
    hashes=dict(parent['source_hashes'])
    for path in ('scripts/research_beta_wall.py','scripts/audit_legacy_stereo_trace.py','scripts/audit_legacy_start_controls.py'):
        hashes[path]=sha(Path(path).read_bytes())
    write(output,dict(schema='dvbfixer.private-beta-wall-plan.v1',cases=cases,inputs=inputs,source_hashes=hashes,
        protocol=dict(arms=ARMS,platform='Reference',phases=2,tolerance_kj_mol_nm=10.,max_iterations=1000,
            scale_nm3=SCALE_NM3,wall_k_kj_mol=WALL_K,wall_floor_a3=FLOOR_A3,drive_k_kj_mol=DRIVE_K,
            drive_target_normalized_volume=-1.,
            source='ALL12 saved complete starts; bad initial CB refused, VAL not applicable; preserve both denominators',
            admission='positive named beta volume and CB R; nearplanar positive source allowed as explicit repair/stress research, not production input admission',
            intervention='target THR only; additional finite one-sided wall, artificial drive is labelled mechanism stress NOT original bug replay',
            phase2='reduce ONLY original positional restraint k by10, leave anchors unchanged, same as production; no emergency reflection or output restoration',
            safety='finite wall is not a mathematical guarantee; independently reject bad/nearplanar CA or CB and other raw/rounded geometry failures',
            limits='no native coordinates/errors, profile choice, new seeds/strength grid or production integration; real observed motion reported'),
        native_coordinates_read=False,production_changed=False))


def run(plan_path,output):
    import openmm
    from openmm.unit import angstrom, kilojoules_per_mole, nanometer

    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    verify(plan)
    if plan['protocol']['wall_k_kj_mol']!=WALL_K or plan['protocol']['drive_k_kj_mol']!=DRIVE_K:
        raise ValueError('declared beta-force strength differs')
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results=[]
    for case in plan['cases']:
        source=Path(case['source'])
        prefix=case['prefix']
        template=load_atoms(source/f'{prefix}-atoms.json.gz')
        for arm in plan['protocol']['arms']:
            folder=output/case['id']/f"seed-{case['seed']}"/arm
            folder.mkdir(parents=True,exist_ok=False)
            row=dict(id=case['id'],seed=case['seed'],arm=arm,status=case['admission'])
            if case['admission']=='admitted_positive_R_beta':
                context=None
                try:
                    system=openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
                    integrator=openmm.XmlSerializer.deserialize((source/f'{prefix}-integrator.xml').read_text())
                    if any(f.getForceGroup() in {30,31} for f in system.getForces()):
                        raise ValueError('reserved research force group already used')
                    original_count=system.getNumForces()
                    target=case['target']
                    key=[target['chain_id'],target['sequence_number'],target['insertion_code']]
                    indices=[]
                    for name in ('CA','CB','OG1','CG2'):
                        found=[i for i,a in enumerate(template) if a['key'][:3]==key and a['key'][-1]==name]
                        if len(found)!=1:
                            raise ValueError('ambiguous target beta force indices')
                        indices.append(found[0])
                    if arm in ('wall','drive_wall'):
                        system.addForce(beta_force(indices,'wall'))
                    if arm in ('drive','drive_wall'):
                        system.addForce(beta_force(indices,'drive'))
                    native_restraint=restraint(system)
                    (folder/'system-before.xml').write_text(openmm.XmlSerializer.serialize(system))
                    context=openmm.Context(system,integrator,openmm.Platform.getPlatformByName('Reference'))
                    context.setState(openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()))
                    for stage in ('initial','phase1','phase2'):
                        if stage=='phase2':
                            for i in range(native_restraint.getNumParticles()):
                                index,params=native_restraint.getParticleParameters(i)
                                native_restraint.setParticleParameters(i,index,[params[0]/10,*params[1:]])
                            native_restraint.updateParametersInContext(context)
                        if stage!='initial':
                            openmm.LocalEnergyMinimizer.minimize(context,10*kilojoules_per_mole/nanometer,1000)
                        state=context.getState(positions=True,energy=True)
                        positions=state.getPositions(asNumpy=True).value_in_unit(angstrom).tolist()
                        compressed(folder/f'{stage}.json.gz',dict(atoms=[dict(a,xyz=p) for a,p in zip(template,positions,strict=True)],
                            total_energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)))
                    (folder/'system-after.xml').write_text(openmm.XmlSerializer.serialize(system))
                    row.update(status='completed',original_force_count=original_count,research_force_count=system.getNumForces()-original_count)
                except Exception as error:
                    row.update(status='execution_blocker',error_type=type(error).__name__,error=str(error))
                finally:
                    if context is not None:
                        del context
            row['artifacts']={p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            write(folder/'result.json',row)
            results.append(row)
    verify(plan)
    write(output/'evidence.json',dict(schema='dvbfixer.private-beta-wall-run.v1',cases=results,plan_sha256=sha(raw),
        native_coordinates_read=False,production_changed=False,coordinates_restored=False))
    print(json.dumps(dict(cells=len(results),completed=sum(r['status']=='completed' for r in results),
        refused=sum(r['status'].startswith('refused') for r in results),
        not_applicable=sum(r['status'].startswith('not_applicable') for r in results),
        execution_blockers=sum(r['status']=='execution_blocker' for r in results)),indent=2))


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
