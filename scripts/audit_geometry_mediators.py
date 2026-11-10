"""Retrospective saved-system disaggregation, not a new minimization experiment."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify


def compare_systems(first,second,rows):
    import openmm

    a=openmm.XmlSerializer.deserialize(first)
    b=openmm.XmlSerializer.deserialize(second)
    if a.getNumParticles()!=b.getNumParticles() or a.getNumForces()!=b.getNumForces():
        raise ValueError('system inventory differs across geometry arms')
    forces=[]
    for i in range(a.getNumForces()):
        x,y=a.getForce(i),b.getForce(i)
        if type(x) is not type(y):
            raise ValueError('force ordering/types differ across geometry arms')
        item=dict(type=type(x).__name__,exact=openmm.XmlSerializer.serialize(x)==openmm.XmlSerializer.serialize(y))
        if isinstance(x,openmm.CustomExternalForce):
            if x.getEnergyFunction()!=y.getEnergyFunction() or x.getNumParticles()!=y.getNumParticles():
                raise ValueError('positional restraint definition differs')
            anchors=[]
            weights=[]
            for j in range(x.getNumParticles()):
                p,px=x.getParticleParameters(j)
                q,py=y.getParticleParameters(j)
                if p!=q:
                    raise ValueError('positional restraint identities differ')
                if px[0]!=py[0]:
                    weights.append(rows[p]['key'])
                delta=10*math.dist(px[1:],py[1:])
                if delta>0:
                    anchors.append(dict(atom=rows[p]['key'],delta_angstrom=delta))
            item.update(changed_weights=weights,changed_anchors=anchors)
        if isinstance(x,openmm.CustomCompoundBondForce):
            if x.getEnergyFunction()!=y.getEnergyFunction() or x.getNumBonds()!=y.getNumBonds():
                raise ValueError('CA guard definition differs')
            references=[]
            for j in range(x.getNumBonds()):
                px,vx=x.getBondParameters(j)
                py,vy=y.getBondParameters(j)
                if list(px)!=list(py):
                    raise ValueError('CA guard identities differ')
                if list(vx)!=list(vy):
                    references.append(dict(atoms=[rows[k]['key'] for k in px],before=list(vx),after=list(vy)))
            item['changed_reference_volumes']=references
        forces.append(item)
    return forces


def audit(plan_path,root,output):
    plan_raw=plan_path.read_bytes()
    plan=json.loads(plan_raw)
    verify(plan)
    run_raw=(root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(plan_raw):
        raise ValueError('source geometry run differs')
    for case in run['cases']:
        folder=root/case['id']/f"seed-{case['seed']}"/case['arm']
        for name,digest in case['artifacts'].items():
            if sha((folder/name).read_bytes())!=digest:
                raise ValueError('saved geometry system/state changed')
    comparisons=[]
    for case in plan['cases']:
        for seed in plan['protocol']['seeds']:
            base=root/case['id']/f'seed-{seed}'
            first=load_atoms(base/'raw/phase-1-before.json.gz')
            for arm in ('beta_only','full_generated'):
                second=load_atoms(base/arm/'phase-1-before.json.gz')
                delta=coordinate_delta(first,second)
                forces=compare_systems((base/'raw/phase-1-before-system.xml').read_text(),
                    (base/arm/'phase-1-before-system.xml').read_text(),first)
                comparisons.append(dict(id=case['id'],seed=seed,residue=case['residue'],arm=arm,
                    starting_coordinate_delta=delta,forces=forces,entire_force_list_exact=all(f['exact'] for f in forces)))
    write(output,dict(schema='dvbfixer.start-geometry-mediators.v1',contrasts=comparisons,
        plan_sha256=sha(plan_raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        analysis_timing='retrospective saved-system analysis after source-geometry results; no new optimizer calls',
        native_coordinates_read=False,native_accuracy_reported=False,production_changed=False,
        interpretation='exact force-parameter identity vs geometry-derived anchor/reference mediation; H coordinates measured separately'))
    print(json.dumps(dict(contrasts=len(comparisons),entire_force_list_exact=sum(r['entire_force_list_exact'] for r in comparisons)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
