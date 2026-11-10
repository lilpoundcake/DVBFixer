"""Private retrospective relative-rotamer protocol; no native-based tuning.

Freeze the protocol before generating candidates. Reuse the unchanged fixed-heavy
runner, exposing only common-geometry-feasible candidates and no native controls.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from dataclasses import asdict
from pathlib import Path

from scripts.audit_template_causality import read_pdb
from scripts.compare_branched_completion import measure
from scripts.research_rosetta_joint import sha
from scripts.research_rosetta_repack import write_json

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition

EXPANSION = Path('tmp/pdbfixer-research/ncac-expansion-plan-01.json')
RUN = Path('tmp/pdbfixer-research/ncac-expansion-run-01')
EVIDENCE = Path('docs/research/ncac-whole-mask-evidence.json')
SCORE_PLAN = Path('tmp/pdbfixer-research/independent-score-plan-01/plan.json')
ANGLES = (0., 120., -120.)


def rotate_branches(points, angle):
    """Proper rigid CA–CB rotation, preserving backbone and CB exactly."""
    if not math.isfinite(angle) or angle not in ANGLES:
        raise ValueError('only the three predeclared relative angles are allowed')
    if not {'N', 'CA', 'C', 'O', 'CB'} <= points.keys():
        raise ValueError('complete observed backbone and generated CB required')
    distal = set(points)-{'N', 'CA', 'C', 'O', 'CB'}
    if distal not in ({'CG1', 'CG2'}, {'OG1', 'CG2'}):
        raise ValueError('only complete VAL/THR distal pairs supported')
    result = {n: list(xyz) for n, xyz in points.items()}
    axis = [b-a for a,b in zip(points['CA'],points['CB'],strict=True)]
    length = math.sqrt(sum(x*x for x in axis))
    if length < 1e-8:
        raise ValueError('degenerate rotation axis')
    if angle == 0.:
        return result
    u = [x/length for x in axis]
    c,s = math.cos(math.radians(angle)),math.sin(math.radians(angle))
    for n in distal:
        v = [x-y for x,y in zip(points[n],points['CB'],strict=True)]
        dot = sum(a*b for a,b in zip(u,v,strict=True))
        cross = [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
        result[n] = [points['CB'][i]+c*v[i]+s*cross[i]+(1-c)*dot*u[i] for i in range(3)]
    return result


def freeze(output):
    expansion = json.loads(EXPANSION.read_bytes())
    evidence = json.loads(EVIDENCE.read_bytes())
    geometry = EVIDENCE.with_name(evidence['geometry_ledger_name'])
    if sha(geometry.read_bytes()) != evidence['geometry_ledger_sha256']:
        raise ValueError('parent geometry ledger changed')
    inputs = {str(p):sha(p.read_bytes()) for p in (EXPANSION,RUN/'evidence.json',EVIDENCE,geometry,SCORE_PLAN)}
    hashes = dict(expansion['source_hashes'])
    for p in ('scripts/research_ncac_rotamers.py','scripts/research_rosetta_fixed_scoring.py',
              'scripts/rosetta_fixed_contract.py','scripts/research_rosetta_repack.py'):
        hashes[p] = sha(Path(p).read_bytes())
    record = dict(schema='dvbfixer.ncac-relative-rotamer-protocol.v1', angles_degrees=ANGLES,
                  cases=expansion['cases'], inputs=inputs, source_hashes=hashes, run_path=str(RUN),
                  prior_score_plan=str(SCORE_PLAN), expansion_plan=str(EXPANSION),
                  intervention='rotate both generated distal branches rigidly together around saved CA-CB; CB/backbone/all observations fixed; no optimization',
                  selection='minimum unchanged ref2015 total among ALL common-feasible rounded candidates; missing feasible score blocks selection; stable index tie',
                  native_policy='no native coordinates opened in generation/scoring; no native controls in runner; read native references only after all selections committed',
                  comparator='unchanged NCAC zero-angle and previous own-scored selections; no native tuning',
                  scope='all eleven previously evaluated whole masks, retrospective; no independent-family or public qualification',
                  public_capability=False)
    output.parent.mkdir(parents=True,exist_ok=True)
    write_json(output,record)


def generate(protocol, output):
    raw = protocol.read_bytes()
    record = json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes()) != d:
            raise ValueError('frozen rotamer input/implementation changed')
    if record['angles_degrees'] != list(ANGLES):
        raise ValueError('relative angles differ from declared protocol')
    prior = json.loads(Path(record['prior_score_plan']).read_bytes())
    old_cases = {c['id']:c for c in prior['cases']}
    run_root = Path(record['run_path'])
    run = json.loads((run_root/'evidence.json').read_bytes())
    output.mkdir(parents=True,exist_ok=False)
    (output/'protocol.json').write_bytes(raw)
    cases,ledger,artifacts = [],{},{}
    for case in record['cases']:
        name = case['id']
        source = Path(case['input_path']).read_bytes()
        if sha(source) != case['input_sha256']:
            raise ValueError('masked source changed')
        target = ComponentInstanceRef(**case['target'])
        request = plan_branched_source(source,target).request
        folder = run_root/(name+'-ncac')
        trace_raw = (folder/'evidence.json').read_bytes()
        if sha(trace_raw) != run['results'][name+'-ncac']['trace_sha256']:
            raise ValueError('saved NCAC trace changed')
        trace = json.loads(trace_raw)
        event = [e for e in trace['events'] if e['stage']=='returned-seed'][-1]
        payload = (folder/event['coordinates']).read_bytes()
        if sha(payload) != event['coordinate_sha256']:
            raise ValueError('saved NCAC full-coordinate payload changed')
        atoms = json.loads(gzip.decompress(payload))['atoms']
        key = (target.chain_id,target.sequence_number,target.insertion_code)
        points = {a['key'][3]:a['xyz'] for a in atoms if tuple(a['key'][:3])==key}
        observed = {tuple(a['key']):a['xyz'] for a in read_pdb(source)}
        complete = {tuple(a['key']):a['xyz'] for a in atoms}
        if any(math.dist(complete[k],p)>1e-6 for k,p in observed.items()):
            raise ValueError('saved proposal moved observations')
        old = old_cases[name]
        if old['target'] != case['target'] or set(request.missing_atoms) != set(case['mask']):
            raise ValueError('scorer/source identities differ')
        rows,cells = [],[]
        for i,angle in enumerate(record['angles_degrees']):
            proposal = rotate_branches(points,angle)
            added = tuple(Addition(ExactAtomRef(target,n),n[0],tuple(proposal[n])) for n in request.missing_atoms)
            audit,pdb = measure(source,request,added,anchored=False)
            filename = f'{name}/own-{i}.pdb'
            path = output/filename
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(pdb)
            artifacts[filename] = sha(pdb)
            cell = dict(id=f'own-{i}',arm='own',index=i,pdb=filename,old_source_accepted=audit.accepted,
                        common_accepted=audit.accepted,common_audit=asdict(audit))
            rows.append(dict(cell,relative_angle_degrees=angle,points=proposal))
            if audit.accepted:
                cells.append(cell)
        source_path = output/name/'source.pdb'
        source_path.write_bytes(source)
        artifacts[f'{name}/source.pdb'] = sha(source)
        # Deliberately no reference path or native-control cell is transported.
        cases.append(dict(id=name,status='admitted',cells=cells,target_key=old['target_key'],
                          target_internal=old['target_internal'],expected_external=old['expected_external'],
                          generated=old['generated']))
        ledger[name] = dict(candidates=rows,ncac_previously_executed=case['ncac_previously_executed'])
    write_json(output/'geometry.json',dict(schema='dvbfixer.ncac-relative-rotamer-geometry.v1',cases=ledger,
                                        protocol_sha256=sha(raw),native_coordinates_read=False))
    artifacts['geometry.json'] = sha((output/'geometry.json').read_bytes())
    artifacts['protocol.json'] = sha(raw)
    settings = {k:prior[k] for k in ('expected_rosetta_revision','score_function','seed','radii','init_flags',
                'hydrogen_protocol','extra_heavy_policy','virtual_atom_policy','topology_scope')}
    write_json(output/'plan.json',dict(schema='dvbfixer.rosetta-fixed-scoring-plan.v2',cases=cases,
                artifacts=artifacts,source_hashes=record['source_hashes'],**settings,
                selection_policy=record['selection'],native_policy=record['native_policy'],
                retrospective=True,public_capability=False))
    print(json.dumps(dict(cases=len(cases),proposals=sum(len(c['candidates']) for c in ledger.values()),
                         feasible=sum(len(c['cells']) for c in cases)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('freeze')
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('generate')
    p.add_argument('--protocol',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='freeze':
        freeze(args.output)
    else:
        generate(args.protocol,args.output)
