"""Saved-only all18 whole-phase CH2 audit with both focal sites and global shape."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.attribute_graph_methyl_ha import frame, projection
from scripts.audit_graph_methyl_control import h_close_pairs, ha_geometry
from scripts.audit_graph_methyl_whole import global_shape
from scripts.audit_legacy_start_controls import coordinate_delta
from scripts.audit_legacy_stereo_trace import determinant
from scripts.audit_methyl_start_control import h_bonds
from scripts.audit_methyl_whole_control import hydrogen_angles
from scripts.preflight_graph_methyl_transport import mapped_start
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify

from dvbfixer.ffutils.dat import DatRecord


def site_geometry(rows,points,target):
    ca,ha,outward = frame(rows,points,[*target,'HA'])
    coordinates = {r['key'][-1]:[10*x for x in points[i]] for i,r in enumerate(rows) if r['key'][:4] == target}
    H = [x-y for x,y in zip(points[ha],points[ca])]
    radius = math.hypot(*H)
    angles = {}
    for name in ('N','C','CB'):
        vector = [x-y for x,y in zip(coordinates[name],coordinates['CA'])]
        cosine = sum(x*y for x,y in zip(H,vector))/(radius*math.hypot(*vector))
        angles[name] = math.degrees(math.acos(max(-1.,min(1.,cosine))))
    return dict(residue=target,CA_volume_angstrom3=determinant(coordinates,'CA','N','C','CB'),
        CA_H_angstrom=10*radius,HA_angles_degrees=angles,outward_cosine=projection(H,outward)/radius)


def audit(plan_path,root,output):
    import openmm

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'],c['seed'],c['source_arm'],arm) for c in plan['cases'] for arm in plan['policy']['arms']]
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw or [(c['id'],c['seed'],c['source_arm'],c['arm']) for c in run['cases']] != expected:
        raise ValueError('all18 native CH2 continuation enrollment differs')
    by_case = {(c['id'],c['seed'],c['source_arm']):c for c in plan['cases']}
    results = []
    for cell in run['cases']:
        case = by_case[(cell['id'],cell['seed'],cell['source_arm'])]
        source,baseline = Path(case['source']),Path(case['baseline'])
        folder = root/cell['id']/f"seed-{cell['seed']}"/cell['source_arm']/cell['arm']
        for name,digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('whole CH2 artifact changed')
        result = {k:cell[k] for k in ('id','seed','source_arm','arm','status','returncode')}
        if cell['status'] != 'completed':
            results.append(result)
            continue
        trace = json.loads((folder/'trace.json').read_bytes())
        if trace['case'] != case or trace['arm'] != cell['arm'] or trace['plan_sha256'] != sha(raw) or trace['emergency_reflection'] or trace['observations_restored']:
            raise ValueError('whole CH2 trace contract differs')
        if [f['stage'] for f in trace['frames']] != ['initial','phase1','phase2']:
            raise ValueError('whole CH2 native phase inventory differs')
        topology = json.loads((source/'whole-topology.json').read_bytes())
        h_topology = json.loads((source/'H-topology.json').read_bytes())
        rows,bonds = topology['atoms'],topology['bonds']
        h_folder = Path(plan['h_root'])/case['id']/f"seed-{case['seed']}"/case['source_arm']/cell['arm']
        h_path = json.loads(gzip.decompress((h_folder/'path.json.gz').read_bytes()))
        initial,repairs = mapped_start(h_topology['atoms'],h_topology['bonds'],h_path['final_nm'],topology['H_to_whole_indices'])
        if trace['post_H_repairs'] != repairs:
            raise ValueError('whole CH2 repair bridge inventory differs')
        initial_rows = [dict(r,xyz=[10*x for x in p]) for r,p in zip(rows,initial,strict=True)]
        added = DatRecord.load(source/'prepared.dat').added_keys()
        observed = [i for i,r in enumerate(rows) if r['element'] not in {'H','D',None} and (r['key'][0],r['key'][1],r['key'][2],r['key'][-1]) not in added]
        result.update(exact_H_bridge=True,frames=[])
        for record in trace['frames']:
            stage = record['stage']
            data_raw,system_raw = (folder/f'{stage}.json.gz').read_bytes(),(folder/f'{stage}-system.xml').read_bytes()
            if sha(data_raw) != record['coordinate_sha256'] or sha(system_raw) != record['system_sha256'] or system_raw != (baseline/f'{stage}-system.xml').read_bytes():
                raise ValueError('whole CH2 snapshot/force/anchor parameters differ')
            data = json.loads(gzip.decompress(data_raw))
            points = data['positions_nm']
            expected_points = json.loads(gzip.decompress((baseline/f'{stage}.json.gz').read_bytes()))['positions_nm']
            exact = points == expected_points
            if exact != record['exact_baseline_nm'] or (cell['arm'] == 'unchanged' and not exact) or (stage == 'initial' and points != initial):
                raise ValueError('whole CH2 exact bridge/unchanged native endpoint differs')
            material = [dict(r,xyz=[10*x for x in p]) for r,p in zip(rows,points,strict=True)]
            system = openmm.XmlSerializer.deserialize(system_raw.decode())
            result['frames'].append(dict(stage=stage,energy_kj_mol=data['energy_kj_mol'],
                raw_global_shape=global_shape(material),rounded_global_shape=global_shape(material,rounded=True),
                ASP49=site_geometry(rows,points,['A','49','','ASP']),ALA151=site_geometry(rows,points,['A','151','','ALA']),
                HA=ha_geometry(rows,bonds,points),H_bonds=h_bonds(system,rows,points),H_angles=hydrogen_angles(system,rows,points),
                close_H_pairs=h_close_pairs(rows,bonds,points),motion_from_full_start=coordinate_delta(initial_rows,material),
                original_heavy_motion=coordinate_delta([initial_rows[i] for i in observed],[material[i] for i in observed]),
                exact_native_baseline_nm=exact))
        results.append(result)
    summary = dict(enrolled=len(results),statuses=dict(Counter(c['status'] for c in results)),arms={})
    for arm in plan['policy']['arms']:
        cells = [c for c in results if c['status'] == 'completed' and c['arm'] == arm]
        final = [c['frames'][-1] for c in cells]
        summary['arms'][arm] = dict(completed=len(cells),exact_native_baseline_cells=sum(all(f['exact_native_baseline_nm'] for f in c['frames']) for c in cells),
            final_negative_CA=sum(len(f['raw_global_shape']['negative_CA']) for f in final),
            final_CA_at_or_below_0_5=sum(len(f['raw_global_shape']['CA_at_or_below_0_5']) for f in final),
            final_THR_beta_at_or_below_0_5=sum(len(f['raw_global_shape']['THR_beta_at_or_below_0_5']) for f in final),
            final_all_acute_HA=sum(len(f['HA']['all_acute']) for f in final))
    write(output,dict(schema='dvbfixer.methylene-whole-control-audit.v1',cases=results,summary=summary,
        plan_sha256=sha(raw),run_sha256=sha(run_raw),policy=plan['policy'],audit_script_sha256=sha(Path(__file__).read_bytes()),
        production_changed=False,emergency_repair_executed=False,scientific_qualification=False,original_thr2_causality_resolved=False))
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    audit(args.plan,args.run,args.output)
