"""Saved-only complete carbon-H whole readback, not a production certificate."""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_graph_methyl_control import h_close_pairs, ha_geometry
from scripts.audit_graph_methyl_whole import global_shape
from scripts.audit_legacy_start_controls import coordinate_delta
from scripts.audit_methyl_start_control import h_bonds
from scripts.audit_methyl_whole_control import hydrogen_angles
from scripts.audit_start_geometry import geometry
from scripts.carbon_H_capture_adapter import adapt, labels, read_json
from scripts.preflight_graph_methyl_transport import mapped_start
from scripts.research_beta_wall import restraint
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify

from dvbfixer.ffutils.dat import DatRecord


def cell_folder(root, cell, baseline=False):
    folder = root/cell['id']/f"seed-{cell['seed']}"/cell['conditioning']
    return folder if baseline else folder/cell['arm']


def check_artifacts(folder, cell):
    for name, digest in cell['artifacts'].items():
        if sha((folder/name).read_bytes()) != digest:
            raise ValueError('whole baseline/continuation artifact changed')


def read_frame(folder, record, count):
    raw = (folder/f"{record['stage']}.json.gz").read_bytes()
    xml = (folder/f"{record['stage']}-system.xml").read_bytes()
    if sha(raw) != record['coordinate_sha256'] or sha(xml) != record['system_sha256']:
        raise ValueError('whole frame hash differs')
    data = read_json(folder/f"{record['stage']}.json.gz")
    points = data['positions_nm']
    if (len(points) != count or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points)
            or not math.isfinite(data['energy_kj_mol'])):
        raise ValueError('whole frame coordinate/energy inventory malformed')
    return data, xml


def trace_contract(trace, case, arm, plan_sha):
    if (trace['case_id'] != [case['id'], case['seed'], case['conditioning']]
            or trace['arm'] != arm or trace['plan_sha256'] != plan_sha
            or trace['emergency_reflection'] or trace['observations_restored'] or trace['production_changed']
            or [r['stage'] for r in trace['frames']] != ['initial', 'phase1', 'phase2']):
        raise ValueError('whole native trace contract differs')


def audit_baselines(plan_path, root):
    import openmm
    from openmm.unit import nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    if openmm.__version__ != plan['openmm_version']:
        raise ValueError('whole readback runtime differs')
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['conditioning']) for c in plan['cases']]
    if (run['plan_sha256'] != sha(raw) or run['mode'] != 'native_baseline'
            or (root/'executed-plan.json').read_bytes() != raw
            or [(c['id'], c['seed'], c['conditioning']) for c in run['cases']] != expected):
        raise ValueError('complete native baseline enrollment differs')
    reports = []
    for case, cell in zip(plan['cases'], run['cases'], strict=True):
        folder = cell_folder(root, cell, baseline=True)
        check_artifacts(folder, cell)
        if case['status'] != 'admitted':
            if cell['status'] != case['status']:
                raise ValueError('whole baseline changed inherited refusal')
            reports.append(dict(id=case['id'], seed=case['seed'], conditioning=case['conditioning'], status=cell['status']))
            continue
        if cell['status'] != 'completed':
            raise ValueError('independent whole baseline blocker; no continuation acceptance')
        spec = case['whole']
        trace = read_json(folder/'trace.json')
        trace_contract(trace, case, 'native_baseline', sha(raw))
        original_state = openmm.XmlSerializer.deserialize(Path(spec['state']).read_text())
        original_nm = [list(p) for p in original_state.getPositions().value_in_unit(nanometer)]
        system = openmm.XmlSerializer.deserialize(Path(spec['system']).read_text())
        positional = restraint(system)
        for frame in trace['frames']:
            stage = frame['stage']
            if stage == 'phase2':
                for i in range(positional.getNumParticles()):
                    particle, parameters = positional.getParticleParameters(i)
                    positional.setParticleParameters(i, particle, [parameters[0]/10, *parameters[1:]])
            data, xml = read_frame(folder, frame, len(spec['rows']))
            if xml.decode() != openmm.XmlSerializer.serialize(system):
                raise ValueError('native baseline changed original force/anchor/phase parameters')
            if stage == 'initial' and data['positions_nm'] != original_nm:
                raise ValueError('independent baseline changed original native start')
            if 'historical_angstrom' in spec:
                prior = read_json(spec['historical_angstrom'][stage])['atoms']
                exact = labels(prior) == spec['rows'] and [r['xyz'] for r in prior] == [[10*x for x in p] for p in data['positions_nm']]
                if not exact or not frame['exact_historical_angstrom']:
                    raise ValueError('legacy historical Angstrom baseline readback differs')
        reports.append(dict(id=case['id'], seed=case['seed'], conditioning=case['conditioning'], status='completed',
            original_nm_start_exact=True, phase_systems_exact=True,
            historical_angstrom_exact='historical_angstrom' in spec))
    return dict(cases=reports, run_sha256=sha(run_raw), plan_sha256=sha(raw),
        audit_script_sha256=sha(Path(__file__).read_bytes()),
        summary=dict(enrolled=len(reports), statuses=dict(Counter(c['status'] for c in reports))))


def historical_replay(case, arm):
    if case['context_kind'] != 'prepared-protein-view':
        return None
    root = Path('tmp/pdbfixer-research')
    if arm == 'all_methyl_zero':
        run_root = root/'graph-methyl-whole-run-01'
        folder = run_root/case['id']/f"seed-{case['seed']}"/'graph_tetrahedral_zero'
        public = Path('docs/research/graph-methyl-whole-centres-evidence.json')
    elif case['id'] == '1TIM' and arm in {'all_methylene', 'combined_zero'}:
        run_root = root/'ASP49-methylene-whole-run-01'
        source = 'unchanged' if arm == 'all_methylene' else 'graph_tetrahedral_zero'
        folder = run_root/'1TIM'/f"seed-{case['seed']}"/source/'all_carbon_CH2'
        public = Path('docs/research/ASP49-methylene-whole-evidence.json')
    else:
        return None
    raw = (run_root/'evidence.json').read_bytes()
    manifest_raw = (folder/'result.json').read_bytes()
    manifest = json.loads(manifest_raw)
    if read_json(public)['run_sha256'] != sha(raw) or manifest not in json.loads(raw)['cases'] or manifest['status'] != 'completed':
        raise ValueError('whole historical replay publication binding differs')
    check_artifacts(folder, manifest)
    return dict(folder=str(folder), run_sha256=sha(raw), manifest_sha256=sha(manifest_raw),
        published_audit_sha256=sha(public.read_bytes()))


def audit(plan_path, root, baseline_root, output):
    import openmm

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    baseline = audit_baselines(plan_path, baseline_root)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['conditioning'], arm) for c in plan['cases'] for arm in plan['policy']['arms']]
    if (run['plan_sha256'] != sha(raw) or run['mode'] != 'continuations'
            or run['baseline_run_sha256'] != baseline['run_sha256'] or (root/'executed-plan.json').read_bytes() != raw
            or [(c['id'], c['seed'], c['conditioning'], c['arm']) for c in run['cases']] != expected):
        raise ValueError('complete whole continuation enrollment differs')
    sources = {name:sha(Path(name).read_bytes()) for name in (
        'scripts/audit_carbon_H_whole.py', 'scripts/audit_graph_methyl_whole.py',
        'scripts/audit_start_geometry.py', 'scripts/audit_methyl_start_control.py',
        'scripts/audit_methyl_whole_control.py', 'scripts/audit_graph_methyl_control.py')}
    declarations = {(c['id'], c['seed'], c['conditioning']):c for c in plan['cases']}
    results = []
    for cell in run['cases']:
        case = declarations[(cell['id'], cell['seed'], cell['conditioning'])]
        folder = cell_folder(root, cell)
        check_artifacts(folder, cell)
        result = {k:cell[k] for k in ('id', 'seed', 'conditioning', 'arm', 'status', 'returncode')}
        result.update(context_kind=case['context_kind'], replay_labels=case.get('replay_labels', []))
        if cell['status'] != 'completed':
            if case['status'] != 'admitted' and cell['status'] != case['status']:
                raise ValueError('whole inherited refusal changed')
            results.append(result)
            continue
        if case['status'] != 'admitted':
            raise ValueError('refused whole case executed')
        spec = case['whole']
        rows, bonds = spec['rows'], spec['bonds']
        trace = read_json(folder/'trace.json')
        trace_contract(trace, case, cell['arm'], sha(raw))
        h_path = read_json(Path(plan['h_root'])/case['id']/f"seed-{case['seed']}"/case['conditioning']/cell['arm']/'path.json.gz')
        h_data = adapt(case)
        initial, repairs = mapped_start(h_data['rows'], h_data['bonds'], h_path['final_nm'], spec['mapping'])
        if trace['post_H_repairs'] != repairs:
            raise ValueError('whole H-repair inventory differs')
        initial_rows = [dict(r, xyz=[10*x for x in p]) for r, p in zip(rows, initial, strict=True)]
        added = DatRecord.load(spec['dat']).added_keys()
        observed = [i for i, r in enumerate(rows) if r['element'] not in {'H', 'D', None}
            and (r['key'][0], r['key'][1], r['key'][2], r['key'][-1]) not in added]
        replay = historical_replay(case, cell['arm'])
        base_folder = cell_folder(baseline_root, cell, baseline=True)
        result.update(exact_H_bridge=True, historical_replay=replay, frames=[])
        for frame in trace['frames']:
            stage = frame['stage']
            data, xml = read_frame(folder, frame, len(rows))
            points = data['positions_nm']
            if xml != (base_folder/f'{stage}-system.xml').read_bytes():
                raise ValueError('whole native baseline phase forces/anchors differ')
            exact = points == read_json(base_folder/f'{stage}.json.gz')['positions_nm']
            if exact != frame['exact_baseline_nm'] or (cell['arm'] == 'unchanged' and not exact):
                raise ValueError('whole unchanged native nm endpoint differs')
            if stage == 'initial' and points != initial:
                raise ValueError('whole initial bridge positions differ')
            if replay is not None:
                prior = Path(replay['folder'])
                if (points != read_json(prior/f'{stage}.json.gz')['positions_nm']
                        or xml != (prior/f'{stage}-system.xml').read_bytes()):
                    raise ValueError('whole known historical replay differs')
            material = [dict(r, xyz=[10*x for x in p]) for r, p in zip(rows, points, strict=True)]
            system = openmm.XmlSerializer.deserialize(xml.decode())
            report = dict(stage=stage, energy_kj_mol=data['energy_kj_mol'], exact_native_baseline_nm=exact,
                raw_global_shape=global_shape(material), rounded_global_shape=global_shape(material, rounded=True),
                HA=ha_geometry(rows, bonds, points), H_bonds=h_bonds(system, rows, points),
                H_angles=hydrogen_angles(system, rows, points), close_H_pairs=h_close_pairs(rows, bonds, points),
                motion_from_full_start=coordinate_delta(initial_rows, material),
                original_heavy_motion=coordinate_delta([initial_rows[i] for i in observed], [material[i] for i in observed]))
            if 'target' in spec:
                report.update(raw_masked_heavy_geometry=geometry(material, bonds, spec['target'], spec['residue']),
                    rounded_masked_heavy_geometry=geometry(material, bonds, spec['target'], spec['residue'], rounded=True))
            result['frames'].append(report)
        results.append(result)
    summary = dict(enrolled=len(results), statuses=dict(Counter(c['status'] for c in results)), strata={})
    complete = {}
    for c in results:
        complete.setdefault((c['id'], c['seed'], c['conditioning']), set())
        if c['status'] == 'completed':
            complete[(c['id'], c['seed'], c['conditioning'])].add(c['arm'])
    matched = {k for k, arms in complete.items() if arms == set(plan['policy']['arms'])}
    summary['complete_matched_contexts'] = len(matched)
    for kind in ('geometry-reconstruction-stress', 'prepared-protein-view'):
        summary['strata'][kind] = {}
        for arm in plan['policy']['arms']:
            cells = [c for c in results if c['status'] == 'completed' and c['context_kind'] == kind and c['arm'] == arm]
            final = [c['frames'][-1] for c in cells]
            summary['strata'][kind][arm] = dict(completed=len(cells),
                matched_contexts=sum((c['id'], c['seed'], c['conditioning']) in matched for c in cells),
                exact_native_baseline_cells=sum(all(f['exact_native_baseline_nm'] for f in c['frames']) for c in cells),
                exact_historical_replays=sum(c['historical_replay'] is not None for c in cells),
                final_negative_CA=sum(len(f['raw_global_shape']['negative_CA']) for f in final),
                final_CA_at_or_below_0_5=sum(len(f['raw_global_shape']['CA_at_or_below_0_5']) for f in final),
                final_THR_beta_at_or_below_0_5=sum(len(f['raw_global_shape']['THR_beta_at_or_below_0_5']) for f in final),
                final_all_acute_HA=sum(len(f['HA']['all_acute']) for f in final),
                masked_gate_scope='existing branched-target actual-frame gate' if kind == 'geometry-reconstruction-stress' else 'not applicable',
                final_raw_masked_heavy_passes=sum(f['raw_masked_heavy_geometry']['accepted'] for f in final)
                    if kind == 'geometry-reconstruction-stress' else None,
                final_rounded_masked_heavy_passes=sum(f['rounded_masked_heavy_geometry']['accepted'] for f in final)
                    if kind == 'geometry-reconstruction-stress' else None)
    verify(plan)
    verify(dict(inputs={}, source_hashes=sources))
    write(output, dict(schema='dvbfixer.carbon-H-whole-audit.v1', cases=results, summary=summary,
        baseline=baseline, plan_sha256=sha(raw), run_sha256=sha(run_raw), policy=plan['policy'],
        audit_source_hashes=sources, production_changed=False, emergency_repair_executed=False,
        scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('baseline', 'continuations'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--run', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'continuations':
            p.add_argument('--baseline-root', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'baseline':
        result = audit_baselines(args.plan, args.run)
        write(args.output, result)
        print(json.dumps(result['summary'], indent=2))
    else:
        audit(args.plan, args.run, args.baseline_root, args.output)
