"""Independent saved-only checks and complete-cohort phase-switch reporting."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_stereo_force_projection import historical_verify, verify
from scripts.research_dunbrack_prior import sha, write


def close(first, second):
    if not math.isfinite(first) or not math.isfinite(second) or not math.isclose(first, second, rel_tol=1e-9, abs_tol=1e-6):
        raise ValueError('saved force/energy identity does not hold')


def audit(plan_path, root, output):
    plan_raw = plan_path.read_bytes()
    plan = json.loads(plan_raw)
    verify(plan)
    historical_verify(json.loads(Path(plan['parent_plan']).read_bytes()), Path(plan['historical_root']))
    raw = (root/'evidence.json').read_bytes()
    report = json.loads(raw)
    if report['plan_sha256'] != sha(plan_raw) or (root/'executed-plan.json').read_bytes() != plan_raw:
        raise ValueError('force readback executed plan differs')
    expected = [(r['id'], r['seed'], r['arm']) for r in plan['executions']]
    if [(r['id'], r['seed'], r['arm']) for r in report['cases']] != expected:
        raise ValueError('force readback omitted an enrolled cell')
    switches = []
    for i, cell in enumerate(report['cases'], 1):
        if json.loads((root/f'cell-{i:02d}.json').read_bytes()) != cell:
            raise ValueError('cell report and complete ledger differ')
        if cell['status'] != 'evaluated':
            continue
        if [(s['phase'], s['side']) for s in cell['states']] != [(1, 'before'), (1, 'after'), (2, 'before'), (2, 'after')]:
            raise ValueError('partial/out-of-order state readback')
        for state in cell['states']:
            components = state['force_components']
            if [f['index'] for f in components] != list(range(len(components))):
                raise ValueError('force groups are not exhaustive and unique')
            close(sum(f['energy_kj_mol'] for f in components), state['total_energy_kj_mol'])
            for centre in state['centres']:
                values = centre['component_projections_kj_mol_nm']
                if len(values) != len(components):
                    raise ValueError('force projection inventory differs')
                close(sum(values), centre['total_projection_kj_mol_nm'])
                guard = centre['guard']
                if guard is not None:
                    q = centre['volume_angstrom3'] / 1000 / guard['vref_nm3']
                    close(q, guard['q'])
                    deficit = max(guard['qmin']-q, 0)
                    if (deficit > 0) != guard['active']:
                        raise ValueError('guard activation differs from volume ratio')
                    close(0.5*guard['k']*deficit**2, guard['energy_kj_mol'])
                    close(values[guard['group']], guard['analytic_projection_kj_mol_nm'])
        source = Path(plan['run_root'])/cell['id']/f"seed-{cell['seed']}"/cell['arm']
        identical = (source/'phase-1-after.json.gz').read_bytes() == (source/'phase-2-before.json.gz').read_bytes()
        if not identical:
            raise ValueError('phase switch also changed saved coordinates')
        before, after = cell['states'][1:3]
        if len(before['force_components']) != len(after['force_components']):
            raise ValueError('phase switch force inventory differs')
        external = []
        for left, right in zip(before['force_components'], after['force_components'], strict=True):
            if (left['type'], left['name'], left['index']) != (right['type'], right['name'], right['index']):
                raise ValueError('phase switch force identities differ')
            scale = 0.1 if left['type'] == 'CustomExternalForce' else 1
            if scale == 0.1:
                external.append(left['index'])
            close(scale*left['energy_kj_mol'], right['energy_kj_mol'])
        if len(external) != 1:
            raise ValueError('positional restraint not uniquely identified')
        centres = []
        for left, right in zip(before['centres'], after['centres'], strict=True):
            if left['label'] != right['label'] or left['residue'] != right['residue']:
                raise ValueError('phase switch focal identities differ')
            close(left['volume_angstrom3'], right['volume_angstrom3'])
            for j, (first, second) in enumerate(zip(left['component_projections_kj_mol_nm'], right['component_projections_kj_mol_nm'], strict=True)):
                close(first*(0.1 if j in external else 1), second)
            delta = right['total_projection_kj_mol_nm']-left['total_projection_kj_mol_nm']
            expected_delta = -0.9*left['component_projections_kj_mol_nm'][external[0]]
            close(delta, expected_delta)
            centres.append(dict(label=left['label'], residue=left['residue'],
                volume_angstrom3=left['volume_angstrom3'],
                projection_before_kj_mol_nm=left['total_projection_kj_mol_nm'],
                projection_after_kj_mol_nm=right['total_projection_kj_mol_nm'],
                restraint_projection_before_kj_mol_nm=left['component_projections_kj_mol_nm'][external[0]],
                net_change_kj_mol_nm=delta))
        switches.append(dict(id=cell['id'], seed=cell['seed'], arm=cell['arm'], coordinates_identical=identical, centres=centres))
    centres = [f for c in report['cases'] for s in c.get('states', []) for f in s['centres']]
    ca = [c for c in centres if c['names'][0] == 'CA']
    summary = dict(report['summary'], phase_switches_checked=len(switches),
        ca_centre_states=len(ca), ca_guard_active=sum(c['guard'] is not None and c['guard']['active'] for c in ca),
        ca_nearplanar_angle_states=sum(c['angle_degrees'] is not None and abs(c['angle_degrees']) <= 10 for c in ca))
    write(output, dict(schema='dvbfixer.stereo-force-readback-audit.v1', cases=report['cases'],
        phase_switches=switches, summary=summary, plan_sha256=sha(plan_raw), readback_sha256=sha(raw),
        audit_script_sha256=sha(Path(__file__).read_bytes()), static_evaluation_only=True,
        native_coordinates_read=False, production_changed=False, scientific_qualification=False,
        original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--readback', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.readback, args.output)
