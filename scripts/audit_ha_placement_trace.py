"""Independent saved-only provenance audit of original H-placement replays."""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_ha_start_control import ha_geometry, inventory
from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_placement_trace import verify


def descriptor(angles):
    if set(angles) != {'N', 'C', 'CB'} or not all(math.isfinite(v) for v in angles.values()):
        raise ValueError('incomplete HA triad')
    if max(angles.values()) < 90:
        return 'all_acute'
    if min(angles.values()) > 90:
        return 'all_obtuse'
    return 'mixed_or_boundary'


def existing_fixed(original, before, after):
    source, first, second = [inventory(rows) for rows in (original, before, after)]
    if not source.keys() <= first.keys() or first.keys() != second.keys():
        raise ValueError('original H-placement atom inventory differs')
    for key, row in source.items():
        for mapped in (first, second):
            if mapped[key]['element'] != row['element'] or mapped[key]['xyz'] != row['xyz']:
                raise ValueError('H-only placement moved an existing observation')
    return dict(existing_atoms=len(source), added_atoms=len(first)-len(source), existing_coordinates_exact=True)


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import dalton

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw:
        raise ValueError('H-placement trace plan differs')
    expected = [(c['id'], c['seed'], c['source_arm']) for c in plan['cases']]
    if [(r['id'], r['seed'], r['source_arm']) for r in run['cases']] != expected:
        raise ValueError('H-placement cohort differs')
    results = []
    for case, row in zip(plan['cases'], run['cases'], strict=True):
        folder = root/case['id']/f"seed-{case['seed']}"/case['source_arm']
        source = Path(case['source'])
        for name, digest in row['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('H-placement artifact changed')
        report = {k: row[k] for k in ('id', 'seed', 'source_arm', 'status', 'returncode')}
        if row['status'] != 'completed':
            results.append(report)
            continue
        trace = json.loads((folder/'placement-trace.json').read_bytes())
        if trace['plan_sha256'] != sha(raw) or trace['case'] != case or not trace['phase1_not_executed']:
            raise ValueError('H-placement trace contract differs')
        if len(trace['h_calls']) != 1 or trace['h_calls'] != [dict(platform='Reference', existing_atoms=trace['h_calls'][0]['existing_atoms'],
                final_atoms=trace['h_calls'][0]['final_atoms'], original_call_forwarded=True, tolerance=1., maximum_iterations=50)]:
            raise ValueError('native H-only minimizer settings differ')
        original, before, after = [load_atoms(folder/name) for name in ('h-input.json.gz', 'h-internal-before.json.gz', 'h-internal-after.json.gz')]
        report['observations'] = existing_fixed(original, before, after)
        system = openmm.XmlSerializer.deserialize((folder/'h-internal-system.xml').read_text())
        original_keys = set(inventory(original))
        if system.getNumParticles() != len(after):
            raise ValueError('H-stage particle count differs')
        for i, atom in enumerate(after):
            mass = system.getParticleMass(i).value_in_unit(dalton)
            if tuple(atom['key']) in original_keys:
                if mass != 0:
                    raise ValueError('existing particle was not immobilized in original H System')
            elif atom['element'] not in {'H', 'D'} or mass <= 0:
                raise ValueError('added particle is not a mobile H')
        if tuple(['A', '16', '', 'GLN', 'HA']) in original_keys:
            raise ValueError('GLN16 HA was not generated')
        endpoint = coordinate_delta(load_atoms(source/'h-1-after.json.gz'), after)
        first_start = coordinate_delta(load_atoms(source/'phase-1-before.json.gz'), load_atoms(folder/'phase-1-before.json.gz'))
        if any(not delta[kind]['exact'] for delta in (endpoint, first_start) for kind in ('heavy', 'hydrogen')):
            raise ValueError('faithful H endpoint/first whole-system start not reproduced exactly')
        serialize, deserialize = openmm.XmlSerializer.serialize, openmm.XmlSerializer.deserialize
        if serialize(deserialize((source/'phase-1-before-system.xml').read_text())) != serialize(deserialize((folder/'phase-1-before-system.xml').read_text())):
            raise ValueError('first whole-system parameters differ from saved original')
        first, second = ha_geometry(before), ha_geometry(after)
        report.update(initial=first, final=second, initial_descriptor=descriptor(first['ha_angles_degrees']),
            final_descriptor=descriptor(second['ha_angles_degrees']),
            H_stage_motion=coordinate_delta(before, after), saved_H_endpoint_difference=endpoint,
            saved_whole_system_start_difference=first_start, original_minimizer_settings=True,
            original_whole_system_parameters_exact=True, whole_system_minimization_executed=False)
        results.append(report)
    completed = [r for r in results if r['status'] == 'completed']
    summary = dict(enrolled=len(results), statuses=dict(Counter(r['status'] for r in results)),
        initial_descriptors=dict(Counter(r['initial_descriptor'] for r in completed)),
        final_descriptors=dict(Counter(r['final_descriptor'] for r in completed)),
        exact_reproductions=len(completed))
    write(output, dict(schema='dvbfixer.ha-placement-trace-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        driver_sha256=plan['source_hashes'][str(Path('scripts/research_ha_placement_trace.py').resolve())],
        openmm_version=plan['openmm_version'], protocol=plan['policy'],
        production_changed=False, native_coordinates_read=False, scientific_qualification=False,
        original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
