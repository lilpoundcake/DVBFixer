"""All-cohort initial GLN16 local heavy/angle-parameter identity comparison.

Retrospective explanatory readback, not a new acceptance gate or intervention.
Only hashes of local parameter/coordinate inventories are published; no table or
full-force-field export. Every initial context and every H-angle triad is kept.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_legacy_start_controls import load_atoms
from scripts.audit_stereo_force_projection import resolve, verify
from scripts.research_dunbrack_prior import sha, write


def digest(record):
    return sha(json.dumps(record, sort_keys=True, separators=(',', ':'), allow_nan=False).encode())


def triad_descriptor(angles):
    if set(angles) != {'N', 'C', 'CB'} or not all(math.isfinite(v) for v in angles.values()):
        raise ValueError('complete finite N/C/CB to HA triad required')
    if max(angles.values()) < 90:
        return 'all_three_acute'
    if min(angles.values()) > 90:
        return 'all_three_obtuse'
    return 'mixed_or_boundary'


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import kilojoules_per_mole, radian

    plan_raw = plan_path.read_bytes()
    plan = json.loads(plan_raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(plan_raw):
        raise ValueError('context readback plan differs')
    prior = json.loads(Path(plan['prior_readback']).read_bytes())
    if [(c['id'], c['seed'], c['arm']) for c in run['cases']] != [(c['id'], c['seed'], c['arm']) for c in prior['cases']]:
        raise ValueError('context readback cohort differs')
    results = []
    for cell in run['cases']:
        if cell['status'] != 'evaluated':
            raise ValueError('context comparison requires the complete retained angle readback')
        raw_bytes = (root/cell['raw_artifact']).read_bytes()
        if sha(raw_bytes) != cell['raw_sha256']:
            raise ValueError('raw angle terms changed')
        raw = json.loads(raw_bytes)
        state = next(s for s in raw['states'] if (s['phase'], s['side']) == (1, 'before'))
        centre = next(c for c in state['centres'] if c['label'] == 'context_ca')
        source = Path(plan['run_root'])/cell['id']/f"seed-{cell['seed']}"/cell['arm']
        rows = load_atoms(source/'phase-1-before.json.gz')
        system = openmm.XmlSerializer.deserialize((source/'phase-1-before-system.xml').read_text())
        focal = set(resolve(rows, centre))
        parameters, local_atoms = [], set()
        for force in system.getForces():
            if not isinstance(force, openmm.HarmonicAngleForce):
                continue
            for i in range(force.getNumAngles()):
                a, b, c, reference, stiffness = force.getAngleParameters(i)
                if not focal.intersection((a, b, c)):
                    continue
                first, last = sorted((rows[a]['key'], rows[c]['key']))
                parameters.append(dict(atoms=[first, rows[b]['key'], last],
                    reference_radians=reference.value_in_unit(radian),
                    stiffness_kj_mol_rad2=stiffness.value_in_unit(kilojoules_per_mole/radian**2)))
                local_atoms.update((a, b, c))
        parameters.sort(key=lambda record: json.dumps(record, sort_keys=True))
        heavy = sorted((dict(key=rows[i]['key'], xyz=rows[i]['xyz'], element=rows[i]['element'])
            for i in local_atoms if rows[i]['element'] not in {'H', 'D'}), key=lambda r: r['key'])
        angles = {}
        for term in centre['terms']:
            a, b, c = term['atoms']
            if b != ['A', '16', '', 'GLN', 'CA']:
                continue
            partner = c[-1] if a[-1] == 'HA' else a[-1] if c[-1] == 'HA' else None
            if partner in {'N', 'C', 'CB'}:
                if partner in angles:
                    raise ValueError('ambiguous HA angle term')
                angles[partner] = term['observed_degrees']
        results.append(dict(id=cell['id'], seed=cell['seed'], arm=cell['arm'],
            local_heavy_atom_count=len(heavy), local_heavy_coordinate_sha256=digest(heavy),
            local_angle_parameter_sha256=digest(parameters), local_angle_count=len(parameters),
            ha_angles_degrees=angles, ha_descriptor=triad_descriptor(angles),
            H_projection_kj_mol_nm=centre['categories']['H_containing']['projection_kj_mol_nm'],
            heavy_projection_kj_mol_nm=centre['categories']['heavy_only']['projection_kj_mol_nm']))
    summary = dict(cells=len(results), unique_local_heavy_coordinate_sets=len({r['local_heavy_coordinate_sha256'] for r in results}),
        unique_local_angle_parameter_sets=len({r['local_angle_parameter_sha256'] for r in results}),
        ha_descriptors=dict(Counter(r['ha_descriptor'] for r in results)))
    write(output, dict(schema='dvbfixer.stereo-angle-context.v1', cases=results, summary=summary,
        plan_sha256=sha(plan_raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        scope='retrospective ALL24 initial GLN16 local contexts; angle-force attribution, not a trajectory intervention',
        production_changed=False, native_coordinates_read=False, scientific_qualification=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--readback', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.readback, args.output)
