"""Strict saved-capture adapters for the private carbon-H factorial.

No Context, H generation or coordinate construction. Legacy labels must match
the retained full topology exactly; native harmonic bonds independently bind its
graph. New/mobile ownership is derived from pre-H identities and original masses.
"""
from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

from scripts.graph_methyl_geometry import inventory


def read_json(path):
    path = Path(path)
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def labels(rows):
    result = [dict(key=r['key'], element=r['element']) for r in rows]
    inventory(result, [])
    return result


def validate_native(rows, bonds, pre_rows, initial, final, system):
    import openmm
    from openmm.unit import dalton

    adjacent = inventory(rows, bonds)
    count = len(rows)
    if (system.getNumParticles() != count or system.getNumConstraints()
            or system.usesPeriodicBoundaryConditions()
            or any(system.isVirtualSite(i) for i in range(count))):
        raise ValueError('native H particle/constraint/periodicity contract differs')
    for points in (initial, final):
        if len(points) != count or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
            raise ValueError('native H nm coordinate inventory is malformed')
    forces = [f for f in system.getForces() if isinstance(f, openmm.HarmonicBondForce)]
    if len(forces) != 1:
        raise ValueError('unique native harmonic bond graph unavailable')
    force = forces[0]
    native_edges = [tuple(sorted((int(force.getBondParameters(i)[0]), int(force.getBondParameters(i)[1]))))
        for i in range(force.getNumBonds())]
    if len(native_edges) != len(set(native_edges)) or set(native_edges) != {tuple(sorted(b)) for b in bonds}:
        raise ValueError('captured topology differs from actual native System bond graph')
    prior_labels = labels(pre_rows)
    prior = {tuple(r['key']): r for r in pre_rows}
    indices = {tuple(r['key']): i for i, r in enumerate(rows)}
    if any(tuple(r['key']) not in indices for r in prior_labels):
        raise ValueError('pre-H identity missing; aliases are not permitted')
    for key, row in prior.items():
        i = indices[key]
        if row['element'] != rows[i]['element'] or row['xyz'] != [10*x for x in initial[i]]:
            raise ValueError('pre-H observations differ in element or exact coordinate')
    new = {i for i, r in enumerate(rows) if tuple(r['key']) not in prior}
    mobile = {i for i in range(count) if system.getParticleMass(i).value_in_unit(dalton) > 0}
    if any(rows[i]['element'] not in {'H', 'D'} or len(adjacent[i]) != 1 for i in new):
        raise ValueError('new ownership is not singly bonded hydrogen')
    if mobile != new:
        raise ValueError('native mobile masses differ from newly added H ownership')
    if any(initial[i] != final[i] for i in range(count) if i not in mobile):
        raise ValueError('native H endpoint moved a fixed observation')
    return sorted(new), sorted(mobile)


def adapt(case):
    import openmm

    if case['status'] == 'refused_source':
        raise ValueError('inherited source refusal must not load native coordinates')
    path = read_json(case['baseline_path'])
    atoms = read_json(case['atoms'])['atoms']
    rows = labels(atoms)
    bonds = read_json(case['topology'])['bonds']
    pre = read_json(case['pre_H'])['atoms']
    initial, final = path['initial_nm'], path['final_nm']
    source = Path(case['system']).parent
    if case['context_kind'] == 'geometry-reconstruction-stress':
        full = read_json(source/'phase-1-before.json.gz')['atoms']
        after = read_json(source/'h-internal-after.json.gz')['atoms']
        if labels(full) != rows or labels(after) != rows:
            raise ValueError('legacy native/full topology label bijection unavailable; no guessed aliases')
        if ([r['xyz'] for r in atoms] != [[10*x for x in p] for p in initial]
                or [r['xyz'] for r in after] != [[10*x for x in p] for p in final]):
            raise ValueError('legacy original nm replay start/endpoint differs from capture')
    elif case['context_kind'] != 'prepared-protein-view':
        raise ValueError('unknown capture adapter kind')
    system = openmm.XmlSerializer.deserialize(Path(case['system']).read_text())
    new, mobile = validate_native(rows, bonds, pre, initial, final, system)
    if case['context_kind'] == 'prepared-protein-view':
        ownership = read_json(case['ownership'])
        if (ownership['rows'] != rows or ownership['bonds'] != bonds
                or sorted(ownership['new_h']) != new or sorted(ownership['mobile']) != mobile):
            raise ValueError('V2 native topology/ownership capture differs')
    return dict(rows=rows, bonds=bonds, new_h=new, mobile=mobile,
        initial_nm=initial, final_nm=final, system=system)
