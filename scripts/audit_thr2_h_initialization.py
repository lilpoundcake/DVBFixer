"""Reconstruct saved Python RNG draws, not new H sampling or minimization.

All 315 initial H positions in each of the three fresh captures are checked
against the already-bound installed Modeller initializer formula. Report actual
OpenCL transport error, never assert double/single precision bit equality.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

from scripts.audit_stereo_force_projection import volume_gradient
from scripts.audit_thr2_native_capture import read, verify_files
from scripts.research_dunbrack_prior import sha, write
from scripts.research_thr2_native_capture import check


def tuple_tree(value):
    return tuple(tuple_tree(v) for v in value) if isinstance(value, list) else value


def initializer(parent, neighbors, draws):
    if not neighbors or len(draws) != 3 or any(not 0 <= x < 1 for x in draws):
        raise ValueError('only the observed bonded-parent initializer branch is admitted')
    raw = [sum(parent[k]-point[k] for point in neighbors) for k in range(3)]
    length = math.sqrt(sum(x*x for x in raw))
    if length <= 1e-12:
        raise ValueError('singular initializer; do not invent a replacement')
    base = [0.1*x/length for x in raw]
    perturbed = [x+0.05*y for x, y in zip(base, draws, strict=True)]
    length = math.sqrt(sum(x*x for x in perturbed))
    if length <= 1e-12:
        raise ValueError('singular perturbed initializer')
    return [parent[k]+base[k] for k in range(3)], [parent[k]+0.1*perturbed[k]/length for k in range(3)]


def target_geometry(rows):
    points = {r['key'][-1]:[v/10 for v in r['xyz']] for r in rows if r['key'][:4] == ['A', '2', '', 'THR']}
    heavy, _, _ = volume_gradient([points[name] for name in ('CB', 'CA', 'OG1', 'CG2')])
    ligand, _, _ = volume_gradient([points[name] for name in ('HB', 'CA', 'OG1', 'CG2')])
    return dict(heavy_center_volume_angstrom3=1000*heavy, four_ligand_volume_angstrom3=1000*ligand,
        heavy_vs_four_ligand_sign_conflict=heavy*ligand < 0)


def audit(plan_path, root, output):
    parent = read(plan_path)
    check(parent)
    execution = read(root/'evidence.json')
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or [r['label'] for r in execution['native']] != parent['protocol']['labels']:
        raise ValueError('complete three-start source required')
    reports, inputs = [], {}
    for result in execution['native']:
        source = root/result['label']
        verify_files(source, result)
        inputs.update({str(source/name):digest for name, digest in result['artifacts'].items()})
        capture = read(source/'capture.json')
        if len(capture['h_calls']) != 1 or capture['calls'][0]['stage'] != 'hydrogen':
            raise ValueError('exactly one original native H call required')
        h = capture['h_calls'][0]
        original = read(source/'h-1-before.json.gz')['atoms']
        native = read(source/'call-1-before-atoms.json.gz')['atoms']
        topology = read(source/'call-1-topology.json')
        if any(r['element'] in ('H', 'D') for r in original):
            raise ValueError('this saved RNG accounting requires the observed strip-all-H input')
        original_map = {tuple(r['key']):[v/10 for v in r['xyz']] for r in original}
        adjacent = {i:[] for i in range(len(native))}
        for a, b in topology['bonds']:
            adjacent[a].append(b)
            adjacent[b].append(a)
        rng = random.Random()
        rng.setstate(tuple_tree(h['initialized_rng']))
        errors, hb_base = [], None
        for index, row in enumerate(native):
            if row['element'] not in ('H', 'D'):
                continue
            parents = [j for j in adjacent[index] if native[j]['element'] not in ('H', 'D')]
            if len(parents) != 1:
                raise ValueError('initial H does not have exactly one heavy parent')
            heavy = parents[0]
            neighbors = [j for j in adjacent[heavy] if native[j]['element'] not in ('H', 'D')]
            base, predicted = initializer(original_map[tuple(native[heavy]['key'])],
                [original_map[tuple(native[j]['key'])] for j in neighbors], [rng.random() for _ in range(3)])
            errors.append(10*math.dist(predicted, [v/10 for v in row['xyz']]))
            if row['key'] == ['A', '2', '', 'THR', 'HB']:
                hb_base = base
        exact_rng = rng.getstate() == tuple_tree(h['after_rng'])
        if not exact_rng or hb_base is None:
            raise ValueError('saved random draw accounting differs; no added draws or replay')
        base_rows = [dict(r, xyz=[10*v for v in hb_base]) if r['key'] == ['A', '2', '', 'THR', 'HB'] else r for r in native]
        stages = {}
        for name in ('call-1-before-atoms', 'call-1-after-atoms', 'h-1-after', 'call-2-before-atoms', 'call-2-after-atoms', 'call-3-after-atoms'):
            stages[name] = target_geometry(read(source/(name+'.json.gz'))['atoms'])
        reports.append(dict(label=result['label'], initialized_hydrogens=len(errors), python_draws=3*len(errors),
            exact_final_python_rng_state=exact_rng, maximum_initializer_vs_native_transport_error_angstrom=max(errors),
            rms_initializer_vs_native_transport_error_angstrom=math.sqrt(sum(x*x for x in errors)/len(errors)),
            deterministic_unperturbed_hb_formula_geometry=target_geometry(base_rows), stages=stages))
    write(output, dict(schema='dvbfixer.thr2-h-initialization.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()), audit_script_sha256=sha(Path(__file__).read_bytes()),
        inputs=inputs, captures=reports, extra_minimization=False, new_random_sampling=False, production_changed=False,
        volume_definition='det(CA-origin, OG1-origin, CG2-origin), origin=CB for heavy-centre or HB for four-ligand measure; these are distinct geometric observables, not interchangeable CIP claims',
        limits='saved RNG/formula reconstruction of the fresh cohort, not the lost original trajectory or a qualified deterministic H initializer'))
    print(json.dumps(reports, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
