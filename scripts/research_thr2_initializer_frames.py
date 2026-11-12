"""Finite saved-RNG frame controls: formula algebra only, never native sampling."""
from __future__ import annotations

import argparse
import itertools
import math
import random
import time
from pathlib import Path

from scripts.audit_thr2_h_initialization import initializer, tuple_tree
from scripts.audit_thr2_native_capture import read, verify_files
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_thr2_native_capture import check

ROOT = Path(__file__).resolve().parents[1]
HB = ['A', '2', '', 'THR', 'HB']
ARMS = ('fixed-lab', 'covariant-offset', 'zero-offset')


def determinant(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
            -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))


def rotations():
    result = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            matrix = [[signs[i] if j == perm[i] else 0 for j in range(3)] for i in range(3)]
            if determinant(matrix) == 1:
                result.append(matrix)
    return sorted(result)


def rotate(matrix, point):
    return [sum(row[i]*point[i] for i in range(3)) for row in matrix]


def inverse(matrix, point, shift):
    return [sum(matrix[j][i]*(point[j]-shift[j]) for j in range(3)) for i in range(3)]


def initialize_offset(parent, neighbors, offset):
    if not neighbors or len(offset) != 3:
        raise ValueError('Only the already-observed bonded-parent branch is admitted')
    raw = [sum(parent[k]-p[k] for p in neighbors) for k in range(3)]
    radius = math.sqrt(sum(x*x for x in raw))
    if radius <= 1e-12:
        raise ValueError('Singular heavy-neighbor initializer')
    combined = [0.1*x/radius+y for x, y in zip(raw, offset, strict=True)]
    radius = math.sqrt(sum(x*x for x in combined))
    if radius <= 1e-12:
        raise ValueError('Singular perturbed initializer')
    return [parent[k]+0.1*combined[k]/radius for k in range(3)]


def geometry(poses, center, hydrogen, neighbors, epsilon=1e-10):
    heavy = 1000*determinant([[poses[j][k]-poses[center][k] for k in range(3)] for j in neighbors])
    four = 1000*determinant([[poses[j][k]-poses[hydrogen][k] for k in range(3)] for j in neighbors])
    degenerate = abs(heavy) <= epsilon or abs(four) <= epsilon
    return dict(heavy_volume_angstrom3=heavy, four_ligand_volume_angstrom3=four,
                degenerate=degenerate, sign_conflict=None if degenerate else heavy*four < 0)


def capture_data(source):
    captured = read(source/'capture.json')
    if len(captured['h_calls']) != 1 or captured['calls'][0]['stage'] != 'hydrogen':
        raise ValueError('Original one-call H source is required')
    original = read(source/'h-1-before.json.gz')['atoms']
    atoms = read(source/'call-1-before-atoms.json.gz')['atoms']
    bonds = read(source/'call-1-topology.json')['bonds']
    mapping = {tuple(a['key']):[v/10 for v in a['xyz']] for a in original}
    if len(mapping) != len(original) or any(a['element'] in ('H', 'D') for a in original):
        raise ValueError('Unique H-free source required')
    heavy = {i:mapping[tuple(a['key'])] for i, a in enumerate(atoms) if a['element'] not in ('H', 'D')}
    adjacent = {i:[] for i in range(len(atoms))}
    for a, b in bonds:
        adjacent[a].append(b)
        adjacent[b].append(a)
    rng = random.Random()
    rng.setstate(tuple_tree(captured['h_calls'][0]['initialized_rng']))
    specs = []
    for i, row in enumerate(atoms):
        if row['element'] not in ('H', 'D'):
            continue
        parents = [j for j in adjacent[i] if j in heavy]
        if len(parents) != 1:
            raise ValueError('Unique heavy parent required')
        parent = parents[0]
        neighbors = [j for j in adjacent[parent] if j in heavy]
        if not neighbors:
            raise ValueError('Unobserved unbonded-parent RNG branch')
        specs.append(dict(index=i, parent=parent, neighbors=neighbors, draws=[rng.random() for _ in range(3)]))
    if rng.getstate() != tuple_tree(captured['h_calls'][0]['after_rng']):
        raise ValueError('Saved draw accounting differs')
    hb = [s for s in specs if atoms[s['index']]['key'] == HB]
    if len(hb) != 1:
        raise ValueError('Exact target HB is not uniquely admitted')
    centers = []
    for spec in specs:
        parent = spec['parent']
        if atoms[parent]['element'] != 'C' or len(spec['neighbors']) != 3:
            continue
        if sum(atoms[j]['element'] in ('H', 'D') for j in adjacent[parent]) != 1:
            continue
        item = dict(spec)
        if spec['index'] == hb[0]['index']:
            names = {atoms[j]['key'][-1]:j for j in spec['neighbors']}
            item['neighbors'] = [names[name] for name in ('CA', 'OG1', 'CG2')]
        centers.append(item)
    return dict(atoms=atoms, heavy=heavy, specs=specs, centers=centers, hb=hb[0],
                exact_rng=True, native=[[v/10 for v in row['xyz']] for row in atoms])


def freeze(output):
    parent_path = ROOT/'docs/research/thr2-native-capture-plan.json'
    parent = read(parent_path)
    check(parent)
    source = ROOT/'tmp/pdbfixer-research/thr2-native-capture-run-01'
    execution = read(source/'evidence.json')
    if execution['plan_sha256'] != sha(parent_path.read_bytes()):
        # The public plan and original private plan must be byte-identical.
        raise ValueError('Original execution plan binding differs')
    inputs = {str(parent_path):sha(parent_path.read_bytes()), str(source/'evidence.json'):sha((source/'evidence.json').read_bytes())}
    for name in ('docs/research/thr2-h-initialization-evidence.json', 'docs/research/thr2-hswap-evidence.json',
                 'docs/plans/thr2-initializer-frame-controls.md'):
        p = ROOT/name
        inputs[str(p)] = sha(p.read_bytes())
    sources, heavy_reference = [], None
    for row in execution['native']:
        folder = source/row['label']
        verify_files(folder, row)
        inputs.update({str(folder/name):digest for name, digest in row['artifacts'].items()})
        data = capture_data(folder)
        if len(data['specs']) != 315 or len(data['heavy']) != 327:
            raise ValueError('Original 315-H/327-heavy source inventory differs')
        identity = [a['key'] for a in data['atoms']], data['heavy']
        if heavy_reference is None:
            heavy_reference = identity
        if identity != heavy_reference:
            raise ValueError('Sources are not matched heavy/identity')
        sources.append(dict(label=row['label'], source=str(folder), hydrogens=315, heavy=327,
                            graph_selected_carbon_centers=len(data['centers'])))
    if [s['label'] for s in sources] != ['unseeded', 'seed-1', 'seed-2']:
        raise ValueError('All three saved labels required')
    shifts = [[0., 0., 0.], [1.25, -0.75, 0.5]]
    matrices = rotations()
    script_paths = [Path(__file__).resolve(), ROOT/'scripts/audit_thr2_initializer_frames.py',
                    ROOT/'scripts/audit_thr2_h_initialization.py', ROOT/'tests/test_thr2_initializer_frames.py']
    write(output, dict(schema='dvbfixer.thr2-initializer-frames-plan.v1', parent_plan=str(parent_path),
        inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in script_paths}, sources=sources,
        rotations=matrices, translations_nm=shifts, arms=list(ARMS),
        cases=[dict(label=s['label'], rotation=r, translation=t) for s in sources
               for r in range(len(matrices)) for t in range(len(shifts))],
        policy=dict(frames=144, formula_arms=432, transformed_h_formula_evaluations=136080,
            identity_h_formula_checks=945, maximum_seconds=120, algebra_error_bound_nm=1e-12,
            volume_degenerate_bound_angstrom3=1e-10, new_random_sampling=False, native_calls=0,
            limits='Correlated deterministic frame controls on one structure and three saved draw streams; not native trajectory, probability, CIP or production initializer qualification'),
        production_changed=False, main_merge_authorized=False))


def run(plan_path, output):
    plan = read(plan_path)
    verify(plan)
    check(read(Path(plan['parent_plan'])))
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    data = {s['label']:capture_data(Path(s['source'])) for s in plan['sources']}
    reference, initialization_checks = {}, []
    for label, source in data.items():
        poses = [list(source['heavy'].get(i, (0., 0., 0.))) for i in range(len(source['atoms']))]
        zero = [list(p) for p in poses]
        for spec in source['specs']:
            base, point = initializer(source['heavy'][spec['parent']], [source['heavy'][j] for j in spec['neighbors']], spec['draws'])
            poses[spec['index']], zero[spec['index']] = point, base
        reference[label] = dict(native=poses, zero=zero)
        initialization_checks.append(dict(label=label, exact_saved_rng=True,
            maximum_formula_vs_saved_opencl_angstrom=max(10*math.dist(poses[s['index']], source['native'][s['index']]) for s in source['specs'])))
    started, results = time.monotonic(), []
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for index, case in enumerate(plan['cases']):
        if time.monotonic()-started > plan['policy']['maximum_seconds']:
            results.append(dict(index=index, **case, status='budget_not_executed'))
            continue
        source = data[case['label']]
        matrix, shift = plan['rotations'][case['rotation']], plan['translations_nm'][case['translation']]
        transformed = {i:[x+y for x, y in zip(rotate(matrix, p), shift, strict=True)] for i, p in source['heavy'].items()}
        frames = []
        for arm in plan['arms']:
            row = dict(arm=arm, status='completed')
            try:
                poses = [list(transformed.get(i, (0., 0., 0.))) for i in range(len(source['atoms']))]
                for spec in source['specs']:
                    offset = [0.05*d for d in spec['draws']]
                    if arm == 'covariant-offset':
                        offset = rotate(matrix, offset)
                    elif arm == 'zero-offset':
                        offset = [0., 0., 0.]
                    poses[spec['index']] = initialize_offset(transformed[spec['parent']], [transformed[j] for j in spec['neighbors']], offset)
                restored = [inverse(matrix, p, shift) for p in poses]
                wanted = reference[case['label']]['zero' if arm == 'zero-offset' else 'native']
                errors = [math.dist(a, b) for a, b in zip(restored, wanted, strict=True)]
                hb_index = source['hb']['index']
                centers = [(s, geometry(poses, s['parent'], s['index'], s['neighbors'])) for s in source['centers']]
                target = next(g for s, g in centers if s['index'] == hb_index)
                exact_identity = restored == wanted if matrix == identity and case['translation'] == 0 else None
                row.update(target=target, maximum_full_h_backtransform_error_nm=max(errors[s['index']] for s in source['specs']),
                    hb_backtransform_error_nm=errors[hb_index], exact_identity_formula=exact_identity,
                    maximum_heavy_backtransform_error_nm=max(errors[i] for i in source['heavy']),
                    conflicting_carbon_hydrogens=[source['atoms'][s['index']]['key'] for s, g in centers if g['sign_conflict']],
                    degenerate_carbon_hydrogens=[source['atoms'][s['index']]['key'] for s, g in centers if g['degenerate']])
                compressed(output/f'frame-{index:03d}-{arm}.json.gz', dict(positions_nm=poses))
            except Exception as error:
                row.update(status='error', error_type=type(error).__name__, error=str(error))
            frames.append(row)
        results.append(dict(index=index, **case, status='completed', arms=frames))
    verify(plan)
    write(output/'evidence.json', dict(schema='dvbfixer.thr2-initializer-frames-execution.v1',
        plan_sha256=sha(plan_path.read_bytes()), initialization_checks=initialization_checks, cases=results,
        artifacts={p.name:sha(p.read_bytes()) for p in sorted(output.glob('*.json.gz'))},
        native_calls=0, new_random_sampling=False, production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    else:
        run(args.plan.resolve(), args.output.resolve())
