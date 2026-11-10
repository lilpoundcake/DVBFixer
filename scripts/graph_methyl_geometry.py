"""Private graph-defined methyl geometry; no production call sites or selection.

Only explicitly new, mobile H are owned. Native topology indices and full
case-sensitive identities are retained, never inferred from hydrogen aliases.
"""
from __future__ import annotations

import math

from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS

PHASES = tuple(range(0, 360, 15))


def inventory(rows, bonds):
    identities = [tuple(row['key']) for row in rows]
    if any(len(key) != 5 or any(not isinstance(part, str) for part in key) for key in identities) or len(set(identities)) != len(rows):
        raise ValueError('ambiguous graph atom identity')
    adjacency = [set() for _ in rows]
    for a, b in bonds:
        if (type(a) is not int or type(b) is not int or min(a, b) < 0 or max(a, b) >= len(rows)
                or a == b or b in adjacency[a]):
            raise ValueError('invalid/duplicate native topology bond')
        adjacency[a].add(b)
        adjacency[b].add(a)
    return adjacency


def discover(rows, bonds, newly_added, mobile):
    adjacency = inventory(rows, bonds)
    count = len(rows)
    for selection in (newly_added, mobile):
        if any(type(i) is not int or not 0 <= i < count for i in selection):
            raise ValueError('invalid new/mobile ownership index')
    if any(rows[i]['element'] not in {'H', 'D'} for i in newly_added):
        raise ValueError('newly added H ownership contains a heavy atom')
    groups, exclusions = [], []
    for parent in sorted(range(count), key=lambda i: tuple(rows[i]['key'])):
        row = rows[parent]
        if row['element'] != 'C':
            continue
        neighbors = adjacency[parent]
        hydrogens = sorted((i for i in neighbors if rows[i]['element'] in {'H', 'D'}), key=lambda i: tuple(rows[i]['key']))
        heavy = [i for i in neighbors if rows[i]['element'] not in {'H', 'D'}]
        if len(hydrogens) != 3 or len(heavy) != 1:
            continue
        if not set(hydrogens) <= set(newly_added) & set(mobile):
            exclusions.append(dict(parent=parent, key=row['key'], reason='measured/kept/nonmobile H; no rewrite'))
            continue
        axis = heavy[0]
        residue = tuple(row['key'][:4])
        name = VARIANT_TO_PARENT.get(row['key'][3], row['key'][3])
        if (name not in CANONICAL_HEAVY_ATOMS or row['key'][-1] not in CANONICAL_HEAVY_ATOMS[name]
                or tuple(rows[axis]['key'][:4]) != residue
                or any(tuple(rows[i]['key'][:4]) != residue or adjacency[i] != {parent} for i in hydrogens)):
            raise ValueError('generated methyl outside canonical intra-residue context')
        frames = [i for i in adjacency[axis] if i != parent and rows[i]['element'] not in {'H', 'D'}]
        frames.sort(key=lambda i: (tuple(rows[i]['key'][:4]) != residue, tuple(rows[i]['key'])))
        if not frames:
            raise ValueError('undefined second-shell heavy methyl frame')
        groups.append(dict(parent=parent, axis=axis, frame=frames[0], hydrogens=hydrogens,
            parent_key=row['key'][:], frame_key=rows[frames[0]]['key'][:]))
    owned = [i for group in groups for i in group['hydrogens']]
    if len(set(owned)) != len(owned):
        raise ValueError('overlapping graph methyl ownership')
    return dict(groups=groups, exclusions=exclusions)


def unit(vector):
    norm = math.hypot(*vector)
    if not math.isfinite(norm) or norm <= 1e-12:
        raise ValueError('undefined methyl frame or radius')
    return [x/norm for x in vector]


def canonical(points, group, phase):
    if phase not in PHASES or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('undeclared methyl phase/nonfinite coordinates')
    indices = [group['parent'], group['axis'], group['frame'], *group['hydrogens']]
    if len(group['hydrogens']) != 3 or len(set(indices)) != 6 or any(type(i) is not int or not 0 <= i < len(points) for i in indices):
        raise ValueError('invalid methyl geometry inventory')
    origin = points[group['parent']]
    axis = unit([x-y for x, y in zip(points[group['axis']], origin)])
    reference = [x-y for x, y in zip(points[group['frame']], origin)]
    projection = sum(x*y for x, y in zip(reference, axis))
    first = unit([x-projection*y for x, y in zip(reference, axis)])
    second = [axis[1]*first[2]-axis[2]*first[1], axis[2]*first[0]-axis[0]*first[2], axis[0]*first[1]-axis[1]*first[0]]
    result = [list(p) for p in points]
    for ordinal, index in enumerate(group['hydrogens']):
        radius = math.dist(points[index], origin)
        if radius <= 1e-12:
            raise ValueError('undefined original C-H radius')
        theta = math.radians(phase+120*ordinal)
        direction = [-a/3+math.sqrt(8/9)*(math.cos(theta)*u+math.sin(theta)*v) for a, u, v in zip(axis, first, second)]
        result[index] = [o+radius*d for o, d in zip(origin, direction)]
    return result


def choose(energies):
    if set(energies) != set(PHASES) or any(not math.isfinite(e) for e in energies.values()):
        raise ValueError('complete finite initial-energy phase grid required')
    minimum = min(energies.values())
    def signed(angle):
        return angle if angle <= 180 else angle-360
    return min((a for a in PHASES if energies[a] <= minimum+1e-8), key=lambda a: (abs(signed(a)), signed(a)))


def initialize(points, groups, selected=False, energy=None):
    """Zero all groups, then one deterministic pass; separate scoring Context."""
    if groups != sorted(groups, key=lambda g: tuple(g['parent_key'])):
        raise ValueError('methyl parent order differs from enrolled lexicographic order')
    indices = [i for group in groups for i in group['hydrogens']]
    if len(indices) != len(set(indices)):
        raise ValueError('overlapping methyl initialization ownership')
    result = [list(p) for p in points]
    # Each construction uses original radii, never a previously rounded candidate.
    for group in groups:
        zero = canonical(points, group, 0)
        for i in group['hydrogens']:
            result[i] = zero[i]
    grids = []
    if selected:
        if energy is None:
            raise ValueError('initial System energy evaluator missing')
        for group in groups:
            values, candidates = {}, {}
            for phase in PHASES:
                own = canonical(points, group, phase)
                candidate = [p[:] for p in result]
                for i in group['hydrogens']:
                    candidate[i] = own[i]
                values[phase] = float(energy(candidate))
                candidates[phase] = [candidate[i][:] for i in group['hydrogens']]
            phase = choose(values)
            for i, point in zip(group['hydrogens'], candidates[phase], strict=True):
                result[i] = point
            grids.append(dict(parent_key=group['parent_key'], selected_phase=phase,
                energies=[dict(phase=a, energy_kj_mol=values[a]) for a in PHASES]))
    return result, grids
