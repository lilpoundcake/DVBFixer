"""Private heavy-compatible carbon-CH2 initial geometry, not production repair.

Each new H makes tetrahedral angles with both fixed heavy bonds. Their original
heavy-heavy angle is preserved; the H-H angle is implied, not independently
forced. The first lexicographic H takes +normal of lexicographic heavy neighbors.
"""
from __future__ import annotations

import math

from scripts.graph_methyl_geometry import inventory, unit

from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def discover(rows, bonds, newly_added, mobile):
    adjacent = inventory(rows, bonds)
    for selected in (newly_added, mobile):
        if any(type(i) is not int or not 0 <= i < len(rows) for i in selected):
            raise ValueError('invalid CH2 ownership index')
    if any(rows[i]['element'] not in {'H', 'D'} for i in newly_added):
        raise ValueError('new H ownership includes a heavy atom')
    groups, exclusions = [], []
    for parent in sorted(range(len(rows)), key=lambda i: tuple(rows[i]['key'])):
        row = rows[parent]
        if row['element'] != 'C':
            continue
        hydrogens = sorted((i for i in adjacent[parent] if rows[i]['element'] in {'H', 'D'}), key=lambda i: tuple(rows[i]['key']))
        heavy = sorted(adjacent[parent]-set(hydrogens), key=lambda i: tuple(rows[i]['key']))
        if len(hydrogens) != 2 or len(heavy) != 2:
            continue
        if not set(hydrogens) <= set(newly_added) & set(mobile):
            exclusions.append(dict(parent_key=row['key'], reason='measured/kept/nonmobile H; no rewrite'))
            continue
        name = VARIANT_TO_PARENT.get(row['key'][3], row['key'][3])
        if (name not in CANONICAL_HEAVY_ATOMS or row['key'][-1] not in CANONICAL_HEAVY_ATOMS[name]
                or any(rows[i]['key'][:4] != row['key'][:4] for i in hydrogens+heavy)
                or any(adjacent[i] != {parent} for i in hydrogens)):
            raise ValueError('CH2 outside canonical intra-residue context')
        groups.append(dict(parent=parent, parent_key=row['key'][:], hydrogens=hydrogens, heavy=heavy))
    return dict(groups=groups, exclusions=exclusions)


def canonical(points, group):
    if any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('invalid CH2 coordinates')
    indices = [group['parent'], *group['heavy'], *group['hydrogens']]
    if len(group['heavy']) != 2 or len(group['hydrogens']) != 2 or len(set(indices)) != 5 or any(type(i) is not int or not 0 <= i < len(points) for i in indices):
        raise ValueError('invalid CH2 geometry inventory')
    origin = points[group['parent']]
    first, second = [unit([x-y for x, y in zip(points[i], origin)]) for i in group['heavy']]
    normal = unit([first[1]*second[2]-first[2]*second[1], first[2]*second[0]-first[0]*second[2], first[0]*second[1]-first[1]*second[0]])
    bisector = unit([x+y for x, y in zip(first, second)])
    cosine = sum(x*y for x, y in zip(first, bisector))
    if cosine <= 1/3:
        raise ValueError('fixed heavy angle cannot support distinct tetrahedral C-H-heavy slots')
    axial = -1/(3*cosine)
    perpendicular = math.sqrt(1-axial*axial)
    if perpendicular <= 1e-8:
        raise ValueError('degenerate CH2 H-H separation')
    result = [list(p) for p in points]
    for sign, index in zip((1., -1.), group['hydrogens'], strict=True):
        radius = math.dist(points[index], origin)
        if radius <= 1e-12:
            raise ValueError('undefined original CH2 C-H radius')
        direction = [axial*b+sign*perpendicular*n for b, n in zip(bisector, normal)]
        result[index] = [o+radius*d for o, d in zip(origin, direction)]
    return result


def initialize(points, groups):
    if groups != sorted(groups, key=lambda g: tuple(g['parent_key'])):
        raise ValueError('CH2 parent identity order differs')
    owned = [i for group in groups for i in group['hydrogens']]
    if len(owned) != len(set(owned)):
        raise ValueError('overlapping CH2 ownership')
    result = [list(p) for p in points]
    for group in groups:
        local = canonical(points, group)
        for i in group['hydrogens']:
            result[i] = local[i]
    return result
