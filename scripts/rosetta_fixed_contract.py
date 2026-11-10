"""Research-only atom-fixed scoring contract; no Rosetta imports or new gate.

The source graph is deliberately partial outside the target. Missing short paths
are labelled unresolved, never asserted to be proven nonbonded interactions.
"""
from __future__ import annotations

import math
from collections import deque

from scripts.research_rosetta_repack import identity


def indexed(rows: list[dict]) -> dict:
    result = {}
    for row in rows:
        key = identity(row)
        if key in result:
            raise ValueError('duplicate snapshot atom identity')
        xyz = row['xyz']
        if len(xyz) != 3 or any(isinstance(v, bool) or not isinstance(v, (float, int)) or not math.isfinite(v) for v in xyz):
            raise ValueError('nonfinite or malformed coordinates')
        result[key] = row
    return result


def pair_key(a, b) -> tuple:
    a, b = tuple(a), tuple(b)
    if a == b:
        raise ValueError('self bond')
    return tuple(sorted((a, b)))


def import_audit(expected: list[dict], actual: list[dict], expected_external: list, actual_external: list,
                 target_internal: list, actual_target_internal: list, allowed_extra: list) -> dict:
    """Exact heavy-coordinate/identity and independently declared topology checks."""
    before, after = indexed(expected), indexed(actual)
    missing = sorted(before.keys()-after.keys())
    moved = sorted(k for k in before.keys() & after.keys() if before[k]['xyz'] != after[k]['xyz'])
    extras = sorted(after.keys()-before.keys())
    disallowed = sorted(set(extras)-{tuple(k) for k in allowed_extra})
    external = {pair_key(*p) for p in expected_external}
    loaded_external = {pair_key(*p) for p in actual_external}
    internal = {pair_key(*p) for p in target_internal}
    loaded_internal = {pair_key(*p) for p in actual_target_internal}
    return dict(accepted=not (missing or moved or disallowed or external != loaded_external or internal != loaded_internal),
                missing=missing, moved=moved, scorer_owned_extras=extras, disallowed_extras=disallowed,
                missing_external=sorted(external-loaded_external), unexpected_external=sorted(loaded_external-external),
                missing_target_internal=sorted(internal-loaded_internal), unexpected_target_internal=sorted(loaded_internal-internal))


def short_paths(start: tuple, edges: list, *, max_bonds: int = 3) -> dict:
    graph = {}
    for a, b in edges:
        a, b = pair_key(a, b)
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
    distances = {start: 0}
    queue = deque([start])
    while queue:
        a = queue.popleft()
        if distances[a] == max_bonds:
            continue
        for b in sorted(graph.get(a, ())):
            if b not in distances:
                distances[b] = distances[a]+1
                queue.append(b)
    return distances


def contacts(rows: list[dict], generated: list, edges: list, radii: dict[str, float]) -> list[dict]:
    """Replay the old -0.6 Å clearance diagnostic with an explicit path label.

This is descriptive only: it neither returns acceptance nor removes a finding.
Bonded and angle-connected pairs match the existing comparator's exclusions.
"""
    atoms = indexed(rows)
    declared_generated = [tuple(a) for a in generated]
    generated = set(declared_generated)
    if len(generated) != len(declared_generated):
        raise ValueError('duplicate generated identity')
    if not generated <= atoms.keys():
        raise ValueError('generated identity absent')
    if any(tuple(a) not in atoms or tuple(b) not in atoms for a, b in edges):
        raise ValueError('graph endpoint absent')
    checked, findings = set(), []
    for a in sorted(generated):
        paths = short_paths(a, edges)
        for b in sorted(atoms):
            pair = pair_key(a, b) if a != b else None
            if paths.get(b, 4) < 3 or pair in checked:
                continue
            checked.add(pair)
            length = math.dist(atoms[a]['xyz'], atoms[b]['xyz'])
            clearance = length-radii[atoms[a]['element']]-radii[atoms[b]['element']]
            if clearance < -.6:
                findings.append(dict(left=a, right=b, distance_angstrom=length, clearance_angstrom=clearance,
                                     partial_graph_bond_distance=paths.get(b),
                                     topology='known_1_4' if paths.get(b) == 3 else 'no_short_path_in_partial_graph',
                                     both_generated=b in generated, old_finding_retained=True))
    return findings


def choose_scored_own(rows: list[dict]) -> dict:
    """Choose only within the complete old-gate-feasible own pool; no truth input."""
    indices = [r['index'] for r in rows]
    if len(set(indices)) != len(indices):
        raise ValueError('duplicate candidate index')
    pool = [r for r in rows if r['old_source_accepted'] and r['common_accepted']]
    if not pool:
        return dict(status='no_old_gate_feasible_candidate', selected_index=None)
    if any(r['score_status'] != 'scored' or not isinstance(r.get('total_score'), (float, int))
           or isinstance(r['total_score'], bool) or not math.isfinite(r['total_score']) for r in pool):
        return dict(status='selection_blocked_incomplete_scores', selected_index=None)
    winner = min(pool, key=lambda r: (r['total_score'], r['index']))
    return dict(status='selected', selected_index=winner['index'], total_score=winner['total_score'])
