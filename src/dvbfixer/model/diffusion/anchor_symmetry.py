"""Research-only heavy-atom label sensitivity from canonical FF graph symmetry.

No observations are renamed or moved. The exact-name metric remains mandatory;
the minimum possible maximum under allowed graph automorphisms is additional
analysis, not a production movement exception or a general chemistry oracle.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

REVISION = "forcefield-heavy-graph-label-sensitivity-v1"
_FIXED_CORE = frozenset({"N", "CA", "C", "O", "OXT", "CB"})


@dataclass(frozen=True)
class SymmetryMeasurement:
    exact_maximum_angstrom: float
    symmetry_minimax_angstrom: float
    mobile_name_for_reference: tuple[tuple[str, str], ...]
    automorphism_count: int
    forcefield_template: str
    graph_sha256: str
    revision: str = REVISION


@lru_cache(maxsize=1)
def _forcefield():
    from openmm.app import ForceField

    from dvbfixer.ffutils import FF_ALIASES

    return ForceField(*FF_ALIASES["amber"])


@lru_cache(maxsize=256)
def residue_label_permutations(resname: str, atom_names: tuple[str, ...]):
    """Conservative exact atom-type/charge/H-degree graph automorphisms.

    Protonation-ambiguous HIS is not guessed into HIE/HID. Unknown templates,
    incomplete axes and excessive symmetry return an explicit failure, never
    a falsely improved zero. Backbone and CB names are fixed throughout.
    """
    if not atom_names or len(set(atom_names)) != len(atom_names):
        raise ValueError("symmetry analysis requires a unique complete heavy axis")
    templates = _forcefield()._templates
    template = None
    for key in (resname, "N" + resname, "C" + resname):
        candidate = templates.get(key)
        if candidate is not None and {a.name for a in candidate.atoms
                if a.element is not None and a.element.symbol not in {"H", "D"}} == set(atom_names):
            template = candidate
            break
    if template is None:
        raise ValueError("no unambiguous force-field template matches this heavy axis")
    neighbors = {i: set() for i in range(len(template.atoms))}
    for first, second in template.bonds:
        neighbors[first].add(second)
        neighbors[second].add(first)
    heavy = {i: a for i, a in enumerate(template.atoms) if a.name in atom_names}
    adjacency = {a.name: {template.atoms[j].name for j in neighbors[i] if j in heavy} for i, a in heavy.items()}
    labels = {a.name: (a.element.symbol, a.type, a.parameters.get("charge"),
        sum(template.atoms[j].element is not None and template.atoms[j].element.symbol in {"H", "D"} for j in neighbors[i]),
        a.name if a.name in _FIXED_CORE else "") for i, a in heavy.items()}
    allowed = {name: sorted(other for other in atom_names if labels[other] == labels[name]
                            and len(adjacency[other]) == len(adjacency[name])) for name in atom_names}
    order = sorted(atom_names, key=lambda name: (len(allowed[name]), -len(adjacency[name]), name))
    permutations = []

    def visit(mapping, used):
        if len(mapping) == len(order):
            permutations.append(tuple((name, mapping[name]) for name in sorted(atom_names)))
            if len(permutations) > 64:
                raise ValueError("symmetry analysis exceeds the declared automorphism bound")
            return
        name = order[len(mapping)]
        for other in allowed[name]:
            if other in used or any((old in adjacency[name]) != (new in adjacency[other]) for old, new in mapping.items()):
                continue
            mapping[name] = other
            visit(mapping, used | {other})
            del mapping[name]

    visit({}, set())
    if not permutations:
        raise ValueError("force-field heavy graph has no valid identity mapping")
    graph = {"template": template.name, "labels": labels,
             "neighbors": {name: sorted(value) for name, value in adjacency.items()}, "fixed_core": sorted(_FIXED_CORE)}
    graph_sha = hashlib.sha256(json.dumps(graph, sort_keys=True, allow_nan=False).encode()).hexdigest()
    return template.name, graph_sha, tuple(sorted(permutations))


def measure_label_sensitivity(resname: str, reference, mobile) -> SymmetryMeasurement:
    """Coordinates must already share a fitting frame; caller inputs are borrowed."""
    if set(reference) != set(mobile):
        raise ValueError("symmetry analysis requires identical complete heavy axes")
    names = tuple(sorted(reference))
    a = np.asarray([reference[name] for name in names], dtype=float)
    b = np.asarray([mobile[name] for name in names], dtype=float)
    if a.shape != (len(names), 3) or b.shape != a.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("symmetry analysis requires finite three-dimensional coordinates")
    template, graph_sha, permutations = residue_label_permutations(resname, names)
    exact = float(np.linalg.norm(a - b, axis=1).max())
    scores = [(max(float(np.linalg.norm(np.asarray(reference[name]) - np.asarray(mobile[other])))
                   for name, other in permutation), permutation) for permutation in permutations]
    maximum, selected = min(scores)
    return SymmetryMeasurement(exact, maximum, selected, len(permutations), template, graph_sha)
