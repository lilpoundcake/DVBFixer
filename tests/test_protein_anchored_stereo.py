"""Independent observed-partner/RDKit checks, not empirical rotamer tests."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from scripts.branched_stereo_oracle import rdkit_labels

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.anchored_branch import complete_anchored_branch
from dvbfixer.protein_completion.anchored_source import (
    materialize_anchored_source,
    plan_anchored_source,
)
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.geometry import add, cross, dot, scaled, sub, unit


@pytest.mark.parametrize("name,number,missing", (("VAL", 5, "CG1"), ("VAL", 5, "CG2"),
                                                ("THR", 7, "OG1"), ("THR", 7, "CG2")))
def test_missing_either_branch_retains_correct_independent_stereo(name, number, missing):
    root = Path(__file__).resolve().parents[1]
    source = (root/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    target = ComponentInstanceRef(1, "A", str(number))
    masked = b"".join(l for l in source.splitlines(keepends=True) if l.startswith(b"ATOM  ")
                      and int(l[22:26]) <= 72 and not (int(l[22:26]) == number and l[12:16].strip().decode() == missing))
    atoms = _parse_atom_lines(masked.decode().splitlines(keepends=True))
    observations = tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position) for a in atoms)
    req = CompletionRequest(hashlib.sha256(masked).hexdigest(), target, name, observations, (missing,))
    result = complete_anchored_branch(req)
    assert len(result.attempts) == 1
    candidate = result.attempts[0].candidate
    points = {o.atom.name: o.xyz for o in observations if o.atom.component == target}
    points[missing] = candidate.additions[0].xyz
    assert rdkit_labels(points, name) == ({"CA": "S", "CB": "R"} if name == "THR" else {"CA": "S"})
    if name == "THR":
        branch = BRANCHES[name]
        partner = branch.second if missing == branch.first else branch.first
        cb = points["CB"]
        normal = unit(cross(sub(points["CA"], cb), sub(points[partner], cb)))
        tip = points[missing]
        mirrored = add(tip, scaled(normal, -2*dot(sub(tip, cb), normal)))
        assert rdkit_labels(points | {missing: mirrored}, name) == {"CA": "S", "CB": "S"}


@pytest.mark.parametrize("name,number,missing", (("VAL", 5, "CG1"), ("VAL", 5, "CG2"),
                                                ("THR", 7, "OG1"), ("THR", 7, "CG2")))
def test_rounded_singleton_source_has_independent_cip(name, number, missing):
    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    src = b"".join(l for l in original.splitlines(keepends=True) if l.startswith(b"ATOM  ")
                   and int(l[22:26]) <= 72 and not (int(l[22:26]) == number and l[12:16].strip().decode() == missing))
    src += b"TER\nEND\n"
    target = ComponentInstanceRef(1, "A", str(number))
    plan = plan_anchored_source(src, target)
    out = materialize_anchored_source(src, plan, complete_anchored_branch(plan.request).selected)
    points = {a.name: a.position for a in _parse_atom_lines(out.pdb.decode().splitlines(keepends=True))
              if _component_ref(a) == target}
    assert rdkit_labels(points, name) == ({"CA": "S", "CB": "R"} if name == "THR" else {"CA": "S"})
