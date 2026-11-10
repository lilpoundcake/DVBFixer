"""Independent RDKit graph/CIP controls, run in the scientific environment."""
from __future__ import annotations

from pathlib import Path

import pytest
from scripts.branched_stereo_oracle import rdkit_labels

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom, _parse_atom_lines
from dvbfixer.protein_completion.branched import BRANCHES, named_beta_volume
from dvbfixer.protein_completion.branched_source import (
    materialize_branched_source,
    plan_branched_source,
)
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.geometry import add, cross, dot, scaled, sub, unit
from dvbfixer.protein_completion.service import complete_sidechain


@pytest.mark.parametrize("resid,name", ((5, "VAL"), (17, "VAL"), (26, "VAL"), (70, "VAL"),
    (7, "THR"), (9, "THR"), (12, "THR"), (14, "THR"), (22, "THR"), (55, "THR"), (66, "THR")))
def test_real_reference_named_volume_agrees_with_independent_cip(resid, name):
    lines = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1UBQ.pdb").read_text().splitlines()
    points = {l[12:16].strip(): tuple(float(l[i:i+8]) for i in (30, 38, 46)) for l in lines
              if l.startswith("ATOM  ") and l[21] == "A" and int(l[22:26]) == resid}
    labels = rdkit_labels(points, name)
    assert labels == ({"CA": "S", "CB": "R"} if name == "THR" else {"CA": "S"})
    branch = BRANCHES[name]
    assert branch.expected_named_volume_sign*named_beta_volume(points["CA"], points["CB"],
                                                              points[branch.first], points[branch.second]) > .5


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_every_synthetic_chi_proposal_has_independent_stereochemistry(name):
    target = ComponentInstanceRef(1, "A", "1")
    points = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
    branch = BRANCHES[name]
    request = CompletionRequest("a"*64, target, name,
        tuple(Observation(ExactAtomRef(target, n), n[0], p) for n, p in points.items()),
        ("CB", branch.first, branch.second))
    for attempt in complete_sidechain(request).attempts:
        proposal = points | {a.atom.name: a.xyz for a in attempt.candidate.additions}
        assert rdkit_labels(proposal, name) == ({"CA": "S", "CB": "R"} if name == "THR" else {"CA": "S"})
        if name == "THR":
            cb = proposal["CB"]
            normal = unit(cross(sub(proposal["CA"], cb), sub(proposal["OG1"], cb)))
            mirrored = add(proposal["CG2"], scaled(normal, -2*dot(sub(proposal["CG2"], cb), normal)))
            assert rdkit_labels(proposal | {"CG2": mirrored}, name) == {"CA": "S", "CB": "S"}


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_all_three_rounded_source_proposals_have_independent_cip(name):
    points = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
    lines = []
    for i, (atom_name, xyz) in enumerate(points.items(), 1):
        line = _format_generated_atom(serial=i, name=atom_name, element=atom_name[0], position=xyz,
                                      residue_name=name, chain_id="A", sequence_number="1", insertion_code="")
        lines.append("ATOM  "+line[6:])
    src = ("".join(lines)+"TER\nEND\n").encode()
    plan = plan_branched_source(src, ComponentInstanceRef(1, "A", "1"))
    for attempt in complete_sidechain(plan.request).attempts:
        output = materialize_branched_source(src, plan, attempt.candidate)
        readback = {a.name: a.position for a in _parse_atom_lines(output.pdb.decode().splitlines(keepends=True))}
        assert rdkit_labels(readback, name) == ({"CA": "S", "CB": "R"} if name == "THR" else {"CA": "S"})
