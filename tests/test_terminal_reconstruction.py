"""Terminal construction geometry/custody without MODELLER or Rosetta imports."""
import math
from types import SimpleNamespace

import numpy as np
import pytest
from scripts.hybrid_gap_reconstruction import build_backbone, dihedral
from scripts.terminal_reconstruction_research import prepare_backbone, terminal_torsions

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def case():
    residues = [ResidueIdentity("d", str(i+1), "") for i in range(6)]
    angles = [(-1.1, -.6), (-2., 2.), (-1.4, .8), (-1., -.5), (-1.2, 1.), (-1.4, .7)]
    points = build_backbone([[0., 0., 0.], [1.5, 0., 0.], [2., 1.1, 0.]], angles)
    source, proposal, fixed, generated = {}, {}, [], []
    for i, (r, row) in enumerate(zip(residues, points, strict=True)):
        for name, xyz in {**row, "CB": row["CA"] + [0., 0., 1.5]}.items():
            a = AtomIdentity(r.chain, r.residue_number, r.insertion_code, name)
            line = f"ATOM  {len(proposal)+1:5d} {name:^4s} ALA d{i+1:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           {name[0]}  \n"
            k = r.chain, r.residue_number, r.insertion_code, name
            proposal[k] = {"line": line, "xyz": xyz, "resname": "ALA", "element": name[0]}
            (generated if i < 2 else fixed).append(a)
            if i >= 2:
                source[k] = proposal[k]
    request = SimpleNamespace(target_sequences=[SimpleNamespace(chain="d", sequence="AAAAAA")],
        gaps=[SimpleNamespace(chain="d", left_anchor=None, right_anchor=residues[2], generated_residues=residues[:2],
                             target_interval=SimpleNamespace(start=0, stop=2))],
        sequence_placements=[SimpleNamespace(chain="d", observed_target_indices=list(range(2, 6)), observed_residues=residues[2:])],
        fixed_atoms=fixed, generated_atoms=generated, heterogen_contexts=(), retained_explicit_links=())
    return request, source, proposal, points, residues


def test_reverse_terminal_exact_junction_and_ideal_geometry():
    request, source, proposal, points, _ = case()
    text, payload, report = prepare_backbone(request, source, proposal)
    assert report["boundary_C_N_A"] == pytest.approx(1.33)
    for i, row in enumerate(report["ideal_rows"]):
        for name, xyz in row.items():
            assert xyz == pytest.approx(points[i][name], abs=1e-12)
    last = report["ideal_rows"][-1]
    anchor = points[2]
    assert abs(dihedral(last["CA"], last["C"], anchor["N"], anchor["CA"])) == pytest.approx(math.pi)
    assert all(row["line"] in text for row in source.values())
    assert all(a[3] != "CB" for a in [(line[21], line[22:26].strip(), line[26], line[12:16].strip())
                                     for line in text.splitlines()[:8]])
    assert payload["operation"] == "pack-and-cartmin" and payload["candidates"] == 1
    assert payload["max_iterations"] == 500 and report["source_variables"] == 0


def test_proposal_rigid_frame_cannot_change_source_anchored_reconstruction():
    request, source, proposal, _, _ = case()
    before = prepare_backbone(request, source, proposal)[2]
    rotation = np.array([[0., -1, 0], [1., 0, 0], [0., 0, 1]])
    transformed = {a: {**row, "xyz": np.asarray(row["xyz"]) @ rotation + [15., 20., 10.]}
                   for a, row in proposal.items()}
    after = prepare_backbone(request, source, transformed)[2]
    for a, b in zip(before["ideal_rows"], after["ideal_rows"], strict=True):
        for name in a:
            assert a[name] == pytest.approx(b[name], abs=1e-12)


def test_terminal_phi_one_is_undefined_not_read_from_hidden_reference():
    _, _, proposal, _, order = case()
    torsions, anchor_phi = terminal_torsions(order, proposal, 2)
    assert torsions[0][0] == 0.
    assert torsions[1] == pytest.approx([-2., 2.])
    assert anchor_phi == pytest.approx(-1.4)


def test_changed_proposal_sequence_or_incomplete_axis_fails():
    request, source, proposal, _, _ = case()
    proposal["d", "1", "", "N"]["resname"] = "GLY"
    with pytest.raises(ValueError, match="sequence"):
        prepare_backbone(request, source, proposal)
    request, source, proposal, _, _ = case()
    proposal.pop(("d", "1", "", "CB"))
    with pytest.raises(ValueError, match="atom axis"):
        prepare_backbone(request, source, proposal)
