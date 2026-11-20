"""Private hybrid contracts; do not import the licensed Rosetta runtime."""
import math
from types import SimpleNamespace

import numpy as np
import pytest
from scripts.hybrid_gap_reconstruction import (
    audit_native,
    backbone_jacobian,
    build_backbone,
    build_backbone_from_right,
    diagnostic_counts,
    dihedral,
    extend,
    fit_route,
    prepare_payload,
    run,
    validate_inputs,
)
from scripts.hybrid_gap_rosetta import carbonyl_oxygen, closure_xml

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def test_internal_coordinate_roundtrip():
    a, b, c = np.array([0., 0., 0.]), np.array([1., 0., 0.]), np.array([1., 1., 0.])
    for torsion in (-2.4, -.9, 0., 1.3, math.pi):
        d = extend(a, b, c, 1.33, 1.9, torsion)
        assert np.linalg.norm(d - c) == pytest.approx(1.33)
        error = (dihedral(a, b, c, d) - torsion + math.pi) % (2 * math.pi) - math.pi
        assert error == pytest.approx(0., abs=1e-12)


def test_build_ideal_backbone_with_trans_peptides():
    torsions = [(-1.1, -.6), (-2., 2.), (-1.4, .8), (-1., -.5)]
    rows = build_backbone([[0., 0., 0.], [1.5, 0., 0.], [2., 1.1, 0.]], torsions)
    for i, row in enumerate(rows):
        assert np.linalg.norm(row["N"] - row["CA"]) == pytest.approx(1.46)
        assert np.linalg.norm(row["CA"] - row["C"]) == pytest.approx(1.525)
        if i + 1 < len(rows):
            following = rows[i + 1]
            assert np.linalg.norm(row["C"] - following["N"]) == pytest.approx(1.33)
            assert abs(dihedral(row["CA"], row["C"], following["N"], following["CA"])) == pytest.approx(math.pi)
            assert dihedral(row["N"], row["CA"], row["C"], following["N"]) == pytest.approx(torsions[i][1])


def test_degenerate_backbone_is_not_silently_replaced():
    with pytest.raises(ValueError, match="degenerate"):
        dihedral([0, 0, 0], [0, 0, 0], [0, 0, 0], [1, 0, 0])
    with pytest.raises(ValueError, match="degenerate"):
        extend([0, 0, 0], [1, 0, 0], [2, 0, 0], 1., 1., 1.)


def test_reverse_growth_preserves_right_junction_and_torsions():
    angles = [(-1.1, -.6), (-2., 2.), (-1.4, .8), (-1., -.5)]
    forward = build_backbone([[0., 0., 0.], [1.5, 0., 0.], [2., 1.1, 0.]], [*angles, (-1.2, 1.)])
    right = forward[-1]
    reverse = build_backbone_from_right([right[n] for n in ("N", "CA", "C")], angles, -1.2)
    for a, b in zip(forward, reverse):
        for name in ("N", "CA", "C", "O"):
            assert a[name] == pytest.approx(b[name], abs=1e-12)


def test_backbone_jacobian_matches_finite_differences():
    left = [[0., 0., 0.], [1.5, 0., 0.], [2., 1.1, 0.]]
    q = np.array([[-1.1, -.6], [-2., 2.], [-1.4, .8], [-1., -.5]])
    rows = build_backbone(left, q)
    jac = backbone_jacobian(rows)
    for column in range(q.size):
        plus, minus = q.copy(), q.copy()
        plus.flat[column] += 1e-6
        minus.flat[column] -= 1e-6
        x, y = build_backbone(left, plus), build_backbone(left, minus)
        fd = (np.array([[r[n] for n in ("N", "CA", "C", "O")] for r in x]) -
              np.array([[r[n] for n in ("N", "CA", "C", "O")] for r in y])) / 2e-6
        assert jac[..., column] == pytest.approx(fd, abs=5e-8)


def test_route_fit_uses_only_saved_ca_path_and_is_bounded():
    left = [[0., 0., 0.], [1.5, 0., 0.], [2., 1.1, 0.]]
    angles = [(-1.1, -.6), (-2., 2.), (-1.4, .8), (-1., -.5)]
    forward = build_backbone(left, [*angles, (-1.2, 1.)])
    right = [forward[-1][n] for n in ("N", "CA", "C")]
    ca = [r["CA"] for r in forward[:-1]]
    initial = np.array(angles) + .05
    _, report = fit_route(left, right, initial, -1.2, ca, max_nfev=30)
    assert report["evaluations"] <= 30
    assert report["ca_rms_to_diffusion_A"] < .1
    assert report["reference_coordinates_read"] is False
    assert report["observed_variables"] == 0


def test_native_audit_catches_motion_before_copying_and_case_collision():
    fixed = ["D", "82", "A", "CA"]
    new = ["d", "82", "A", "CA"]
    payload = {"fixed": [{"identity": fixed, "xyz": [0, 0, 0]}], "generated": [new]}
    native = {"atoms": [{"identity": fixed, "xyz": [0, 0, 0]}, {"identity": new, "xyz": [1, 2, 3]}]}
    assert audit_native(payload, native)[1] == 0
    native["atoms"][0]["xyz"] = [1e-6, 0, 0]
    with pytest.raises(ValueError, match="before copying"):
        audit_native(payload, native)
    native["atoms"][0]["xyz"] = [0, 0, 0]
    native["atoms"][1]["identity"] = fixed
    with pytest.raises(ValueError, match="axis"):
        audit_native(payload, native)


def test_changed_source_refused_before_any_worker():
    atom = SimpleNamespace(chain="d", residue_number="82", insertion_code="A", atom_name="CA")
    request = SimpleNamespace(fixed_atoms=[atom], generated_atoms=[])
    source = {("d", "82", "A", "CA"): {"line": "original"}}
    proposal = {("d", "82", "A", "CA"): {"line": "changed"}}
    with pytest.raises(ValueError, match="observed"):
        validate_inputs(request, source, proposal)


def test_closure_is_bounded_and_packs_before_selection_without_design():
    xml = closure_xml({"start": 10, "stop": 13, "cut": 11}, range(10, 14), 20, 32)
    assert 'closure_attempts="32"' in xml
    assert 'pre_selection_mover="pack"' in xml
    assert '<RestrictToRepacking' in xml and '<PreventRepackingRLT' in xml
    assert 'resnums="10,11,12,13"' in xml
    assert 'torsion="180"' in xml
    assert 'randomize_alpha_backbone_by_rama' not in xml


def test_readiness_error_not_hidden_or_mislabelled_as_loop_geometry():
    counts = diagnostic_counts([{"severity": "ERROR", "category": "missing_terminals"},
                                {"severity": "ERROR", "category": "non_planar_amide"}])
    assert counts == {"geometry_errors": 1, "other_errors": 1,
                      "all_error_categories": ["missing_terminals", "non_planar_amide"]}


def synthetic_case():
    """Known geometry is an engineering control, not a model accuracy result."""
    residues = [ResidueIdentity("d", str(i + 80), "A" if i == 2 else "") for i in range(5)]
    rows = build_backbone([[0., 0., 0.], [1.5, 0., 0.], [2., 1.1, 0.]], [(-1.1, -.6)] * 5)
    source, proposal, fixed, generated = {}, {}, [], []
    serial = 0
    for i, (residue, points) in enumerate(zip(residues, rows, strict=True)):
        for name, xyz in {**points, "CB": points["CA"] + [0., 0., 1.5]}.items():
            serial += 1
            a = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, name)
            line = (f"ATOM  {serial:5d} {name:^4s} ALA {residue.chain}{int(residue.residue_number):4d}"
                    f"{residue.insertion_code or ' '}   {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00           {name[0]}  \n")
            k = (residue.chain, residue.residue_number, residue.insertion_code, name)
            proposal[k] = {"line": line, "xyz": xyz, "resname": "ALA", "element": name[0]}
            if i in (0, 4):
                source[k] = proposal[k]
                fixed.append(a)
            else:
                generated.append(a)
    gap = SimpleNamespace(chain="d", generated_residues=residues[1:4], left_anchor=residues[0], right_anchor=residues[4],
                          target_interval=SimpleNamespace(start=1, stop=4))
    request = SimpleNamespace(fixed_atoms=fixed, generated_atoms=generated, gaps=[gap],
        target_sequences=[SimpleNamespace(chain="d", sequence="AAAAA")],
        sequence_placements=[SimpleNamespace(chain="d", observed_target_indices=[0, 4], observed_residues=[residues[0], residues[4]])])
    return request, source, proposal


def test_prepare_keeps_observations_and_insertion_codes_and_removes_only_new_sides():
    request, source, proposal = synthetic_case()
    text, payload = prepare_payload(request, source, proposal)
    assert all(row["line"] in text for row in source.values())
    assert len(payload["generated"]) == 15
    assert {tuple(row) for row in payload["generated"]} == {
        (a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    assert payload["residues"][2]["insertion_code"] == "A"
    assert payload["gaps"][0]["cut"] == 2
    generated_side_lines = [r["line"] for k, r in proposal.items() if k not in source and k[3] == "CB"]
    assert not any(line in text for line in generated_side_lines)


def test_partial_residue_and_terminal_masks_are_refused():
    request, source, proposal = synthetic_case()
    a = request.generated_atoms.pop()
    request.fixed_atoms.append(a)
    k = a.chain, a.residue_number, a.insertion_code, a.atom_name
    source[k] = proposal[k]
    with pytest.raises(ValueError, match="partial-residue"):
        validate_inputs(request, source, proposal)
    request, source, proposal = synthetic_case()
    request.gaps[0].left_anchor = None
    with pytest.raises(ValueError, match="internal gaps"):
        validate_inputs(request, source, proposal)


def test_existing_output_is_not_adopted_or_modified(tmp_path):
    sentinel = tmp_path / "user-work.txt"
    sentinel.write_text("preserve")
    args = SimpleNamespace(output=tmp_path)
    with pytest.raises(FileExistsError):
        run(args)
    assert not hasattr(args, "owned_output")
    assert list(tmp_path.iterdir()) == [sentinel]


def test_postclosure_carbonyl_is_planar_without_moving_anchors():
    alpha, carbon = [0., 0., 0.], [1.525, 0., 0.]
    next_n = extend([0, 1, 0], alpha, carbon, 1.33, math.radians(116.2), .7).tolist()
    before = [alpha.copy(), carbon.copy(), next_n.copy()]
    oxygen = np.array(carbonyl_oxygen(alpha, carbon, next_n))
    assert np.linalg.norm(oxygen - carbon) == pytest.approx(1.23)
    u, v = oxygen - carbon, np.array(next_n) - carbon
    angle = math.degrees(math.acos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))))
    assert angle == pytest.approx(123.)
    assert abs(dihedral(alpha, oxygen, carbon, next_n)) == pytest.approx(math.pi)
    assert before == [alpha, carbon, next_n]
