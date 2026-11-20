"""Private IDPForge transport/compatibility tests without scientific model imports."""
from types import SimpleNamespace

import numpy as np
import pytest
from deploy.idpforge.prepare_runtime import patch_text
from deploy.idpforge.research_worker import source_secondary_structure
from scripts.audit_idpforge_research import register_ca
from scripts.idpforge_research import admit

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def case(chain="d"):
    residues = [ResidueIdentity(chain, str(i+1), "") for i in range(6)]
    fixed, generated, source = [], [], {}
    for i, r in enumerate(residues):
        for name in CANONICAL_HEAVY_ATOMS["ALA"]:
            a = AtomIdentity(chain, r.residue_number, "", name)
            if i < 2:
                generated.append(a)
            else:
                fixed.append(a)
                line = f"ATOM  {len(fixed):5d} {name:^4s} ALA {chain}{i+1:4d}    {1.:8.3f}{2.:8.3f}{3.:8.3f}  1.00  0.00           {name[0]}  \n"
                source[chain, r.residue_number, "", name] = {"line": line, "xyz": [1., 2., 3.], "resname": "ALA"}
    request = SimpleNamespace(target_sequences=[SimpleNamespace(chain=chain, sequence="AAAAAA")],
        gaps=[SimpleNamespace(chain=chain, left_anchor=None, right_anchor=residues[2],
              generated_residues=residues[:2], target_interval=SimpleNamespace(start=0, stop=2))],
        sequence_placements=[SimpleNamespace(chain=chain, observed_target_indices=[2, 3, 4, 5],
                                            observed_residues=residues[2:])],
        heterogen_contexts=(), retained_explicit_links=(), fixed_atoms=fixed, generated_atoms=generated)
    return request, source


def test_known_sequence_source_only_and_case_sensitive_terminal_axis():
    request, source = case()
    result = admit(request, source)
    assert result["chain"] == "d"
    assert result["sequence"] == "AAAAAA"
    assert result["gap_length"] == 2
    assert result["fixed_atoms"] == 20
    assert result["generated_atoms"] == 10
    assert not result["cropped"] and not result["reference_coordinates_read"]
    assert [r["generated"] for r in result["residues"]] == [True, True, False, False, False, False]


def test_partial_source_and_generated_ownership_refused():
    request, source = case()
    a = request.fixed_atoms.pop()
    del source[a.chain, a.residue_number, a.insertion_code, a.atom_name]
    with pytest.raises(ValueError, match="complete canonical"):
        admit(request, source)
    request, source = case()
    request.generated_atoms.pop()
    with pytest.raises(ValueError, match="whole-gap ownership"):
        admit(request, source)


def test_cannot_feed_internal_or_multichain_context_to_terminal_control():
    request, source = case()
    request.gaps[0].left_anchor = request.gaps[0].right_anchor
    with pytest.raises(ValueError, match="N-terminal gap"):
        admit(request, source)
    request, source = case()
    request.target_sequences.append(SimpleNamespace(chain="D", sequence="AAA"))
    with pytest.raises(ValueError, match="one chain"):
        admit(request, source)


def test_unsupported_context_and_partial_occupancy_refused():
    request, source = case()
    request.retained_explicit_links = ["LINK"]
    with pytest.raises(ValueError, match="explicit links"):
        admit(request, source)
    request, source = case()
    row = next(iter(source.values()))
    row["line"] = row["line"][:54] + "  0.50" + row["line"][60:]
    with pytest.raises(ValueError, match="full-occupancy"):
        admit(request, source)


@pytest.mark.parametrize("number,icode", [("3", "A"), ("99", "")])
def test_native_parser_insertion_codes_and_discontinuous_numbers_fail_closed(number, icode):
    request, source = case()
    request.sequence_placements[0].observed_residues[0] = ResidueIdentity("d", number, icode)
    request.gaps[0].right_anchor = request.sequence_placements[0].observed_residues[0]
    with pytest.raises(ValueError, match="contiguous numbers"):
        admit(request, source)


def test_compatibility_only_modifies_expected_import_and_dataclass():
    text = 'import torch\nattn_core_inplace_cuda = importlib.import_module("attn_core_inplace_cuda")\ndef forward(): return 1\n'
    patched = patch_text("openfold/model/structure_module.py", text)
    assert "if torch.cuda.is_available() else None" in patched
    assert patched.endswith("def forward(): return 1\n")
    with pytest.raises(ValueError, match="unexpected pinned source"):
        patch_text("openfold/model/structure_module.py", "wrong source")
    assert patch_text("idpforge/model.py", "unchanged recon") == "unchanged recon"


def test_registration_is_one_proper_rigid_transform_for_all_atoms():
    moving = np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1]])
    rotation = np.array([[0., -1, 0], [1., 0, 0], [0., 0, 1]])
    target = moving @ rotation + [3., 4., 5.]
    xyz = np.stack((moving, moving + [2., 0, 0]), axis=1)
    result, report = register_ca(xyz, moving, target)
    assert np.allclose(result, xyz @ rotation + [3., 4., 5.])
    assert report["observed_CA_rms_A"] < 1e-10
    assert np.linalg.det(report["rotation_row_vector"]) == pytest.approx(1.)
    assert not report["reflection"] and not report["per_residue_fit"]


def test_registration_does_not_hide_nonrigid_motion_or_reflect_mirror():
    moving = np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1]])
    target = moving.copy()
    target[0] += .5
    _, report = register_ca(moving, moving, target)
    assert report["observed_CA_rms_A"] > .1
    _, report = register_ca(moving, moving, moving * [-1., 1, 1])
    assert np.linalg.det(report["rotation_row_vector"]) == pytest.approx(1.)
    assert report["observed_CA_rms_A"] > .1
    with pytest.raises(ValueError, match="degenerate"):
        register_ca(moving, np.zeros((4, 3)), target)


def test_secondary_structure_passes_radians_without_double_conversion():
    md = SimpleNamespace(load=lambda path: path,
        compute_dssp=lambda traj, simplified: np.array([["C", "H", "E"]]),
        compute_phi=lambda traj: (None, np.array([[.5, -1.]])),
        compute_psi=lambda traj: (None, np.array([[1., -.5]])))

    def assign(angles):
        assert np.array_equal(angles[:, 0], [-np.pi, .5, -1.])
        assert np.array_equal(angles[:, 1], [1., -.5, np.pi])
        return "BPL"

    assert source_secondary_structure(md, assign, "source-only.pdb", 2) == "CCBHE"
