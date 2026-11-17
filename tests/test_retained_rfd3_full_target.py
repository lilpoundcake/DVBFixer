import subprocess
import sys
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest
from scripts.preflight_retained_rfd3_full_target import check_array, conditioning

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def fixture(chain="A"):
    residues = tuple(ResidueIdentity(chain, str(i), "") for i in (1, 2, 3))
    fixed = [AtomIdentity(chain, str(i), "", name) for i in (1, 3) for name in CANONICAL_HEAVY_ATOMS["CYS"]]
    oxt = AtomIdentity(chain, "3", "", "OXT")
    fixed.append(oxt)
    restored = [AtomIdentity(chain, "2", "", name) for name in CANONICAL_HEAVY_ATOMS["GLY"]]
    lines = []
    for serial, atom in enumerate(fixed, 1):
        field = f" {atom.atom_name:<3}"
        lines.append(f"ATOM  {serial:5d} {field} CYS {chain}{int(atom.residue_number):4d}    "
                     f"{float(serial):8.3f}{2.:8.3f}{3.:8.3f}{1:6.2f}{20.:6.2f}          {atom.atom_name[0]:>2s}\n")
    request = SimpleNamespace(target_sequences=(SimpleNamespace(chain=chain, sequence="CGC"),),
        sequence_placements=(SimpleNamespace(chain=chain, observed_target_indices=(0, 2), observed_residues=(residues[0], residues[2])),),
        gaps=(SimpleNamespace(chain=chain, target_interval=SimpleNamespace(start=1, stop=2), generated_residues=(residues[1],)),),
        fixed_atoms=tuple(fixed), generated_atoms=tuple(restored),
        retained_explicit_links=(SimpleNamespace(atom1=AtomIdentity(chain, "1", "", "SG"), atom2=AtomIdentity(chain, "3", "", "SG")),))
    return request, "".join(lines)


class NativeArray(SimpleNamespace):
    def __len__(self):
        return len(self.coord)


def native(request, manifest, translation=(0., 0., 0.)):
    rows = manifest["axis"]
    indices = {(r["atom"]["residue_number"], r["atom"]["atom_name"]): i for i, r in enumerate(rows)}
    pairs = [(indices["1", "C"], indices["2", "N"], 1),
             (indices["2", "C"], indices["3", "N"], 1),
             (indices["1", "SG"], indices["3", "SG"], 1)]
    bonds = np.array(pairs)
    return NativeArray(src_component=np.array([r["atom"]["chain"] + r["atom"]["residue_number"] for r in rows]),
        atom_name=np.array([r["atom"]["atom_name"] for r in rows]), res_name=np.array([r["residue_name"] for r in rows]),
        coord=np.array([np.asarray(r["coordinates"]) + translation if r["fixed"] else np.zeros(3) for r in rows]),
        is_motif_atom_with_fixed_coord=np.array([r["fixed"] for r in rows]),
        is_motif_atom_with_fixed_seq=np.ones(len(rows), dtype=bool),
        chain_id=np.array([r["atom"]["chain"] for r in rows]),
        res_id=np.array([int(r["atom"]["residue_number"]) for r in rows]), ins_code=np.array([""] * len(rows)),
        bonds=SimpleNamespace(as_array=lambda: bonds))


def test_complete_target_mask_and_terminal_OXT_are_separate_from_native_axis():
    request, source = fixture("a")
    original = source
    text, manifest = conditioning(request, source)
    assert source == original
    assert manifest["contig"] == "a1-3"
    assert manifest["select_fixed_atoms"] == {"a1": "ALL", "a2": "", "a3": "ALL"}
    assert len(manifest["axis"]) == 16
    assert len(manifest["preserved_only_atoms"]) == 1
    assert " OXT" not in text
    assert sum(row["fixed"] for row in manifest["axis"]) == 12
    assert all(row["coordinates"] == (0., 0., 0.) for row in manifest["axis"] if not row["fixed"])


def test_input_annotations_are_preserved_and_conect_is_reindexed():
    request, source = fixture()
    ssbond = "SSBOND   1 CYS A    1    CYS A    3\n"
    text, _ = conditioning(request, ssbond + source + "CONECT    6   12\n")
    assert ssbond in text
    assert "CONECT    6   16\n" in text
    with pytest.raises(ValueError, match="dangling"):
        conditioning(request, source + "CONECT    6 9999\n")


def test_common_centering_translation_does_not_mask_source_distortion():
    request, source = fixture()
    _, manifest = conditioning(request, source)
    result = check_array(native(request, manifest, (17., -8., 4.)), manifest, request)
    assert result["graph_transport_passed"]
    assert result["expected_peptide_bonds"] == result["expected_observed_missing_junctions"] == 2
    assert result["expected_retained_links"] == 1
    assert result["translation_angstrom"] == [17., -8., 4.]
    assert result["max_observed_translation_residual_angstrom"] < 1e-12


@pytest.mark.parametrize("change", ["sequence", "fixed-mask", "coordinates", "dummy-leak", "identity", "residue-collapse", "chain-split"])
def test_native_transport_changes_fail_closed(change):
    request, source = fixture()
    _, manifest = conditioning(request, source)
    array = native(request, manifest)
    if change == "sequence":
        array.is_motif_atom_with_fixed_seq[0] = False
    elif change == "fixed-mask":
        array.is_motif_atom_with_fixed_coord[0] = False
    elif change == "coordinates":
        array.coord[0, 0] += 1.
    elif change == "dummy-leak":
        array.coord[~array.is_motif_atom_with_fixed_coord, 0] = 1.
    elif change == "identity":
        array.src_component[-1] = array.src_component[0]
    elif change == "residue-collapse":
        array.res_id[array.res_id == 3] = 1
    else:
        array.chain_id[array.res_id == 2] = "B"
    with pytest.raises(ValueError):
        check_array(array, manifest, request)


@pytest.mark.parametrize("remove", ["peptide", "link"])
def test_missing_native_bond_is_reported_as_failure_not_silently_accepted(remove):
    request, source = fixture()
    _, manifest = conditioning(request, source)
    array = native(request, manifest)
    bonds = array.bonds.as_array()
    bonds = bonds[1:] if remove == "peptide" else bonds[:-1]
    array.bonds = SimpleNamespace(as_array=lambda: bonds)
    result = check_array(array, manifest, request)
    assert not result["graph_transport_passed"]
    assert result["missing_junction_bonds"] == (1 if remove == "peptide" else 0)
    assert result["missing_retained_links"] == (1 if remove == "link" else 0)


def test_extra_non_neighbor_bond_is_not_a_valid_continuity_bridge():
    request, source = fixture()
    _, manifest = conditioning(request, source)
    array = native(request, manifest)
    lookup = {(r["atom"]["residue_number"], r["atom"]["atom_name"]): i for i, r in enumerate(manifest["axis"])}
    bonds = np.concatenate((array.bonds.as_array(), [[lookup["1", "C"], lookup["3", "N"], 1]]))
    array.bonds = SimpleNamespace(as_array=lambda: bonds)
    result = check_array(array, manifest, request)
    assert result["unexpected_inter_residue_bonds"] == 1
    assert not result["graph_transport_passed"]


def test_source_inventory_or_ownership_cannot_be_rewritten_to_fit():
    request, source = fixture()
    changed = deepcopy(request)
    changed.fixed_atoms = request.fixed_atoms[:-1]
    with pytest.raises(ValueError, match="inventory"):
        conditioning(changed, source)


def test_duplicate_serials_cannot_retarget_conditioning_conect():
    request, source = fixture()
    lines = source.splitlines(keepends=True)
    lines[1] = lines[1][:6] + lines[0][6:11] + lines[1][11:]
    with pytest.raises(ValueError, match="repeated source atom serial"):
        conditioning(request, "".join(lines))


def test_missing_native_annotations_are_a_transport_refusal():
    request, source = fixture()
    _, manifest = conditioning(request, source)
    array = native(request, manifest)
    del array.src_component
    with pytest.raises(ValueError, match="lost required transport annotations"):
        check_array(array, manifest, request)


def test_native_module_import_does_not_require_parent_science_executor():
    code = '''
import importlib.abc
import sys
class RejectParentScience(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname in {"gemmi", "openmm", "scripts.run_mac_modeller_comparison"}:
            raise RuntimeError("native data module imported parent-only scientific dependencies")
sys.meta_path.insert(0, RejectParentScience())
import scripts.preflight_retained_rfd3_full_target
'''
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
