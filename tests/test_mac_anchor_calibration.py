from pathlib import Path

import pytest
from scripts.prepare_reference_reconstruction_case import prepare
from scripts.run_mac_anchor_calibration import affected_anchors

from dvbfixer.model.diffusion.contract import DiffusionRequest


@pytest.fixture
def mask_request(tmp_path):
    case = tmp_path / "case"
    prepare(Path("tests/fixtures/protein_context/1UBQ.pdb"), case, chain="A", first=1, last=72,
            hidden_first=10, hidden_last=13, case_id="control")
    return DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())


def finding(category, *numbers, actionable=True):
    return {"issues": [{"category": category, "actionable": actionable,
        "residues": [{"chain": "A", "residue_number": str(n), "insertion_code": ""} for n in numbers]}]}


def test_only_residual_junction_findings_implicating_observed_flanks_allow_release(mask_request):
    assert [r.residue_number for r in affected_anchors(mask_request, finding("non_planar_amide", 9, 10))] == ["9"]
    assert [r.residue_number for r in affected_anchors(mask_request, finding("junction-peptide-connectivity", 13, 14))] == ["14"]
    assert not affected_anchors(mask_request, finding("non_planar_amide", 11, 12))
    assert not affected_anchors(mask_request, finding("clash", 9, 11))
    assert not affected_anchors(mask_request, finding("non_planar_amide", 9, 10, actionable=False))


@pytest.mark.parametrize("pdb,chain,first,last", [("1PGA", "A", 1, 56), ("2CI2", "I", 19, 83)])
def test_observed_terminal_oxt_is_kept_not_deleted_or_mislabeled_generated(tmp_path, pdb, chain, first, last):
    from scripts.prepare_mac_comparison_expansion import prepare as prepare_oxt
    spec = {"id": "terminal-oxygen-control", "fixture": f"tests/fixtures/protein_context/{pdb}.pdb",
        "spans": [[chain, first, last]], "hidden": [[chain, first+9, first+12]]}
    output = tmp_path / "case"
    result = prepare_oxt(spec, output, allow_observed_terminal_oxt=True)
    request = DiffusionRequest.from_json((output / "evaluation-mask.json").read_text())
    assert [a.residue_number for a in request.fixed_atoms if a.atom_name == "OXT"] == [str(last)]
    assert not any(a.atom_name == "OXT" for a in request.generated_atoms)
    assert result["observed_terminal_oxt_preserved"]
    original_oxt = next(line for line in Path(spec["fixture"]).read_text().splitlines() if line.startswith("ATOM  ") and line[12:16].strip() == "OXT")
    assert original_oxt in (output / "source.pdb").read_text().splitlines()
