from pathlib import Path

import pytest
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_conditioning_contracts import materialize, prepare
from scripts.prepare_mac_comparison_expansion import CASES
from scripts.prepare_mac_comparison_expansion import prepare as prepare_case

from dvbfixer.model.diffusion.contract import AtomIdentity


@pytest.fixture
def tiny(tmp_path):
    from scripts.prepare_reference_reconstruction_case import prepare as prepare_tiny
    case = tmp_path / "tiny"
    prepare_tiny(Path("tests/fixtures/protein_context/1UBQ.pdb"), case, chain="A", first=1, last=72,
                 hidden_first=10, hidden_last=13, case_id="test-control")
    return case


def test_context_is_source_only_and_keeps_the_same_generated_mask(tiny, tmp_path):
    receipts = [prepare(tiny, tmp_path / context, gaps=(0,), context=context, information="all-heavy")
                for context in ("flanks8", "spatial6", "full")]
    assert all(r["target_generated_atoms"] == 28 for r in receipts)
    counts = [len(r["axis"]) for r in receipts]
    assert counts[0] <= counts[1] <= counts[2] == 574
    assert receipts[0]["source_target_spans"] == [[1, 21]]
    assert all(not r["reference_coordinates_read"] for r in receipts)


def test_backbone_ablation_withholds_sidechains_without_changing_output_ownership(tiny, tmp_path):
    output = tmp_path / "bb"
    receipt = prepare(tiny, output, gaps=(0,), context="full", information="backbone")
    assert receipt["sampler_fixed_atoms"] == 68*4
    assert receipt["sampler_unfixed_atoms"] == 574-68*4
    assert receipt["target_generated_atoms"] == 28
    rows = {AtomIdentity(**r["atom"]): r for r in receipt["axis"]}
    for atom, value in inventory(output / "conditioning.pdb").items():
        row = rows[atom]
        if not row["sampler_fixed"]:
            assert (value[2] == 0).all()
        if row["original_observation"]:
            assert not row["target_generated"]


def test_independent_projection_omits_other_gap_not_observations(tmp_path):
    case = tmp_path / "pair"
    prepare_case(CASES[1], case)
    receipt = prepare(case, tmp_path / "input", gaps=(0,), context="full", information="all-heavy")
    represented = {(r["atom"]["chain"], r["atom"]["residue_number"]) for r in receipt["axis"]}
    assert not represented & {("A", str(i)) for i in range(46, 50)}
    assert receipt["sampler_fixed_atoms"] == receipt["full_source_fixed_atoms"]
    assert receipt["source_target_spans"] == [[0, 45], [49, 72]]


def test_materializer_restores_all_observed_atom_records_and_rejects_overlap(tiny):
    native = inventory(tiny / "reference-projection-NOT-CONDITIONING.pdb")
    source = inventory(tiny / "source.pdb")
    generated = {a: v[2] for a, v in native.items() if a not in source}
    text = materialize(tiny, [generated])
    original_lines = [line for line in (tiny / "source.pdb").read_text().splitlines() if line.startswith("ATOM  ")]
    assert all(line in text.splitlines() for line in original_lines)
    with pytest.raises(ValueError, match="overlapping"):
        materialize(tiny, [generated, generated])
    with pytest.raises(ValueError, match="cover all and only"):
        materialize(tiny, [])
