import pytest
from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS
from scripts.summarize_mac_expansion_panel import aggregate


def controls():
    return [{"case": "synthetic-mask", "backend": "synthetic-control", "seed": seed,
        "raw_restored_geometry_passed": index < 2, "final_geometry_passed": index < 5,
        "final_admission_failures": [], "global_all_l": True, "observed_atom_records_byte_identical": True,
        "repaired": index >= 2, "final_heldout_heavy_rmsd": [1., 2., 3., 4., 5., 100.][index],
        "final_heldout_backbone_rmsd": float(index), "final_sha256": str(seed),
        "native_observed_drift": {"max": 2.}}
        for index, seed in enumerate(DEVELOPMENT_SEEDS)]


def test_summary_preserves_failures_and_all_outcome_accuracy():
    group = aggregate(controls())[0]
    assert group["raw_restored_geometry_passes"] == 2
    assert group["final_source_geometry_and_l_passes"] == 5
    assert group["heavy_rmsd"]["median"] == 3.5
    assert group["heavy_rmsd"]["max"] == 100.
    assert group["distinct_final_outputs"] == 6
    assert group["native_observed_max_drift"] == 2.


@pytest.mark.parametrize("duplicate", [False, True])
def test_summary_refuses_incomplete_or_repeated_seeds(duplicate):
    rows = controls()
    if duplicate:
        rows[-1]["seed"] = rows[0]["seed"]
    else:
        rows.pop()
    with pytest.raises(ValueError, match="every development seed exactly once"):
        aggregate(rows)
