import json

import pytest
from scripts.audit_mac_expansion_results import audit as readback
from scripts.audit_mac_reconstruction_smoke import audit
from scripts.mac_reconstruction_comparison import digest, save
from scripts.prepare_mac_comparison_expansion import CASES, prepare
from scripts.prepare_mac_expansion_inputs import prepare as encode


@pytest.fixture
def study(tmp_path):
    preparation, case = tmp_path / "preparation", tmp_path / "case"
    prepare(CASES[0], preparation)
    encode(preparation, case)
    root = tmp_path / "study"
    work = root / "runs/synthetic-control"
    work.mkdir(parents=True)
    audit(case, case / "reference-projection-NOT-CONDITIONING.pdb", work / "audit", seed=7)
    final = work / "audit/source-restored.pdb"
    save(root / "runs/plan.json", {"jobs": [{"case": str(case), "work": str(work), "arm": "synthetic-test-control", "seed": 7}]})
    save(root / "final-results.json", [{"case": case.name, "backend": "synthetic-test-control", "seed": 7,
        "final_path": str(final), "sha256": digest(final), "heldout_heavy_rmsd": 0., "repaired": False}])
    return root, final


def test_independent_readback_scores_and_checks_observations(study, tmp_path):
    root, _ = study
    output = tmp_path / "readback"
    readback(root, output)
    row = json.loads((output / "independent-results.json").read_text())["rows"][0]
    assert row["observed_atom_records_byte_identical"] and row["global_all_l"]
    assert row["final_heldout_heavy_rmsd"] < 1e-10
    assert row["per_gap_heavy_rmsd"] == [0.]


def test_independent_readback_rejects_fixed_atom_mutation(study, tmp_path):
    root, final = study
    lines = final.read_text().splitlines(keepends=True)
    index = next(i for i, line in enumerate(lines) if line.startswith("ATOM  ") and int(line[22:26]) == 5)
    line = lines[index]
    lines[index] = line[:30]+f"{float(line[30:38])+0.1:8.3f}"+line[38:]
    final.write_text("".join(lines))
    rows = json.loads((root / "final-results.json").read_text())
    rows[0]["sha256"] = digest(final)
    (root / "final-results.json").write_text(json.dumps(rows))
    with pytest.raises(ValueError, match="observed ATOM bytes"):
        readback(root, tmp_path / "refused")


def test_modeller_summary_uses_actual_cohort_size_and_recovery_metadata(study, tmp_path):
    from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS
    from scripts.summarize_mac_modeller_comparison import summarize

    root, final = study
    original_plan = json.loads((root / "runs/plan.json").read_text())
    case = original_plan["jobs"][0]["case"]
    from pathlib import Path
    case = Path(case)
    plan = {"jobs": [{"seed": seed} for seed in DEVELOPMENT_SEEDS], "technical_amendment": None,
        "source_code_and_input_sha256": {str(case / name): digest(case / name) for name in
            ("source.pdb", "evaluation-mask.json", "reference-projection-NOT-CONDITIONING.pdb")}}
    save(root / "plan.json", plan)
    save(root / "completed.json", {"completed": 6})
    rows = [{"case": str(case), "case_id": case.name, "seed": seed, "completed": True, "returncode": 0,
        "elapsed_seconds": 1., "peak_rss_bytes": 1000, "result": {"seed": seed, "final_path": str(final), "final_sha256": digest(final),
            "final_source_geometry_and_l_passed": True, "raw_restored_geometry_passed": True, "repaired": False,
            "final_heldout_heavy_rmsd": 0., "final_heldout_backbone_rmsd": 0., "native_observed_drift": {"max": 0.}}}
        for seed in DEVELOPMENT_SEEDS]
    save(root / "progress.json", rows)
    output = tmp_path / "modeller-summary"
    summarize(root, output)
    summary = json.loads((output / "summary.json").read_text())
    assert summary["declared_outcomes"] == 6
    assert summary["all_declared_outcomes_included_in_rmsd_summaries"]
    assert "recovered" not in (output / "README.md").read_text()
