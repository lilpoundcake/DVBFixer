import pytest
from scripts.analyze_1non_refinement import callback_summary, classify_clashes


def clash(a, b):
    def atom(a):
        return dict(zip(("chain", "resid", "icode", "atom"), a, strict=True))
    return {"severity": "ERROR", "category": "clash", "details": {
        "atoms": [atom(a), atom(b)], "overlap_angstrom": 1.2}}


def test_clash_ownership_is_case_sensitive_and_counts_contacts_not_residues():
    a, b, c, d = [(chain, "82", "", "CA") for chain in ("A", "a", "B", "D")]
    report = classify_clashes([clash(a, b), clash(a, c), clash(b, d)], {a, c}, {b, d})
    assert report["counts"] == {"generated-fixed": 1, "generated-generated-cross-chain": 1, "fixed-fixed": 1}
    assert sum(r["contacts"] for r in report["residue_contact_counts"]) == 6
    with pytest.raises(ValueError, match="ownership"):
        classify_clashes([clash(a, b)], {a}, set())


def test_callback_summary_does_not_infer_convergence_from_energy_drop():
    rows = [{"stationarity_passed": False, "projected_cartesian_dof_rms_kj_mol_nm": 25.,
             "energy_kj_mol": -i} for i in range(110)]
    summary = callback_summary(rows)
    assert summary["stationarity_pass_callbacks"] == 0
    assert summary["last_100"]["energy_change_kj_mol"] == -99
    assert summary["last_100"]["rms_min"] == 25
