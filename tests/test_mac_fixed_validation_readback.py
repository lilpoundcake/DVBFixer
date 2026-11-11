import pytest
from scripts.audit_mac_fixed_validation import check_keys, check_records


def test_reserved_cohort_rejects_duplicate_missing_and_undeclared_outcomes():
    expected = {("case", "model", 211), ("case", "model", 223)}
    rows = [{"case": "case", "backend": "model", "seed": seed} for seed in (211, 223)]
    check_keys(rows, expected)
    for invalid in (rows[:1], rows + rows[:1], [rows[0], {**rows[1], "seed": 7}]):
        with pytest.raises(ValueError, match="outcomes"):
            check_keys(invalid, expected)


def test_terminal_oxt_record_is_exact_observation_not_just_same_coordinates():
    line = "ATOM      1  OXT ALA a  82A      1.000   2.000   3.000  1.00 20.00           O\n"
    assert check_records(line, line) == 1
    with pytest.raises(ValueError, match="observed"):
        check_records(line, line.replace("20.00", "21.00"))
    with pytest.raises(ValueError, match="duplicate"):
        check_records(line, line + line)
    with pytest.raises(ValueError, match="observed"):
        check_records(line, "")
