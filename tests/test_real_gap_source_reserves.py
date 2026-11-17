from scripts.acquire_real_gap_sources import select_reserves
from scripts.audit_1he8_joint_relaxation import read_restraints


def row(name, length, group, control=False):
    return {"entry_id": name, "longest_gap_metadata": length, "natural_gap_count_metadata": 1,
            "conservative_entry_component": group, "small_control_not_primary": control}


def test_long_reserves_are_not_dropped_and_small_control_not_counted_as_primary():
    rows = [row("1HE8", 92, 1), row("LONG", 1066, 2), row("1BRS", 2, 3, control=True),
            row("SHORT", 4, 4), row("MEDIUM", 30, 5)]
    first, strata = select_reserves(rows)
    second, _ = select_reserves(list(reversed(rows)))
    assert [r["entry_id"] for r in first] == [r["entry_id"] for r in second]
    assert {r["entry_id"] for r in first} == {"1HE8", "LONG", "SHORT", "MEDIUM"}
    assert strata[-1]["selected"] == 1
    assert strata[-1]["longest_span"] == [100, None]


def test_historical_restraint_frame_and_case_are_not_reanchored_to_candidate():
    declaration = {"restraints": [{"atom": {"chain": "a", "residue_number": "82", "insertion_code": "A", "atom_name": "CA"},
        "reference_angstrom": [1.0, 2.0, 3.0], "k_kj_mol_nm2": 10000.0}]}
    restraint, = read_restraints(declaration)
    assert restraint.reference_angstrom == (1.0, 2.0, 3.0)
    assert (restraint.atom.chain, restraint.atom.residue_number, restraint.atom.insertion_code) == ("a", "82", "A")
