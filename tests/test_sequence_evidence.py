"""Conservative observed-backbone hints for shared sequence placement."""

from dvbfixer.pdbutils.sequence_evidence import observed_peptide_continuity


def atom(serial, chain, number, code, name, x, alt=" "):
    return f"ATOM  {serial:5d} {name:>4}{alt}ALA {chain}{number:4d}{code or ' '}   {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C  \n"


def test_exact_case_sensitive_chain_and_insertion_code_geometry():
    lines = [atom(1, "D", 82, "", "C", 0), atom(2, "D", 82, "A", "N", 1.33),
        atom(3, "D", 82, "A", "C", 3), atom(4, "D", 85, "", "N", 15),
        atom(5, "d", 82, "A", "N", 100)]
    assert observed_peptide_continuity(lines, "D", ((82, ""), (82, "A"), (85, ""))) == (True, False)
    assert observed_peptide_continuity(lines, "d", ((82, ""), (82, "A"))) == (None,)


def test_ter_duplicate_altloc_and_missing_atoms_are_not_false_peptide_evidence():
    left, right = atom(1, "A", 1, "", "C", 0), atom(2, "A", 2, "", "N", 1.33)
    keys = ((1, ""), (2, ""))
    assert observed_peptide_continuity([left, "TER\n", right], "A", keys) == (False,)
    assert observed_peptide_continuity([left, right, right], "A", keys) == (None,)
    assert observed_peptide_continuity([left, atom(2, "A", 2, "", "N", 1.33, "A")], "A", keys) == (None,)
    assert observed_peptide_continuity([left], "A", keys) == (None,)
