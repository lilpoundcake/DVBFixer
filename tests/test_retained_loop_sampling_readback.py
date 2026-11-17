import pytest
from scripts.audit_retained_loop_context_sampling import fixed_records

from dvbfixer.model.diffusion.contract import AtomIdentity


def atom(chain="a", icode="A"):
    return f"ATOM      1  CA  ALA {chain}  82{icode}      1.000   2.000   3.000  1.00 20.00           C\n"


def test_raw_endpoint_custody_keeps_case_and_insertion_code():
    fixed = {AtomIdentity("a", "82", "A", "CA")}
    line = atom()
    assert list(fixed_records(line, fixed).values()) == [line]
    with pytest.raises(ValueError, match="omitted"):
        fixed_records(atom("A"), fixed)
    with pytest.raises(ValueError, match="duplicate"):
        fixed_records(line+line, fixed)
