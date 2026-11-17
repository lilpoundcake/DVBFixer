import pytest
from scripts.audit_retained_sampler_return_repairs import check_observed_records

from dvbfixer.model.diffusion.contract import AtomIdentity


def test_repair_readback_allows_only_declared_coordinate_fields():
    atom = AtomIdentity("a", "82", "A", "CA")
    source = "ATOM      1  CA  ALA a  82A      1.000   2.000   3.000  1.00 20.00           C\n"
    moved = source[:30]+"   1.100"+source[38:]
    check_observed_records(source, moved, {atom}, {atom})
    with pytest.raises(ValueError, match="nonmobile"):
        check_observed_records(source, moved, {atom}, set())
    renamed = moved[:17]+"GLY"+moved[20:]
    with pytest.raises(ValueError, match="metadata"):
        check_observed_records(source, renamed, {atom}, {atom})
