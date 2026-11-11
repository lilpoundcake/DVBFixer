import pytest
from scripts.audit_mac_sequence_evidence import check_masks


def test_empty_consumer_masks_are_required_not_just_declared():
    valid = {"no_structural_template_features": True, "masks": {"template_mask": {"nonzero": 0}}}
    check_masks(valid)
    for invalid in ({**valid, "no_structural_template_features": False},
                    {**valid, "masks": {}},
                    {**valid, "masks": {"template_mask": {"nonzero": 1}}}):
        with pytest.raises(ValueError, match="template"):
            check_masks(invalid)
