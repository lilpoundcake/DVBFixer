"""Private approved-topology graph controls; invented identities, no coordinates."""
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError
from dvbfixer.model.diffusion.rfd3_native_topology import validate_native_graph


def data():
    a, b, c, d, e, f = [AtomIdentity(chain, number, code, "SG") for chain, number, code in
                         (("A", "77", ""), ("C", "77", ""), ("A", "7", ""),
                          ("A", "16", ""), ("C", "7", "A"), ("C", "16", ""))]
    return dict(cysteines=frozenset((a, b, c, d, e, f)), observed=frozenset((a, b, d, f)),
                active=frozenset((c, e)), explicit_source_pairs=(), source_inferred_pairs=((a, b),),
                target_pairs=((a, b), (c, d), (e, f)))


def test_native_target_keeps_inferred_observed_link_and_missing_atom_roles():
    result = validate_native_graph(**data())
    assert result.target_pairs == 3 and result.cysteine_endpoints == 6
    assert result.explicit_source_pairs == 0 and result.source_inferred_observed_pairs == 1
    assert len(result.additional_active_pairs) == 2
    assert not result.source_annotations_rewritten and not result.scientific_acceptance


@pytest.mark.parametrize("mutation", ["overlap", "duplicate", "lost-source", "two-partners", "self", "unknown", "inactive", "inferred-generated"])
def test_unapproved_or_inconsistent_graphs_fail_closed(mutation):
    values = data()
    a, b = values["source_inferred_pairs"][0]
    c, d = values["target_pairs"][1]
    e, f = values["target_pairs"][2]
    if mutation == "overlap":
        values["active"] |= {a}
    elif mutation == "duplicate":
        values["target_pairs"] += ((b, a),)
    elif mutation == "lost-source":
        values["target_pairs"] = ((a, d), (c, b), (e, f))
    elif mutation == "two-partners":
        values["target_pairs"] += ((a, c),)
    elif mutation == "self":
        values["target_pairs"] += ((a, a),)
    elif mutation == "unknown":
        values["target_pairs"] += ((a, AtomIdentity("a", "77", "", "SG")),)
    elif mutation == "inactive":
        values["active"] -= {e}
    else:
        values["source_inferred_pairs"] += ((c, d),)
    with pytest.raises(CustodyError):
        validate_native_graph(**values)


def test_additional_observed_pair_requires_separate_source_only_inference():
    values = data()
    a, b = values["source_inferred_pairs"][0]
    values["source_inferred_pairs"] = ()
    with pytest.raises(CustodyError, match="between observations"):
        validate_native_graph(**values)
    values["source_inferred_pairs"] = ((a, b),)
    assert validate_native_graph(**values).source_inferred_observed_pairs == 1
