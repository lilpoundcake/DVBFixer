"""All pairs retain their exact matched denominators, including empty pairs."""
import pytest
from scripts.summarize_prospective_guard_pairs import METRICS, summarize


def value(number):
    return dict(rounded=dict.fromkeys(METRICS,number))


def test_all_pairs_not_only_favorable_pair_are_reported():
    cases=dict(x=dict(accuracies=dict(a=value(1.),b=value(2.))),y=dict(accuracies=dict(a=value(3.),c=value(4.))))
    rows=summarize(['a','b','c'],cases)
    assert [(r['first'],r['second'],r['matched_count']) for r in rows]==[('a','b',1),('a','c',1),('b','c',0)]
    assert rows[-1]['medians']['b'][METRICS[0]] is None


def test_nonfinite_saved_accuracy_is_not_an_implicit_exclusion():
    with pytest.raises(ValueError,match='finite'):
        summarize(['a','b'],dict(x=dict(accuracies=dict(a=value(float('nan')),b=value(1.)))))
