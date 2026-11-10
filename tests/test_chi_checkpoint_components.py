"""Native component text is complete, finite and separate from cutoff certainty."""
import pytest
from scripts.audit_chi_checkpoint_components import component_cutoff, parse_components


def line(i,n,analytic):
    return f'core.optimization: deriv_dev: {i} {n} 1.0 1 6 4 {analytic} 1E-9 0.0\n'


def test_complete_component_vectors_parse_without_sparse_threshold_filtering():
    rows=parse_components(line(1,2,'8E-5')+line(2,2,'6E-5'),2)
    assert [r['analytic_per_degree'] for r in rows[0]]==[8e-5,6e-5]
    assert component_cutoff([8e-5,6e-5],1e-4)['certainly_at_or_below_gmax']


def test_incomplete_or_nonfinite_native_vector_refuses():
    for text in (line(1,2,'0.1'),line(1,2,'nan')+line(2,2,'0.1')):
        with pytest.raises(ValueError):
            parse_components(text,2)


def test_rounding_near_cutoff_is_not_misclassified_as_certainty():
    row=component_cutoff([1e-4],1e-4)
    assert not row['certainly_at_or_below_gmax'] and not row['certainly_above_gmax']
    assert component_cutoff([2e-4],1e-4)['certainly_above_gmax']
