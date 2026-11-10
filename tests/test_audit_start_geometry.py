"""Exact actual-frame identities must not collapse chain case or insertion codes."""
import pytest
from scripts.audit_start_geometry import heavy_map


def test_distinct_chain_case_and_insertion_codes_remain_distinct():
    rows=[dict(key=['D','82','A','THR','CB'],element='C',xyz=(0.,0.,0.)),
          dict(key=['d','82','','THR','CB'],element='C',xyz=(1.,0.,0.))]
    assert len(heavy_map(rows))==2


def test_duplicate_heavy_atoms_refuse_instead_of_overwriting():
    atom=dict(key=['A','2','','THR','CB'],element='C',xyz=(0.,0.,0.))
    with pytest.raises(ValueError,match='duplicate'):
        heavy_map([atom,atom])
