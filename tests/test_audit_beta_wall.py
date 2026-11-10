"""Actual-frame geometry views preserve movement, not original observations."""
from scripts.audit_beta_wall import refs


def test_context_identity_keeps_chain_case_and_insertion_codes():
    rows=[dict(key=['D','82','A','THR','CB']),dict(key=['d','82','','THR','CB'])]
    first,second=refs(rows)
    assert first.component.chain_id=='D' and second.component.chain_id=='d'
    assert first.component.insertion_code=='A' and second.component.insertion_code==''
    assert first!=second
