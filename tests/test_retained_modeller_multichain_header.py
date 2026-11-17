from types import SimpleNamespace

from scripts import run_retained_loop_modeller_controls as controls


def test_private_multichain_pir_uses_case_sensitive_explicit_endpoint_chains(monkeypatch):
    original = ">P1;source\nstructureX:source:FIRST:@:LAST:@::::\nAA/AA*\n>P1;target\nsequence:target::::::::\nAA/AA*\n"
    request = SimpleNamespace(target_sequences=(SimpleNamespace(chain="a"), SimpleNamespace(chain="A")))
    monkeypatch.setattr(controls, "modeller_alignment", lambda request: original)
    changed = controls.multichain_alignment(request)
    assert "structureX:source:FIRST:a:LAST:A::::" in changed
    assert changed.replace("FIRST:a:LAST:A", "FIRST:@:LAST:@") == original
