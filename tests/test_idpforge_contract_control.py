"""Lightweight inspection of private probe; no Torch/upstream/checkpoint import."""
import ast

import pytest
from scripts.idpforge_contract_control import extract_recon


def test_extracts_only_recon_not_module_imports_or_other_methods():
    source = '''
import unavailable_native_module
raise RuntimeError("must not execute module")
class IDPForge:
    def forward(self):
        raise RuntimeError("must not call real forward")
    @torch.no_grad()
    def recon(self, value):
        return value
'''
    method = extract_recon(source)
    assert method.name == "recon"
    assert len(method.body) == 1
    assert isinstance(method.body[0], ast.Return)
    assert ast.unparse(method.decorator_list[0]) == "torch.no_grad()"


@pytest.mark.parametrize("source", [
    "class Other: pass",
    "class IDPForge: pass",
    "class IDPForge:\n def recon(self): pass\n def recon(self): pass",
    "class IDPForge:\n def recon(self): pass\nclass IDPForge: pass",
])
def test_refuses_missing_or_ambiguous_upstream_method(source):
    with pytest.raises(ValueError, match="exactly one"):
        extract_recon(source)
