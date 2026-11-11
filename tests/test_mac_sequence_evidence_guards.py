"""Exercise read-only guard bodies without importing or launching ML runtimes."""
import ast
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

PROJECT = Path(__file__).resolve().parents[1]


class Tensor:
    def __init__(self, values):
        self.values = np.asarray(values)

    def detach(self):
        return self

    def cpu(self):
        return self

    def numpy(self):
        return self.values


def load_guard(backend, original):
    source = PROJECT / "deploy" / backend / "mps_sequence_evidence_smoke.py"
    tree = ast.parse(source.read_text())
    body = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "guarded"]
    assert len(body) == 1
    namespace = {"np": np, "Path": Path, "json": json, "os": os, "ORIGINAL": original}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(source), "exec"), namespace)
    return namespace["guarded"]


@pytest.mark.parametrize("backend", ["boltz-2", "protenix-v1"])
def test_no_template_guard_forwards_original_objects_and_rejects_populated_masks(backend, tmp_path, monkeypatch):
    receipt = tmp_path / "receipt.json"
    monkeypatch.setenv("DVBFIXER_SEQUENCE_EVIDENCE_RECEIPT", str(receipt))
    calls = []

    def original(*args, **kwargs):
        calls.append((args, kwargs))
        return "unchanged model output"

    guard = load_guard(backend, original)
    names = (["template_mask", "template_mask_cb", "template_mask_frame"] if backend == "boltz-2"
             else ["template_pseudo_beta_mask", "template_backbone_frame_mask"])
    features = {name: Tensor(np.zeros((1, 4))) for name in names}
    data = features if backend == "boltz-2" else {"input_feature_dict": features}
    model = SimpleNamespace()
    assert guard(model, data, 17, preserved=True) == "unchanged model output"
    assert calls == [((model, data, 17), {"preserved": True})]
    report = json.loads(receipt.read_text())
    assert report["no_structural_template_features"]
    assert set(report["masks"]) == set(names)
    assert all(row["nonzero"] == 0 for row in report["masks"].values())
    before = receipt.read_bytes()
    for name in names:
        features[name].values[0, 0] = 1
        with pytest.raises(ValueError, match="populated structural template"):
            guard(model, data)
        features[name].values[0, 0] = 0
    assert len(calls) == 1
    assert receipt.read_bytes() == before


@pytest.mark.parametrize("backend", ["boltz-2", "protenix-v1"])
def test_no_template_guard_fails_closed_on_missing_consumed_feature(backend):
    guard = load_guard(backend, lambda *args, **kwargs: pytest.fail("model must not execute"))
    data = {} if backend == "boltz-2" else {"input_feature_dict": {}}
    with pytest.raises(KeyError):
        guard(None, data)
