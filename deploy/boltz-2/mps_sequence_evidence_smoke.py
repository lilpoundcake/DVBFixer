"""Read-only final-consumer template-mask guard around the unchanged stock CLI."""
from __future__ import annotations

import json
import os
import runpy
from pathlib import Path

import numpy as np
from boltz.model.models.boltz2 import Boltz2

ORIGINAL = Boltz2.forward


def guarded(self, feats, *args, **kwargs):
    masks = {}
    for name in ("template_mask", "template_mask_cb", "template_mask_frame"):
        values = feats[name].detach().cpu().numpy()
        count = int(np.count_nonzero(values))
        if count:
            raise ValueError("no-template E arm received a populated structural template mask")
        masks[name] = {"shape": list(values.shape), "nonzero": count}
    Path(os.environ["DVBFIXER_SEQUENCE_EVIDENCE_RECEIPT"]).write_text(json.dumps({"masks": masks,
        "no_structural_template_features": True, "read_only_data_inspection": True,
        "model_numerical_method_unchanged": True}, indent=2)+"\n")
    return ORIGINAL(self, feats, *args, **kwargs)


if __name__ == "__main__":
    Boltz2.forward = guarded
    runpy.run_module("boltz.main", run_name="__main__")
