"""Read-only final-data guard for the private no-structural-template E arm."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from mps_research_smoke import run
from runner.inference import InferenceRunner

ORIGINAL = InferenceRunner.predict


def guarded(self, data, *args, **kwargs):
    masks = {}
    for name in ("template_pseudo_beta_mask", "template_backbone_frame_mask"):
        values = data["input_feature_dict"][name].detach().cpu().numpy()
        count = int(np.count_nonzero(values))
        if count:
            raise ValueError("no-template E arm received a populated structural template mask")
        masks[name] = {"shape": list(values.shape), "nonzero": count}
    Path(os.environ["DVBFIXER_SEQUENCE_EVIDENCE_RECEIPT"]).write_text(json.dumps({"masks": masks,
        "no_structural_template_features": True, "read_only_data_inspection": True,
        "model_numerical_method_unchanged": True}, indent=2)+"\n")
    return ORIGINAL(self, data, *args, **kwargs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "input", "assets", "kalign", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--seed", type=int, choices=(7, 42, 101, 102, 103, 104), required=True)
    args = parser.parse_args()
    InferenceRunner.predict = guarded
    run(args)
