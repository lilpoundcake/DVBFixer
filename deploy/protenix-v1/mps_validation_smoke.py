"""Hash-bound reserved-seed entry to the unchanged private MPS inference method."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from mac_validation_policy import SEEDS, require_policy
from mps_research_smoke import run

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "input", "assets", "kalign", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True, choices=SEEDS)
    args = parser.parse_args()
    permission = require_policy("protenix", args.seed)
    (args.output.parent / "validation-policy.json").write_text(json.dumps(permission, indent=2)+"\n")
    run(args, allowed_seeds=SEEDS)
