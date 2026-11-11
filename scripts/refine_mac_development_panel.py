#!/usr/bin/env python3
"""One declared common generated-only repair arm for retained development outputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import traceback
from pathlib import Path

from scripts.audit_mac_reconstruction_smoke import audit
from scripts.refine_mac_reconstruction_smoke import run as refine


def run(root):
    panel = root / "development-panel"
    rows = json.loads((panel / "results.json").read_text())
    output = root / "development-generated-only-repair"
    output.mkdir(exist_ok=False)
    (output / "declaration.json").write_text(json.dumps({
        "scope": "post-result standard repair arm; historical raw passes/failures unchanged",
        "raw_results_sha256": hashlib.sha256((panel / "results.json").read_bytes()).hexdigest(),
        "all_generated_atoms_mobile": True, "observed_atoms_mobile": False,
        "clean_cases_not_rebuilt": True, "platform": "Reference", "max_iterations": 1000,
        "restart_count": 0, "refinement_seed_all_cases": 7,
        "hidden_reference_supplied_to_refinement": False,
        "repair_script_sha256": hashlib.sha256(Path("scripts/refine_mac_reconstruction_smoke.py").read_bytes()).hexdigest(),
        "reserved_seeds_not_used": True, "qualification": False}, indent=2) + "\n")
    results = []
    for row in rows:
        arm, seed = row["arm"], row["seed"]
        work = panel / f"{arm}-seed{seed}"
        dest = output / work.name
        result = {"arm": arm, "sampling_seed": seed, "refinement_seed": 7}
        try:
            refine(work / "audit", dest)
            result["status"] = json.loads((dest / "status.json").read_text())
            if (dest / "refined.pdb").exists():
                result["independent_audit"] = audit(root / "smoke-input-03", dest / "refined.pdb", dest / "audit", seed=seed)
        except Exception:
            result["failure"] = traceback.format_exc()
        results.append(result)
        (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    (output / "completed.json").write_text(json.dumps({"cases_processed": len(results), "qualified": False}) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    run(parser.parse_args().root.absolute())
