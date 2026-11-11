#!/usr/bin/env python3
"""A separately declared barnase/barstar interface admission control, no sampling.

Does not replace/remove the refused 1A3N cell. No search over alternate masks or
source repairs. Full explicit A3-110 / D3-89 protein projection, hidden A83-90.
"""
from __future__ import annotations

import argparse
import traceback
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save
from scripts.prepare_mac_comparison_expansion import prepare

SOURCE = "tmp/real-complex-cc94-2026-11-06/sources/1BRS.pdb"
EXPECTED = "0076aae54252706a8c8c02f7b1e71ce2394b0e928b413590b4ac12aade9e8f6c"
SPEC = {"id": "1BRS-A3-110-D3-89-gap83-90", "stratum": "barnase/barstar interface candidate",
        "fixture": SOURCE, "spans": [["A", 3, 110], ["D", 3, 89]], "hidden": [["A", 83, 90]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parents[1] / SOURCE
    if digest(source) != EXPECTED:
        raise ValueError("retained experimental 1BRS source changed")
    save(root / "protocol.json", {"spec": SPEC, "source_sha256": EXPECTED, "separate_from_refused_1A3N": True,
        "source_provenance": "retained RCSB experimental download from declared real-complex cohort",
        "mask_selection": "fixed before preparation; no posthoc mask replacement", "sampling_calls": 0,
        "script_sha256": digest(Path(__file__)), "qualification": False})
    (root / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    try:
        manifest = prepare(SPEC, root / SPEC["id"])
        coords = inventory(root / SPEC["id"] / "source.pdb")
        # Only still-observed anchors/context enter this context diagnostic.
        anchors = [v[2] for a, v in coords.items() if a.chain == "A" and a.residue_number in {"82", "91"}]
        partner = [v[2] for a, v in coords.items() if a.chain == "D"]
        distances = np.linalg.norm(np.array(anchors)[:, None, :]-np.array(partner)[None, :, :], axis=-1)
        save(root / "result.json", {"mask_admitted": True, "manifest": manifest,
            "observed_anchor_partner_min_distance_angstrom": float(distances.min()),
            "observed_anchor_partner_atom_pairs_within_6A": int((distances <= 6).sum()),
            "interface_claim": "anchor-context proximity only; hidden coordinates not used", "sampling_executed": False})
    except Exception:
        save(root / "result.json", {"mask_admitted": False, "error": traceback.format_exc(), "sampling_executed": False})
    print((root / "result.json").read_text())


if __name__ == "__main__":
    main()
