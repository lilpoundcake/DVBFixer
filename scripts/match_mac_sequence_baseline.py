#!/usr/bin/env python3
"""CPU-match retained template-arm proposals without resampling or overwriting."""
from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save
from scripts.refine_mac_reconstruction_smoke import run as refine
from scripts.run_mac_modeller_comparison import execute

PROJECT = Path(__file__).resolve().parents[1]
OLD = PROJECT / "tmp/mac-reconstruction-2026-11-11"


def declare(output):
    output.mkdir(parents=True, exist_ok=False)
    jobs, files = [], [Path(__file__), Path(refine.__code__.co_filename), Path(execute.__code__.co_filename)]
    for backend in ("boltz", "protenix"):
        for seed in DEVELOPMENT_SEEDS:
            source = OLD / f"{backend}-audit" if seed == 7 else OLD / "development-panel" / f"{backend}-seed{seed}/audit"
            if not (source / "source-restored.pdb").is_file():
                raise ValueError("retained source-template proposal missing")
            job_id = f"{backend}-source-template"
            work = output / "execution" / f"{job_id}-seed-{seed}"
            jobs.append({"case_id": job_id, "case": str(source), "seed": seed, "backend": backend,
                "command": [sys.executable, "-m", "scripts.match_mac_sequence_baseline", "worker", "--output", str(work), "--audit", str(source)]})
            files += [source / name for name in ("source-restored.pdb", "source.pdb", "evaluation-mask.json", "audit.json")]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)), "jobs": jobs,
        "count": 12, "source_code_and_input_sha256": {str(p): digest(p) for p in files},
        "scope": "necessary-only common CPU repair of existing template-arm proposals; old Reference outputs unchanged",
        "new_sampling": False, "reference_based_selection": False, "qualification": False})
    print("Declared 12 retained baseline readbacks/necessary CPU repairs; zero new ML calls.")


def worker(source, output):
    # The guarded executor owns the exclusive reservation and opens these two
    # diagnostic files before starting us. Never adopt an existing payload.
    if not output.is_dir() or {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("baseline worker requires a fresh executor-owned reservation")
    error = None
    try:
        refine(source, output / "common-repair", platform="CPU")
    except Exception:
        error = traceback.format_exc()
        save(output / "repair-failure.json", {"traceback": error, "scientific_failure_retained": True})
    final = output / "common-repair/refined.pdb"
    repaired = final.is_file()
    if not repaired:
        final = source / "source-restored.pdb"
    save(output / "result.json", {"final_path": str(final), "final_sha256": digest(final), "repaired": repaired,
        "repair_failure": error, "input_sha256": digest(source / "source-restored.pdb"), "resampled": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        worker(args.audit.resolve(), args.output.resolve())
