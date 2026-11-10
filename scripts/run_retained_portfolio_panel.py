#!/usr/bin/env python3
"""Predeclare and execute a fixed-budget P2 workflow on historical saved donors.

This prospectively declares execution, NOT a blinded prospective scientific
cohort. The separately versioned corrected 1N8Z arm is never pooled silently.
No new sampling, threshold changes, P3 reconstruction, or baseline replacement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    args.output = args.output.resolve(strict=True)
    repository = Path(__file__).resolve().parents[1]
    scripts = repository / "scripts"
    cohort = args.cohort.resolve(strict=True)
    original_bytes = (cohort / "cohort-final-outcomes.json").read_bytes()
    original = json.loads(original_bytes)
    cases = [row["pdb"] for row in original["rows"]]
    if len(cases) != 10 or len(set(cases)) != 10:
        raise ValueError("retained panel requires the declared ten unique historical cases")
    numerical_sources = [repository / "src/dvbfixer/model/diffusion" / name for name in
        ("gap_portfolio.py", "boundary_refinement.py", "geometry_conditioning.py", "portfolio_contract.py",
         "portfolio_evidence.py", "portfolio_publication.py")]
    numerical_sources.extend(scripts / name for name in ("research_gap_portfolio.py", "repair_gap_portfolio.py",
        "audit_gap_portfolio.py", "publish_gap_portfolio.py"))
    numerical_sources.append(Path(__file__).resolve())
    digests = {str(p.relative_to(repository)): hashlib.sha256(p.read_bytes()).hexdigest() for p in numerical_sources}
    declaration = {"scope": __doc__, "cases": cases,
        "original_cohort_sha256": hashlib.sha256(original_bytes).hexdigest(),
        "software_arms": {row["pdb"]: row["software_arm"] for row in original["rows"]},
        "selection": {"max_evaluations": 64, "beam_width": 4},
        "optional_issue_local_repair": {"refinement_seed": 7, "max_iterations": 250, "restart_count": 4,
            "platform": "CPU", "only_when_full_selected_assessment_fails": True},
        "new_neural_samples": 0, "P3": False, "source_sha256": digests}
    (args.output / "declaration.json").write_text(json.dumps(declaration, indent=2) + "\n")
    (args.output / "run-script.py").write_bytes(Path(__file__).read_bytes())
    ledger = []
    for case in cases:
        for path, expected in digests.items():
            if hashlib.sha256((repository / path).read_bytes()).hexdigest() != expected:
                raise RuntimeError("retained panel software changed during execution")
        directory = args.output / case
        directory.mkdir()
        row = {"case": case, "software_arm": declaration["software_arms"][case], "status": "running"}
        ledger.append(row)
        (args.output / "ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
        try:
            with (directory / "execution.log").open("x") as log:
                def run(name, *arguments, stream=log):
                    subprocess.run([sys.executable, str(scripts / name), *map(str, arguments)],
                        stdout=stream, stderr=subprocess.STDOUT, check=True, cwd=repository)

                portfolio = directory / "selection"
                run("research_gap_portfolio.py", "--cohort", cohort, "--case", case,
                    "--output", portfolio, "--max-evaluations", 64, "--beam-width", 4)
                selected = json.loads((portfolio / "result.json").read_text())["selected"]["assessment"]
                needs_repair = bool(selected["admission_failures"] or any(i["actionable"] for i in selected["issues"]))
                repair = directory / "repair" if needs_repair else None
                if repair is not None:
                    run("repair_gap_portfolio.py", "--portfolio", portfolio, "--output", repair,
                        "--refinement-seed", 7, "--max-iterations", 250, "--restart-count", 4)
                arguments = ["--portfolio", portfolio, "--report", directory / "readback.json"]
                if repair is not None:
                    arguments.extend(["--repair", repair])
                run("audit_gap_portfolio.py", *arguments)
                audit = json.loads((directory / "readback.json").read_text())
                row.update(private_geometry_and_custody_passed=audit["passed"], issue_local_repair=needs_repair,
                    final_sha256=audit["final_sha256"], status="geometry-rejected")
                if audit["passed"]:
                    arguments = ["--portfolio", portfolio, "--cohort", cohort, "--case", case,
                        "--work", directory / "publication-work", "--bundle", directory / "bundle"]
                    if repair is not None:
                        arguments.extend(["--repair", repair])
                    run("publish_gap_portfolio.py", *arguments)
                    row.update(status="private-bundle-verified", publication=json.loads(
                        (directory / "publication-work/summary.json").read_text()))
                row["original_outcomes_replaced"] = False
        except Exception as exc:
            row.update(status="execution-error", error_type=type(exc).__name__, error=str(exc))
            (args.output / "ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
            raise
        (args.output / "ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
        print(json.dumps(row), flush=True)
    (args.output / "summary.json").write_text(json.dumps({"scope": __doc__, "rows": ledger,
        "original_cohort_unchanged": hashlib.sha256((cohort / "cohort-final-outcomes.json").read_bytes()).hexdigest()
            == declaration["original_cohort_sha256"], "new_neural_samples": 0,
        "scientific_qualification": False, "public_profile_enabled": False}, indent=2) + "\n")


if __name__ == "__main__":
    main()
