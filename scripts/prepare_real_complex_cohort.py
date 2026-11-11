#!/usr/bin/env python3
"""Freeze the approved ten-complex cc94 MPS cohort without neural sampling.

Sources are experimental PDB downloads, never artificially masked structures.
Uses the shipped assembly, altloc, sequence, input preparation and admission
services. A refusal is retained, not repaired by changing scientific policy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from collections import Counter
from dataclasses import asdict
from pathlib import Path

SEEDS = (7, 101, 102, 103, 104)
PROFILE = "protpardelle-1c-cc94-mps-v1"
# Same archived nine originals as the pre-anchor study, plus deposited 9CQD.
SOURCES = {
    "9CQD": "cab9d5dee612be94ab235cc916ef0f1e852e5102ad4c889eda9d4bbe8ed77b9c",
    "1N8Z": "ed934d36503335c6764b08880f9e50933d99eb13b9bf145339cfcaee150eb21e",
    "1JPS": "22cb05ba2059722d2a571337de9b4f0ef6adc8946b13d16e5997d972a9210e14",
    "1AHW": "bc73f266b34041673df285e1e0c66a507d9452642ff64451424167d32e117cdf",
    "1G9M": "b1d4b3cfb26ca9c5bb2ac5c83da0a3255e047ba7812e1f0a1ca285a58396541d",
    "1GC1": "73fb7c47f6b4c5ce0fe51888cf5e32bd5ee3de96bcc57e2fbaf10692119bff60",
    "1BRS": "0076aae54252706a8c8c02f7b1e71ce2394b0e928b413590b4ac12aade9e8f6c",
    "1I2M": "d11d60e277912fc62f2379ae4d2a64fb5d7970dbc9742197a11ea6f4cf4a096c",
    "1KTZ": "54adccd77517b42c9d6019afe4eee8c0a6aa9e7fd6966a0c1121eda7b3e9426e",
    "1HE8": "1814c7e5e4a5104301564d41d01d82482958bf20b46ef9eab16021feebc4b489",
}


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path: Path, value: object) -> None:
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def prepare_case(root: Path, pdb: str) -> dict:
    from dvbfixer.altloc_selection import select_altlocs
    from dvbfixer.biological_assembly import parse_biological_assemblies, render_biological_assembly
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.model.diffusion.input_geometry import assess_input_geometry
    from dvbfixer.model.diffusion.input_preparation import prepare_input
    from dvbfixer.model.diffusion.protein_context import (
        prepare_cc94_request,
        validate_context_capability,
    )
    from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy
    from dvbfixer.model.diffusion.scope import assess_diffusion_scope
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request
    from dvbfixer.model.pipeline import _strip_hetatm_lines

    directory = root / "cases" / pdb
    directory.mkdir(parents=True, exist_ok=False)
    source = root / "sources" / f"{pdb}.pdb"
    row = {"pdb": pdb, "seeds": list(SEEDS), "sampling_started": False}
    try:
        if digest(source) != SOURCES[pdb]:
            raise ValueError("source differs from the prospectively frozen experimental download")
        entry_path = root / "sources" / f"{pdb}.entry.json"
        entry = json.loads(entry_path.read_bytes())
        if entry["rcsb_id"] != pdb or not any(e["method"] == "X-RAY DIFFRACTION" for e in entry["exptl"]):
            raise ValueError("entry identity or experimental method does not match")
        lines = source.read_text().splitlines(keepends=True)
        assemblies = parse_biological_assemblies(lines)
        assembly = render_biological_assembly(lines, assemblies["1"],
            list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"), renumber=False)
        assembly_text = "".join(assembly)
        (directory / "assembly1.pdb").write_text(assembly_text)
        selected, report = select_altlocs(assembly_text,
            source_sha256=hashlib.sha256(assembly_text.encode()).hexdigest(), policy="occupancy")
        save(directory / "altlocs.json", report)
        if selected is None:
            raise ValueError("residue-coherent altloc selection did not produce a structure")
        (directory / "selected.pdb").write_text(selected)
        removed = Counter(l[17:20].strip() for l in selected.splitlines() if l.startswith("HETATM"))
        protein = "".join(_strip_hetatm_lines(selected.splitlines(keepends=True), keep_water=False))
        model_input = directory / "protein.pdb"
        model_input.write_text(protein)
        targets = _target_sequences(model_input, None)
        if not targets:
            raise ValueError("biological assembly lost deposited SEQRES")
        (directory / "target.fasta").write_text("".join(f">{c}\n{s}\n" for c, s in targets.items()))
        row.update(original_sha256=digest(source), entry_sha256=digest(entry_path),
            title=entry["struct"]["title"], assembly="1", removed_heterogen_atom_counts=dict(removed),
            target_lengths={c: len(s) for c, s in targets.items()},
            protein_sha256=digest(model_input), source_derivative="explicit protein-only biological assembly")
        request = build_cli_diffusion_request(model_input, targets, seeds=SEEDS, profile=PROFILE,
            no_terminal=False, refinement_options=RefinementPolicy().backend_options())
        (directory / "request-before-preparation.json").write_text(request.to_json())
        row["gaps"] = [asdict(gap) for gap in request.gaps]
        row["gap_count"] = len(request.gaps)
        # REMARK 465 remains evidence, not an alternative greedy sequence mapper.
        save(directory / "deposited-missing-records.json",
             [line.rstrip() for line in lines if line.startswith("REMARK 465")])
        prepared = prepare_input(request, model_input.read_bytes(), directory)
        if prepared is not None:
            request, model_input = prepared
        row["prepared_input"] = str(model_input.relative_to(root))
        request = prepare_cc94_request(request, model_input.read_bytes(), "auto")
        validate_context_capability(request, "protpardelle")
        scope = assess_diffusion_scope(request, model_input.read_bytes())
        save(directory / "scope.json", asdict(scope))
        if not scope.supported:
            raise ValueError(scope.message)
        findings = assess_input_geometry(request, model_input.read_bytes())
        save(directory / "source-findings.json", [asdict(f) for f in findings])
        if any(f.severity is Severity.ERROR for f in findings):
            raise ValueError("severe source findings remain after supported input preparation")
        request_path = model_input.parent / "request-frozen.json"
        with request_path.open("x") as stream:
            stream.write(request.to_json())
        row.update(status="admitted", request=str(request_path.relative_to(root)),
            request_sha256=digest(request_path), fixed_atoms=len(request.fixed_atoms),
            generated_atoms=len(request.generated_atoms),
            minimum_denoiser_invocations=len(request.gaps) * len(SEEDS))
    except Exception as exc:
        row.update(status="admission_refused", error_type=type(exc).__name__, error=str(exc))
        (directory / "preparation-error.txt").write_text(traceback.format_exc())
    save(directory / "preparation.json", row)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--case", choices=tuple(SOURCES))
    parser.add_argument("--case-timeout", type=int, default=1800)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.case:
        print(json.dumps(prepare_case(root, args.case), ensure_ascii=False, default=str), flush=True)
        return
    repo = Path(__file__).resolve().parents[1]
    protocol = {"protocol": "real-complex-cc94-mps-five-seed-v1", "seeds": list(SEEDS),
        "structures": list(SOURCES), "profile": PROFILE, "source_sha256": SOURCES,
        "assembly": "1", "altloc_policy": "occupancy", "protein_only": True,
        "preserve_all_assembly_protein_chains": True, "full_deposited_seqres": True,
        "artificial_deletions": False, "no_terminal": False, "context": "auto",
        "sampling": {"updates": 500, "step_scale": 1.0, "s_churn": 0},
        "refinement": {"platform": "CPU", "solvent": False, "production_defaults": True},
        "anchor_release": "production geometry-gated; no forced release or threshold tuning",
        "input_preparation": "shipped observed-heavy completion and bounded repair before freezing",
        "seed_replacement": False, "qualification": "case cohort, not accuracy or training-disjoint benchmark",
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
        "case_timeout_seconds": args.case_timeout,
        "script_sha256": digest(Path(__file__)),
        "source_code_sha256": {str(p.relative_to(repo)): digest(p)
            for p in [*(repo / "src/dvbfixer/model/diffusion").glob("*.py"),
                      repo / "src/dvbfixer/model/diffusion_cli.py",
                      repo / "deploy/protpardelle-1c/cc94_runner.py"]}}
    save(root / "protocol.json", protocol)
    shutil.copyfile(__file__, root / "prepare-script.py")
    rows = []
    for pdb in SOURCES:
        started = time.monotonic()
        with (root / f"prepare-{pdb}.log").open("xb") as log:
            try:
                result = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                    "--root", str(root), "--case", pdb], stdout=log, stderr=subprocess.STDOUT,
                    timeout=args.case_timeout, env={**os.environ, "OPENMM_CPU_THREADS": "1"})
                report = root / "cases" / pdb / "preparation.json"
                row = json.loads(report.read_text()) if report.exists() else {
                    "pdb": pdb, "status": "preparation_process_failed", "exit_code": result.returncode}
            except subprocess.TimeoutExpired:
                row = {"pdb": pdb, "status": "preparation_timeout", "sampling_started": False}
        row["preparation_elapsed_seconds"] = round(time.monotonic() - started, 3)
        rows.append(row)
        # Per-case files are immutable; only this in-progress index changes.
        temporary = root / "preparation-summary.next.json"
        temporary.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n")
        temporary.replace(root / "preparation-summary.json")
        print(pdb, row["status"], row.get("gap_count"), row.get("error", ""), flush=True)
    print("Admitted:", sum(r["status"] == "admitted" for r in rows), "/", len(rows), flush=True)


if __name__ == "__main__":
    main()
