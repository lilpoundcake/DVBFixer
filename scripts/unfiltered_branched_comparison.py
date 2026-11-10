"""Plan-before-run controls selected without reconstruction outcomes.

Three explicit chain-A ATOM projections, first two VAL and first two THR in file
order, four masks per residue. Reuse the prior comparator without changing gates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict
from pathlib import Path

from scripts.branched_template_baseline import preflight, template_proposal
from scripts.compare_branched_completion import measure, score

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.anchored_branch import complete_anchored_branch
from dvbfixer.protein_completion.anchored_source import plan_anchored_source
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition, CompletionRequest
from dvbfixer.protein_completion.service import complete_sidechain

FIXTURES = {"1EMV": "tests/fixtures/regressions/1EMV.pdb",
            "1TIM": "tests/fixtures/protein_context/1TIM.pdb",
            "1DQJ": "tests/fixtures/regressions/1DQJ.pdb"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+"\n").encode()


def masks(residue_name: str) -> tuple[tuple[str, ...], ...]:
    branch = BRANCHES[residue_name]
    return ((branch.first, branch.second), ("CB", branch.first, branch.second), (branch.first,), (branch.second,))


def select_targets(lines: list[str]) -> tuple[tuple[str, ComponentInstanceRef], ...]:
    """File-order identity selection only: no geometry/occupancy/own-pass filtering."""
    seen = set()
    selected = []
    counts = {"VAL": 0, "THR": 0}
    for line in lines:
        if not line.startswith("ATOM  ") or line[21] != "A":
            continue
        name = line[17:20]
        if name not in counts:
            continue
        ref = ComponentInstanceRef(1, line[21], line[22:26].strip(), line[26].strip())
        if ref in seen:
            continue
        seen.add(ref)
        if counts[name] < 2:
            selected.append((name, ref))
            counts[name] += 1
    if counts != {"VAL": 2, "THR": 2}:
        raise ValueError("declared structure lacks two observed identities of each residue type")
    return tuple(selected)


def plan(root: Path, destination: Path) -> None:
    """Materialize/commit every masked source before ANY reconstruction is run."""
    destination.mkdir(parents=True, exist_ok=False)
    cases = []
    artifacts = {}

    def save(name, data):
        (destination/name).write_bytes(data)
        artifacts[name] = sha(data)

    for label, filename in FIXTURES.items():
        original = (root/filename).read_bytes()
        save(f"{label}-original.pdb", original)
        lines = [l for l in original.decode("ascii").splitlines(keepends=True)
                 if l.startswith("ATOM  ") and l[21] == "A"]
        reference = "".join(lines)+"TER\nEND\n"
        save(f"{label}-reference.pdb", reference.encode())
        for name, target in select_targets(lines):
            for mask in masks(name):
                identity = f"{label}-{name}{target.sequence_number}{target.insertion_code}-{'-'.join(mask)}"
                source = "".join(l for l in lines if not (
                    l[22:26].strip() == target.sequence_number and l[26].strip() == target.insertion_code
                    and l[12:16].strip() in mask))+"TER\nEND\n"
                save(identity+".pdb", source.encode())
                cases.append({"id": identity, "fixture": label, "target": asdict(target), "residue_name": name,
                              "mask": mask, "mode": "anchored" if len(mask) == 1 else "paired",
                              "source": identity+".pdb", "reference": f"{label}-reference.pdb"})
    files = [Path(__file__).resolve(), root/"scripts/compare_branched_completion.py",
             root/"scripts/branched_template_baseline.py", root/"scripts/branched_stereo_oracle.py",
             *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    record = {"schema": "dvbfixer.unfiltered-branched-plan.v1", "selection": "first two VAL and first two THR in chain-A file order; all four masks",
              "projection": "Every original chain-A ATOM line, including altlocs/partial occupancies; bare TER/END. Not intact complexes.",
              "case_count": len(cases), "cases": cases, "artifacts": artifacts,
              "source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in files},
              "common_gate": "unchanged prior measured-chi gate; illustrative angles, not a calibrated gold standard",
              "template": "shared seeded helper on Reference; no reflection/downstream refinement",
              "independence_limit": "1EMV is homologous to previously tested 1FR2; no family-level generalization",
              "no_own_success_filter": True, "production_enabled": False}
    (destination/"plan.json").write_bytes(json_bytes(record))
    print(f"Frozen {len(cases)} cases before proposal execution: {destination/'plan.json'}")


def reporting_scores(request: CompletionRequest, additions: tuple[Addition, ...], truth: dict) -> dict:
    """Reporting only; symmetry may permute exclusively generated VAL methyls."""
    result = score(additions, truth)
    names = {a.atom.name for a in additions}
    can_swap = request.residue_name == "VAL" and {"CG1", "CG2"} <= names
    sym = result["exact_name_rmsd_angstrom"]
    swapped = False
    if can_swap:
        other = dict(truth)
        other["CG1"], other["CG2"] = truth["CG2"], truth["CG1"]
        alternative = score(additions, other)["exact_name_rmsd_angstrom"]
        swapped = alternative < sym
        sym = min(sym, alternative)
    points = {o.atom.name: o.xyz for o in request.observations if o.atom.component == request.residue}
    points.update({a.atom.name: a.xyz for a in additions})
    primary = BRANCHES[request.residue_name].first
    chi = signed_dihedral_degrees(*(points[n] for n in ("N", "CA", "CB", primary)))
    ref_chi = signed_dihedral_degrees(*(truth[n] for n in ("N", "CA", "CB", primary)))
    return {**result, "generated_only_val_symmetry_rmsd_angstrom": sym,
            "val_swap_admissible": can_swap, "val_swapped_in_reporting": swapped,
            "named_chi1_error_degrees": abs((chi-ref_chi+180) % 360-180) if chi is not None and ref_chi is not None else None,
            "reporting_only_no_reranking": True}


def run(plan_root: Path, destination: Path, summary: Path) -> None:
    if plan_root.is_symlink():
        raise ValueError("regular plan root required")
    plan_root = plan_root.resolve(strict=True)
    data = read_private_file_bytes(plan_root, "plan.json", max_bytes=2_000_000)
    record = json.loads(data)
    if record["schema"] != "dvbfixer.unfiltered-branched-plan.v1" or record["case_count"] != 48 or len(record["cases"]) != 48:
        raise ValueError("complete frozen 48-case plan required")
    root = Path(__file__).resolve().parents[1]
    for name, digest in record["source_hashes"].items():
        if sha((root/name).read_bytes()) != digest:
            raise ValueError("source implementation changed after plan freezing")
    payloads = {}
    for name, digest in record["artifacts"].items():
        if Path(name).name != name:
            raise ValueError("flat plan artifacts required")
        value = read_private_file_bytes(plan_root, name, max_bytes=16*1024*1024)
        if sha(value) != digest:
            raise ValueError("frozen plan input changed")
        payloads[name] = value
    environment = preflight()  # Environment failure is not a science rejection.
    destination.mkdir(parents=True, exist_ok=False)
    report = {"schema": "dvbfixer.unfiltered-branched-results.v1", "plan_sha256": sha(data),
              "plan_root": str(plan_root), "retained_directory": str(destination), "environment": environment,
              "python": platform.python_version(), "results": {}}
    (destination/"executed-plan.json").write_bytes(data)
    for case in record["cases"]:
        folder = destination/case["id"]
        folder.mkdir()
        source = payloads[case["source"]]
        (folder/"input.pdb").write_bytes(source)
        target = ComponentInstanceRef(**case["target"])
        anchored = case["mode"] == "anchored"
        try:
            source_plan = plan_anchored_source(source, target) if anchored else plan_branched_source(source, target)
        except ValueError as exc:
            row = {"status": "shared_admission_refused", "reason": str(exc), "own_executed": False, "template_executed": False}
        else:
            request = source_plan.request
            (folder/"source-plan.json").write_bytes(json_bytes(asdict(source_plan)))
            own = complete_anchored_branch(request) if anchored else complete_sidechain(request)
            (folder/"own-attempts.json").write_bytes(json_bytes(asdict(own)))
            row = {"status": "admitted", "own_executed": True, "template_executed": False,
                   "request_digest": request.digest, "own": {"status": own.status, "reason": own.reason}}
            if own.selected is not None:
                audit, own_pdb = measure(source, request, own.selected.additions, anchored=anchored)
                (folder/"own-proposal.pdb").write_bytes(own_pdb)
                row["own"]["common_audit"] = asdict(audit)

            def retain(snapshot, out=folder):
                (out/f"template-seed-{snapshot.seed}.json").write_bytes(json_bytes(asdict(snapshot)))

            try:
                with (folder/"template.log").open("x") as log, redirect_stdout(log), redirect_stderr(log):
                    template = template_proposal(source, request, retain_snapshot=retain)
                row["template_executed"] = True
                (folder/"template-proposal.json").write_bytes(json_bytes(asdict(template)))
                audit, template_pdb = measure(source, request, template.additions, anchored=anchored)
                (folder/"template-proposal.pdb").write_bytes(template_pdb)
                row["template"] = {"status": "evaluated", "common_audit": asdict(audit), "helper_attempts": len(template.snapshots),
                                   "temporary_observed_motion_angstrom": template.snapshots[-1].observed_max_movement_angstrom}
            except Exception as exc:
                row["template"] = {"status": "execution_or_protocol_blocker", "error_type": type(exc).__name__, "error": str(exc)}
                (folder/"execution-error.json").write_bytes(json_bytes(row["template"]))
            else:
                # No truth is supplied to either proposal. Scoring begins only now.
                truth = {a.name: a.position for a in _parse_atom_lines(payloads[case["reference"]].decode().splitlines(keepends=True))
                         if _component_ref(a) == target}
                row["template"]["scores"] = reporting_scores(request, template.additions, truth)
                if own.selected is not None:
                    row["own"]["scores"] = reporting_scores(request, own.selected.additions, truth)
        (folder/"comparison.json").write_bytes(json_bytes(row))
        row["retained_files"] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir())}
        report["results"][case["id"]] = row
    raw = json_bytes(report)
    (destination/"evidence.json").write_bytes(raw)
    with summary.open("xb") as stream:
        stream.write(json_bytes(report | {"raw_record_sha256": sha(raw)}))
    print(json.dumps({"planned": 48, "shared_admission_refused": sum(r["status"] == "shared_admission_refused" for r in report["results"].values()),
                      "own_common_passes": sum(r.get("own", {}).get("common_audit", {}).get("accepted", False) for r in report["results"].values()),
                      "template_common_passes": sum(r.get("template", {}).get("common_audit", {}).get("accepted", False) for r in report["results"].values()),
                      "template_execution_blockers": sum(r.get("template", {}).get("status") == "execution_or_protocol_blocker" for r in report["results"].values())}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    make = sub.add_parser("plan")
    make.add_argument("--output", type=Path, required=True)
    execute = sub.add_parser("run")
    execute.add_argument("--plan", type=Path, required=True)
    execute.add_argument("--output", type=Path, required=True)
    execute.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "plan":
        plan(Path(__file__).resolve().parents[1], args.output)
    else:
        run(args.plan, args.output, args.summary)


if __name__ == "__main__":
    main()
