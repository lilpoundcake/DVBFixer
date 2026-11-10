"""Retrospective matched controls: own method versus shared seeded PDBFixer.

Both arms are judged by the same measured-chi geometry gate, never by the own
method's staggered grid. Reference coordinates enter reporting only after both
proposals. This is not full legacy Modeller/prepare/minimize comparison.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from scripts.branched_stereo_oracle import rdkit_labels
from scripts.branched_template_baseline import preflight, template_proposal

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.anchored_branch import complete_anchored_branch
from dvbfixer.protein_completion.anchored_source import plan_anchored_source
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.branched_source import _thr_beta_checks, plan_branched_source
from dvbfixer.protein_completion.contracts import Addition, Candidate, CompletionRequest
from dvbfixer.protein_completion.evidence import _targets
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import (
    _global_ca,
    _insert_candidates,
    _layout,
    _readback_additions,
)
from dvbfixer.protein_completion.validation import _audit_candidate


@dataclass(frozen=True)
class ComparisonAudit:
    accepted: bool
    raw: dict
    rounded: dict
    global_stereo_findings: tuple[str, ...]
    rdkit_target_cip: dict[str, str]
    observed_preserved: bool


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def measure(source: bytes, request: CompletionRequest, additions: tuple[Addition, ...], *, anchored: bool):
    """Research-only common geometry, fixed-source and PDB-rounding audit.

    The returned PDB is diagnostic until audit.accepted, not a production
    materialization bypass. Source-bound planning must have admitted this input.
    """
    points = {o.atom.name: o.xyz for o in request.observations if o.atom.component == request.residue}
    points.update({a.atom.name: a.xyz for a in additions})
    branch = BRANCHES[request.residue_name]
    chi = signed_dihedral_degrees(*(points[n] for n in ("N", "CA", "CB", branch.first)))
    candidate = Candidate(request.digest, additions, chi)
    raw = _audit_candidate(request, candidate, anchored_branch=anchored, measured_chi=True)
    output = _insert_candidates(source, (request,), (candidate,))
    rounded_atoms = _readback_additions(source, output, additions)
    rounded = _audit_candidate(request, replace(candidate, additions=rounded_atoms),
                               anchored_branch=anchored, measured_chi=True, pdb_roundtrip=True)
    _, atoms, _ = _layout(output)
    findings = []
    for check in (lambda: _global_ca(atoms, allow_missing_cb=None),
                  lambda: _thr_beta_checks(atoms, pending_target=None)):
        try:
            check()
        except ValueError as exc:
            findings.append(str(exc))
    xyz = {a.name: a.position for a in atoms if _component_ref(a) == request.residue}
    cip = rdkit_labels(xyz, request.residue_name)
    expected = {"CA": "S", "CB": "R"} if request.residue_name == "THR" else {"CA": "S"}
    if cip != expected:
        findings.append("independent_target_CIP_mismatch")
    return ComparisonAudit(raw.accepted and rounded.accepted and not findings,
                            asdict(raw), asdict(rounded), tuple(findings), cip, True), output


def score(additions: tuple[Addition, ...], truth: dict[str, tuple[float, float, float]]) -> dict:
    errors = [math.dist(a.xyz, truth[a.atom.name]) for a in additions]
    if not errors:
        raise ValueError("nonempty generated coordinate set required for comparison")
    return {"exact_name_rmsd_angstrom": math.sqrt(sum(e*e for e in errors)/len(errors)),
            "maximum_error_angstrom": max(errors), "atom_count": len(errors),
            "symmetry_minimized": False, "used_for_selection": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paired", type=Path, required=True)
    parser.add_argument("--anchored", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    # An unavailable stack is an environment blocker: stop before science rather
    # than count it as a failed template proposal or add an optional broad skip.
    environment = preflight()
    cases = []
    for mode, retained in (("paired", args.paired), ("anchored", args.anchored)):
        record = json.loads((retained/"evidence.json").read_bytes())
        declarations = {row["id"]: row for row in record["cases"]}
        for name, result in record["results"].items():
            if result["status"] == "accepted":
                cases.append({"id": mode+"-"+name, "mode": mode, "retained": str(retained),
                              "retained_case": name, "declaration": declarations[name],
                              "input_sha256": result["retained_files"]["input.pdb"],
                              "prior_output_sha256": result["retained_files"]["completed.pdb"]})
    if len(cases) != 32:
        raise ValueError("expected the retained 16 paired plus 16 singleton cases")
    plan = {"schema": "dvbfixer.retrospective-branched-template-comparison.v1", "cases": cases,
            "environment": environment, "case_selection": "all 32 previously accepted source controls; retrospective, not blind",
            "template": "shared seeded rebuild helper 1..5 on OpenMM Reference; preserve all returned states",
            "common_gate": "same physical geometry/stereo/clash/rounded checks, measured chi; no own staggered-grid requirement",
            "observations": "unchanged complete masked source for both arms; no reference coordinates in proposals",
            "no_refinement_or_reflection": True, "not_full_legacy_pipeline": True}
    (args.output/"plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    root = Path(__file__).resolve().parents[1]
    files = [Path(__file__).resolve(), root/"scripts/branched_stereo_oracle.py", root/"scripts/branched_template_baseline.py",
             *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    report = plan | {"retained_directory": str(args.output), "python": platform.python_version(),
                     "source_hashes": {str(f.relative_to(root)): sha(f.read_bytes()) for f in files}, "results": {}}
    for case in cases:
        folder = args.output/case["id"]
        folder.mkdir()
        source = (Path(case["retained"])/case["retained_case"]/"input.pdb").read_bytes()
        if sha(source) != case["input_sha256"]:
            raise ValueError("retained masked input changed")
        (folder/"input.pdb").write_bytes(source)
        # Target identity comes from the saved source request, independently
        # rebound by the current planner; held-out input is not supplied here.
        saved = json.loads((Path(case["retained"])/case["retained_case"]/"source-plan.json").read_bytes())
        target = _targets(json.dumps([saved["request"]["residue"]]).encode())[0]
        anchored = case["mode"] == "anchored"
        admitted = plan_anchored_source(source, target) if anchored else plan_branched_source(source, target)
        request = admitted.request
        (folder/"source-plan.json").write_text(json.dumps(asdict(admitted), indent=2)+"\n")
        own = complete_anchored_branch(request) if anchored else complete_sidechain(request)
        (folder/"own-attempts.json").write_text(json.dumps(asdict(own), indent=2)+"\n")
        if own.selected is None:
            raise ValueError("previously accepted own control no longer selects a proposal")
        own_audit, own_pdb = measure(source, request, own.selected.additions, anchored=anchored)
        (folder/"own-proposal.pdb").write_bytes(own_pdb)
        if sha(own_pdb) != case["prior_output_sha256"]:
            raise ValueError("own comparison changed the prior accepted output")

        def retain(snapshot, destination=folder):
            (destination/f"template-seed-{snapshot.seed}.json").write_text(json.dumps(asdict(snapshot), indent=2, allow_nan=False)+"\n")

        try:
            with (folder/"template.log").open("x") as log, redirect_stdout(log), redirect_stderr(log):
                template = template_proposal(source, request, retain_snapshot=retain)
            (folder/"template-proposal.json").write_text(json.dumps(asdict(template), indent=2, allow_nan=False)+"\n")
            template_audit, template_pdb = measure(source, request, template.additions, anchored=anchored)
            (folder/"template-proposal.pdb").write_bytes(template_pdb)
        except Exception as exc:
            (folder/"execution-error.json").write_text(json.dumps({"type": type(exc).__name__, "error": str(exc)}, indent=2)+"\n")
            raise  # An execution/protocol failure is not a scientific rejection.
        # Only after both proposals and audits: held-out coordinates for reporting.
        original = Path(case["retained"])/f"original-{case['declaration']['fixture']}.pdb"
        truth = {a.name: a.position for a in _parse_atom_lines(original.read_text().splitlines(keepends=True))
                 if _component_ref(a) == target}
        row = {"request_digest": request.digest, "own": {"common_audit": asdict(own_audit),
                   "scores": score(own.selected.additions, truth)},
               "template": {"common_audit": asdict(template_audit), "scores": score(template.additions, truth),
                   "helper_attempts": len(template.snapshots),
                   "observed_temporary_max_movement_angstrom": template.snapshots[-1].observed_max_movement_angstrom,
                   "unpublished_extras": [asdict(a) for a in template.temporary_extra_atoms]},
               "output_label": "diagnostic proposals, only common_audit.accepted denotes scientific gate pass"}
        (folder/"comparison.json").write_text(json.dumps(row, indent=2, allow_nan=False)+"\n")
        row["retained_files"] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir())}
        report["results"][case["id"]] = row
    data = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(data)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(data)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({arm: {"common_passes": sum(c[arm]["common_audit"]["accepted"] for c in report["results"].values()),
                           "total": len(cases)} for arm in ("own", "template")}, indent=2))


if __name__ == "__main__":
    main()
