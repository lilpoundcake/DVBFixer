"""Independent readback of private hybrid evidence; never repairs coordinates."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def audit(root):
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    if (root / "prepared-bindings.json").is_file():
        prepared = json.loads((root / "prepared-bindings.json").read_text())
        if prepared["plan_sha256"] != digest(plan_path):
            raise ValueError("prepared inputs bind a different plan")
        verify_bindings(prepared["files"])
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    if digest(root / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("archived source does not match the frozen request")
    source, proposal = atoms(root / "source.pdb"), atoms(root / "proposal.pdb")
    fixed = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.fixed_atoms}
    generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    if set(source) != fixed or set(proposal) != fixed | generated:
        raise ValueError("archived ownership mismatch")
    assess = assessment_adapter(request, (root / "source.pdb").read_text())
    summary = json.loads((root / "summary.json").read_text())
    rows = []
    for row in summary["candidates"]:
        path = root / row["candidate"]
        if digest(path) != row["sha256"]:
            raise ValueError("candidate bytes changed")
        final = atoms(path)
        if set(final) != fixed | generated:
            raise ValueError("candidate changed inventory")
        if any(final[a]["line"] != source[a]["line"] for a in fixed):
            raise ValueError("observed record changed")
        old_lines = (root / "proposal.pdb").read_text().splitlines(keepends=True)
        new_lines = path.read_text().splitlines(keepends=True)
        if len(old_lines) != len(new_lines):
            raise ValueError("physical records changed")
        for before, after in zip(old_lines, new_lines, strict=True):
            a = (before[21], before[22:26].strip(), before[26].strip(), before[12:16].strip()) if before.startswith("ATOM  ") else None
            if a in generated:
                if before[:30] != after[:30] or before[54:] != after[54:]:
                    raise ValueError("generated metadata changed")
            elif before != after:
                raise ValueError("non-generated record changed")
        native_rows = json.loads((path.parent / "native.json").read_text())["atoms"]
        native = {tuple(r["identity"]): r["xyz"] for r in native_rows}
        if len(native) != len(native_rows) or set(native) != fixed | generated:
            raise ValueError("native atom inventory changed")
        if any(len(xyz) != 3 or not all(math.isfinite(v) for v in xyz) for xyz in native.values()):
            raise ValueError("nonfinite native coordinates")
        maximum = max(math.dist(source[a]["xyz"], native[a]) for a in fixed)
        if not math.isfinite(maximum) or maximum > 1e-7:
            raise ValueError("native fixed atoms moved before serialization")
        for a in generated:
            expected = tuple(float(f"{v:8.3f}") for v in native[a])
            if expected != final[a]["xyz"]:
                raise ValueError("serialized generated coordinates differ from native proposal")
        result = assess(path.read_text())
        if result.geometry_passed != row["geometry_passed"] or list(result.admission_failures) != row["admission_failures"]:
            raise ValueError("assessment does not reproduce")
        rows.append({"candidate": row["candidate"], "fixed_records_exact": True,
                     "native_fixed_max_A": maximum, "generated_count": len(generated),
                     "geometry_passed": result.geometry_passed})
    return {"plan_sha256": digest(plan_path), "candidates": rows,
            "new_diffusion_calls": 0, "experimental_control_only": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    report = audit(args.root.resolve())
    save(args.root / "independent-readback.json", report)
    print(json.dumps(report))
