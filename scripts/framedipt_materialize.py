"""Private generated-only sidechain construction for a retained FrameDiPT endpoint.

Native backbone drift remains retained; original observed records are restored
explicitly for this separate materialization, never claimed natively preserved.
No closure, heavy-atom minimization, source repair or sequence design occurs.
"""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_anchor_calibration_grid import diagnose
from scripts.audit_framedipt_research import audit
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import audit_native, diagnostic_counts, key, launch_worker
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, ATOM37_ORDER
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    args.owned_output = True
    output, sample = args.output.resolve(), args.sample.resolve()
    parent = audit(sample)
    if not parent["sampling_completed"]:
        raise ValueError("materialization requires a completed audited endpoint")
    request = DiffusionRequest.from_json((sample / "request.json").read_text())
    source = atoms(sample / "source.pdb")
    source_text = (sample / "source.pdb").read_text()
    endpoint = dict(np.load(sample / "native-endpoint.npz", allow_pickle=False))
    axis = target_axis(request)
    order = [r for t in request.target_sequences for r in axis[t.chain]]
    generated = set(map(key, request.generated_atoms))
    names = {(r.chain, r.residue_number, r.insertion_code): _CANONICAL_ONE_TO_THREE[letter]
             for t in request.target_sequences for r, letter in zip(axis[t.chain], t.sequence, strict=True)}
    lines = []
    serial = max(int(row["line"][6:11]) for row in source.values())
    frame_derived = []
    for i, residue in enumerate(order):
        rk = residue.chain, residue.residue_number, residue.insertion_code
        existing = [a for a in source if a[:3] == rk]
        if existing:
            lines.extend(source[a]["line"] for a in existing)
        else:
            for name in ("N", "CA", "C", "O", "CB"):
                a = (*rk, name)
                if a not in generated:
                    continue
                frame_derived.append(a)
                serial += 1
                xyz = endpoint["backbone"][i, ATOM37_ORDER.index(name)]
                lines.append(f'ATOM  {serial:5d} {name:^4s} {names[rk]} {residue.chain}{int(residue.residue_number):4d}{residue.insertion_code or " "}   '
                             f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           {name[0]}  \n')
        if i == len(order)-1 or order[i+1].chain != residue.chain:
            lines.append("TER\n")
    (output / "input-backbone.pdb").write_text("".join(lines) + "END\n")
    (output / "source.pdb").write_bytes((sample / "source.pdb").read_bytes())
    (output / "request.json").write_bytes((sample / "request.json").read_bytes())
    payload = {"operation": "pack-only", "gaps": [], "residues": [asdict(r) for r in order],
               "fixed": [{"identity": list(a), "xyz": row["xyz"]} for a, row in source.items()],
               "generated": [list(key(a)) for a in request.generated_atoms],
               "candidates": 1, "seed": 23}
    save(output / "payload.json", payload)
    worker = Path("scripts/hybrid_gap_rosetta.py").resolve()
    files = [sample / name for name in ("plan.json", "request.json", "source.pdb", "native-endpoint.npz", "sampling.json", "independent-readback.json")]
    files += [Path(__file__).resolve(), worker, Path("scripts/audit_framedipt_research.py").resolve(),
              Path("scripts/hybrid_gap_reconstruction.py").resolve(), Path(args.rosetta_python).resolve(),
              output / "input-backbone.pdb", output / "payload.json"]
    save(output / "plan.json", {"protocol": "private-framedipt-backbone-plus-rosetta-sides-v1",
         "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "seed": 23,
         "new_diffusion_calls": 0, "fixed_heavy_atoms": len(source), "added_heavy_atoms": len(generated),
         "frame_derived_backbone_CB": [list(a) for a in frame_derived],
         "rosetta_constructed_side_atoms": [list(a) for a in sorted(generated-set(frame_derived))],
         "native_atomic_preservation_passed": False, "explicit_source_record_restoration": True,
         "backbone_closure": False, "heavy_minimization": False, "sequence_design": False,
         "public_qualification": False, "timeout_seconds": args.timeout})
    freeze_code(output / "plan.json")
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "PYTHONNOUSERSITE": "1"}
    execution = launch_worker([str(args.rosetta_python.absolute()), str(worker), str(output)], output, env, args.timeout)
    verify_bindings(json.loads((output / "plan.json").read_text())["source_code_and_input_sha256"])
    if execution["returncode"] or execution["guard_stop_reason"] is not None:
        raise ValueError("private sidechain worker did not complete")
    native = json.loads((output / "candidate-00/native.json").read_text())
    coordinates, maximum = audit_native(payload, native)
    text = materialize_candidate_pdb(source_text, request.generated_atoms,
        [coordinates[key(a)] for a in request.generated_atoms],
        tuple(names[(a.chain, a.residue_number, a.insertion_code)] for a in request.generated_atoms), request)
    candidate = output / "candidate.pdb"
    candidate.write_text(text)
    actual = atoms(candidate)
    if set(actual) != set(source) | generated or any(actual[a]["line"] != source[a]["line"] for a in source):
        raise ValueError("source records/inventory changed on materialization")
    assessment = assessment_adapter(request, source_text)(text)
    save(output / "assessment.json", asdict(assessment))
    findings = diagnose(candidate, output / "diagnose.json")
    report = {"source_records_exact": True, "rosetta_native_fixed_max_A": maximum,
              "native_diffusion_atomic_preservation": False, "geometry_passed": assessment.geometry_passed,
              "diagnose": diagnostic_counts(findings), "candidate_sha256": digest(candidate),
              "frame_derived_atom_count": len(frame_derived), "rosetta_constructed_atom_count": len(generated)-len(frame_derived),
              "public_output_published": False, "downstream_preparation_tested": False}
    save(output / "summary.json", report)
    print(json.dumps(report))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rosetta-python", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600, choices=range(1, 1801))
    args = parser.parse_args()
    try:
        run(args)
    except Exception as exc:
        import traceback
        if getattr(args, "owned_output", False):
            save(args.output / "failure.json", {"error": str(exc), "traceback": traceback.format_exc()})
        raise
