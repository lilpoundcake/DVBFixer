"""Read-only localization of saved 1NON refinement failures; no FF/optimizer calls."""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from openmm import XmlSerializer
from scripts.anchor_calibration_grid import read_temporary
from scripts.anchor_constrained_coordinates import ReducedCoordinates
from scripts.audit_anchor_constrained_control import independent_projection
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import verify_bindings

from dvbfixer.diagnose.geometry import stereocentre_geometry
from dvbfixer.model.diffusion.contract import DiffusionRequest


def key(atom):
    return atom.chain, atom.residue_number, atom.insertion_code, atom.atom_name


def finding_key(atom):
    return atom["chain"], str(atom["resid"]), atom.get("icode", ""), atom["atom"]


def classify_clashes(findings, generated, fixed):
    groups, pairs, residues = Counter(), [], Counter()
    for finding in findings:
        if finding["severity"] != "ERROR" or finding["category"] != "clash":
            continue
        a, b = [finding_key(a) for a in finding["details"]["atoms"]]
        if a not in generated | fixed or b not in generated | fixed:
            raise ValueError("clash endpoint outside frozen ownership")
        if a in fixed and b in fixed:
            group = "fixed-fixed"
        elif a in fixed or b in fixed:
            group = "generated-fixed"
        else:
            group = "generated-generated-same-chain" if a[0] == b[0] else "generated-generated-cross-chain"
        groups[group] += 1
        for residue in {a[:3], b[:3]}:
            residues[residue] += 1
        pairs.append({"atoms": [a, b], "ownership": group,
                      "overlap_angstrom": finding["details"]["overlap_angstrom"]})
    return {"counts": dict(groups), "pairs": sorted(pairs, key=lambda p: -p["overlap_angstrom"]),
            "residue_contact_counts": [{"residue": r, "contacts": n} for r, n in residues.most_common()]}


def geometry(points):
    centres = {}
    for chain in "ABCD":
        for resid in range(73, 91):
            keys = [(chain, str(resid), "", name) for name in ("CA", "N", "C", "CB")]
            if all(k in points for k in keys):
                measured = stereocentre_geometry(*(points[k] for k in keys))
                centres[f"{chain}{resid}"] = None if measured is None else measured.out_of_plane_degrees
    return {"ca_out_of_plane_degrees": centres,
            "near_planar_ca": {k: v for k, v in centres.items() if v is not None and abs(v) <= 10},
            "junction_cn_angstrom": {c: [float(np.linalg.norm(points[(c, l, "", "C")]-points[(c, r, "", "N")]))
                                         for l, r in (("72", "73"), ("90", "91"))] for c in "ABCD"}}


def callback_summary(rows):
    if not rows:
        raise ValueError("missing saved callbacks")
    result = {"count": len(rows), "stationarity_pass_callbacks": sum(r["stationarity_passed"] for r in rows)}
    for n in (100, 250):
        tail = rows[-n:]
        rms = [r["projected_cartesian_dof_rms_kj_mol_nm"] for r in tail]
        result[f"last_{n}"] = {"count": len(tail), "rms_min": min(rms), "rms_median": float(np.median(rms)),
                              "rms_max": max(rms), "energy_change_kj_mol": tail[-1]["energy_kj_mol"]-tail[0]["energy_kj_mol"]}
    return result


def run(root, output):
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    job = plan["jobs"][0]
    case, prepared = Path(job["case"]), Path(job["snapshot"])
    work = prepared.parent
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    generated, fixed = set(map(key, request.generated_atoms)), set(map(key, request.fixed_atoms))
    paths = [plan_path, Path(__file__), case / "request.json", case / "proposal.pdb", work / "system.xml",
             work / "full-cli.json", prepared / "temporary-axis.json", prepared / "initial-state.npz",
             prepared / "base-system.xml"]
    paths += sorted(work.glob("phase-*/*.npz"))+sorted(work.glob("phase-*/completed.json"))+sorted(work.glob("phase-*/callbacks.jsonl"))
    bindings = {str(p.resolve()): digest(p) for p in paths}
    output.mkdir(parents=True, exist_ok=False)
    topology, _, initial = read_temporary(prepared)
    system = XmlSerializer.deserialize((work / "system.xml").read_text())
    chart = ReducedCoordinates(topology, system, initial)
    identities = [(a.residue.chain.id, a.residue.id, a.residue.insertionCode.strip(), a.name) for a in topology.atoms()]
    heavy = np.array([a.element.symbol != "H" for a in topology.atoms()])
    raw = {k: np.array(v["xyz"]) for k, v in atoms(case / "proposal.pdb").items()}
    frames = [{"stage": "raw", **geometry(raw)}]
    for path in sorted(work.glob("phase-*/*.npz")):
        with np.load(path, allow_pickle=False) as data:
            points = data["positions_nm"].copy()
        if not np.array_equal(points[~chart.mobile], initial[~chart.mobile]):
            raise ValueError("saved state changed fixed atoms")
        frame = {"stage": f"{path.parent.name}/{path.stem}",
                 **geometry({k: points[i]*10 for i, k in enumerate(identities) if heavy[i]})}
        frames.append(frame)
    phases = []
    for directory in sorted(work.glob("phase-*")):
        with np.load(directory / "returned-state.npz", allow_pickle=False) as data:
            points, forces, q = (data[k].copy() for k in ("positions_nm", "raw_forces", "reduced_coordinates"))
        if not np.array_equal(chart.positions(q), points):
            raise ValueError("saved chart coordinates do not reconstruct the returned state")
        projected = chart.project_forces(points, forces)
        rms, residual = independent_projection(system, chart.bonds, points, forces)
        contributions = defaultdict(float)
        residue_forces = defaultdict(float)
        parent_by_h = dict(zip(chart.h, chart.parent, strict=True))
        for i, a in enumerate(identities):
            if not chart.mobile[i]:
                continue
            group = "generated-heavy" if heavy[i] else ("H-on-fixed-heavy" if not chart.mobile[parent_by_h[i]] else "H-on-generated-heavy")
            square = float(np.sum(projected[i]**2))
            contributions[group] += square
            residue_forces[a[:3]] += square
        total = sum(contributions.values())
        uv = q[chart.heavy_size:].reshape(-1, 2)/(2*chart.length[:, None])
        contraction = 1/(1+np.sum(uv*uv, axis=1))
        completed = json.loads((directory / "completed.json").read_text())
        callbacks = [json.loads(line) for line in (directory / "callbacks.jsonl").read_text().splitlines()]
        phases.append({"phase": directory.name, "termination": completed["scipy_message"],
            "iterations": completed["iterations"], "independent_projected_dof_rms": rms,
            "independent_tangency_residual": residual,
            "force_squared_fraction": {k: v/total for k, v in contributions.items()},
            "top_force_residues": [{"residue": r, "squared_force_fraction": s/total}
                                   for r, s in sorted(residue_forces.items(), key=lambda p: -p[1])[:15]],
            "hydrogen_chart_contraction_min_median": [float(contraction.min()), float(np.median(contraction))],
            "hydrogen_chart_contraction_below_1e_minus3": int(np.sum(contraction < .001)),
            "callbacks": callback_summary(callbacks)})
    clashes = classify_clashes(json.loads((work / "full-cli.json").read_text())["findings"], generated, fixed)
    verify_bindings(bindings)
    save(output / "analysis.json", {"plan_sha256": digest(plan_path), "input_sha256": bindings,
         "new_optimizer_calls": 0, "new_diffusion_calls": 0, "new_FF_evaluations": 0,
         "final_clashes": clashes, "saved_geometry_frames": frames, "phases": phases,
         "limitations": ["Only saved checkpoints localize geometry onset; missing iterations were not replayed.",
                         "Force groups differ between phases; energies cannot be compared across phases.",
                         "No force-term ablation or new trajectory: budget versus basin causality remains unresolved."]})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
