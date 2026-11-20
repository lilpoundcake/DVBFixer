"""Private, reference-blind replay of saved diffusion backbones through GenKIC.

This is not a public backend. Internal canonical whole-residue gaps only; source
atoms never become optimization variables. PyRosetta runs in its own existing
environment. Failures and every candidate remain private evidence, not outputs.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_anchor_calibration_grid import diagnose, record_audit
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import freeze_code, verify_bindings
from scripts.run_mac_reconstruction_guarded import resource_stop_reason, stop_owned_process

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

GEOMETRY_CATEGORIES = frozenset({
    "clash", "chirality", "ca_planarity", "coincident_atoms", "bond_length",
    "bond_angle", "non_planar_amide", "ramachandran", "rotamer", "chain_break",
})


def diagnostic_counts(findings):
    errors = [f for f in findings if f["severity"] == "ERROR"]
    return {
        "geometry_errors": sum(f["category"] in GEOMETRY_CATEGORIES for f in errors),
        "other_errors": sum(f["category"] not in GEOMETRY_CATEGORIES for f in errors),
        "all_error_categories": sorted({f["category"] for f in errors}),
    }


def key(atom):
    return atom.chain, atom.residue_number, atom.insertion_code, atom.atom_name


def dihedral(a, b, c, d):
    """Signed radians; degenerate native backbones are refused, not guessed."""
    a, b, c, d = map(lambda p: np.asarray(p, dtype=float), (a, b, c, d))
    axis = c - b
    if np.linalg.norm(axis) < 1e-8:
        raise ValueError("degenerate backbone axis")
    axis /= np.linalg.norm(axis)
    v, w = a - b, d - c
    v -= np.dot(v, axis) * axis
    w -= np.dot(w, axis) * axis
    if min(np.linalg.norm(v), np.linalg.norm(w)) < 1e-8:
        raise ValueError("degenerate backbone plane")
    return float(np.arctan2(np.dot(np.cross(axis, v), w), np.dot(v, w)))


def extend(a, b, c, length, angle, torsion):
    """Internal-coordinate atom placement in Angstrom/radians."""
    a, b, c = map(lambda p: np.asarray(p, dtype=float), (a, b, c))
    u = c - b
    normal = np.cross(b - a, u)
    if min(np.linalg.norm(u), np.linalg.norm(normal)) < 1e-8:
        raise ValueError("degenerate construction frame")
    u /= np.linalg.norm(u)
    normal /= np.linalg.norm(normal)
    side = np.cross(normal, u)
    return c + length * (-math.cos(angle) * u + math.sin(angle) * (
        math.cos(torsion) * side + math.sin(torsion) * normal))


def build_backbone(left, torsions):
    """Build an open ideal chain; closure, not a distance penalty, follows later.

    The preceding observed carbonyl fixes the first generated N. All subsequent
    peptide bonds are trans. Proline cis sampling is intentionally not claimed.
    """
    rad = math.radians
    ca, carbon, oxygen = map(np.asarray, left)
    nitrogen = extend(oxygen, ca, carbon, 1.33, rad(116.2), math.pi)
    alpha = extend(ca, carbon, nitrogen, 1.46, rad(121.7), math.pi)
    previous_c = carbon
    result = []
    for phi, psi in torsions:
        carbon = extend(previous_c, nitrogen, alpha, 1.525, rad(111.2), phi)
        next_n = extend(nitrogen, alpha, carbon, 1.33, rad(116.2), psi)
        oxygen = extend(nitrogen, alpha, carbon, 1.23, rad(120.8), psi + math.pi)
        result.append(dict(zip(("N", "CA", "C", "O"), (nitrogen, alpha, carbon, oxygen))))
        next_ca = extend(alpha, carbon, next_n, 1.46, rad(121.7), math.pi)
        previous_c, nitrogen, alpha = carbon, next_n, next_ca
    return result


def build_backbone_from_right(right, torsions, right_phi):
    """Reverse growth fixes the right peptide plane before an INTERNAL cut closes.

    GenKIC perturbers cannot set a bond to an outside-loop residue. Both boundary
    bonds must therefore already have correct geometry in the open input chain.
    """
    nitrogen, alpha, carbon = map(np.asarray, right)
    result = []
    phi_next = right_phi
    for phi, psi in reversed(torsions):
        previous_c = extend(carbon, alpha, nitrogen, 1.33, math.radians(121.7), phi_next)
        previous_ca = extend(alpha, nitrogen, previous_c, 1.525, math.radians(116.2), math.pi)
        previous_n = extend(nitrogen, previous_c, previous_ca, 1.46, math.radians(111.2), psi)
        oxygen = extend(previous_n, previous_ca, previous_c, 1.23, math.radians(120.8), psi + math.pi)
        result.append(dict(zip(("N", "CA", "C", "O"), (previous_n, previous_ca, previous_c, oxygen))))
        nitrogen, alpha, carbon, phi_next = previous_n, previous_ca, previous_c, phi
    return list(reversed(result))


def backbone_jacobian(rows):
    """Cartesian derivative w.r.t. phi/psi, without finite-difference rebuilds."""
    points = np.array([[r[n] for n in ("N", "CA", "C", "O")] for r in rows])
    jac = np.zeros((len(rows), 4, 3, 2 * len(rows)))
    for i, row in enumerate(rows):
        for offset, (origin, endpoint) in enumerate(((row["N"], row["CA"]), (row["CA"], row["C"]))):
            axis = endpoint - origin
            axis /= np.linalg.norm(axis)
            column = 2 * i + offset
            jac[i + 1:, :, :, column] = np.cross(axis, points[i + 1:] - origin)
            if offset == 0:
                jac[i, 2:, :, column] = np.cross(axis, points[i, 2:] - origin)
            else:
                jac[i, 3, :, column] = np.cross(axis, points[i, 3] - origin)
    return jac


def fit_route(left, right, torsions, right_phi, raw_ca, max_nfev=200):
    """Inverse kinematics to a saved C-alpha route; never a reference structure.

    Endpoint penalties initialize the search only. Acceptance still requires
    subsequent exact GenKIC closure and independent full-atom geometry checks.
    No endpoint penalty is substituted for the final closure condition.
    """
    from scipy.optimize import least_squares

    initial = np.asarray(torsions, dtype=float).ravel()
    terminal = build_backbone_from_right(right, torsions, right_phi)[-1]
    target = np.array([terminal[n] for n in ("N", "CA", "C")])
    raw_ca = np.asarray(raw_ca)

    def residual(q):
        rows = build_backbone(left, q.reshape(-1, 2))
        ca = np.array([r["CA"] for r in rows])
        end = np.array([rows[-1][n] for n in ("N", "CA", "C")])
        return np.concatenate(((ca - raw_ca).ravel(), 5 * (end - target).ravel(),
                               .02 * np.sin(q - initial), .02 * (np.cos(q - initial) - 1)))

    def jacobian(q):
        jac = backbone_jacobian(build_backbone(left, q.reshape(-1, 2)))
        return np.vstack((jac[:, 1].reshape(-1, len(q)), 5 * jac[-1, :3].reshape(-1, len(q)),
                          np.diag(.02 * np.cos(q - initial)), np.diag(-.02 * np.sin(q - initial))))

    fit = least_squares(residual, initial, jac=jacobian, max_nfev=max_nfev, ftol=1e-8, xtol=1e-8, gtol=1e-8)
    rows = build_backbone(left, fit.x.reshape(-1, 2))
    report = {"success": bool(fit.success), "evaluations": fit.nfev, "cost": float(fit.cost),
              "endpoint_max_A": float(np.linalg.norm(np.array([rows[-1][n] for n in ("N", "CA", "C")]) - target, axis=1).max()),
              "ca_rms_to_diffusion_A": float(np.sqrt(np.mean(np.sum((np.array([r["CA"] for r in rows]) - raw_ca)**2, axis=1)))),
              "max_nfev": max_nfev, "observed_variables": 0, "reference_coordinates_read": False}
    return fit.x.reshape(-1, 2).tolist(), report


def validate_inputs(request, source, proposal):
    fixed, generated = set(map(key, request.fixed_atoms)), set(map(key, request.generated_atoms))
    if set(source) != fixed or set(proposal) != fixed | generated or fixed & generated:
        raise ValueError("source/proposal inventory differs from the frozen ownership")
    if any(source[a]["line"] != proposal[a]["line"] for a in fixed):
        raise ValueError("proposal changed an observed record")
    axis = target_axis(request)
    represented = {(r.chain, r.residue_number, r.insertion_code) for chain in axis.values() for r in chain}
    if any(a[:3] not in represented for a in fixed):
        raise ValueError("full protein context must be declared on the target axis; omitted chains are refused")
    for rows in (source, proposal):
        if any(r["line"].startswith("HETATM") or r["element"] in {"H", "D"} for r in rows.values()):
            raise ValueError("private hybrid requires protein-only heavy atoms")
    expected_generated = set()
    for target in request.target_sequences:
        for ordinal, (residue, letter) in enumerate(zip(axis[target.chain], target.sequence, strict=True)):
            expected = _CANONICAL_ONE_TO_THREE[letter]
            names = {a[3] for a in proposal if a[:3] == (residue.chain, residue.residue_number, residue.insertion_code)}
            if not set(CANONICAL_HEAVY_ATOMS[expected]) <= names or names - set(CANONICAL_HEAVY_ATOMS[expected]) - {"OXT"}:
                raise ValueError("incomplete or noncanonical target residue")
            if "OXT" in names and ordinal != len(target.sequence) - 1:
                raise ValueError("OXT is only allowed on a true target terminus")
            if any(r["resname"] != expected for a, r in proposal.items() if a[:3] == (residue.chain, residue.residue_number, residue.insertion_code)):
                raise ValueError("proposal sequence differs from request")
    for gap in request.gaps:
        if gap.left_anchor is None or gap.right_anchor is None or not 3 <= len(gap.generated_residues) <= 25:
            raise ValueError("initial hybrid supports internal gaps of 3..25 residues only")
        for residue in gap.generated_residues:
            owned = {a for a in proposal if a[:3] == (residue.chain, residue.residue_number, residue.insertion_code)}
            if not owned <= generated:
                raise ValueError("partial-residue ownership is not a whole-gap reconstruction")
            expected_generated.update(owned)
    if expected_generated != generated:
        raise ValueError("generated axis contains atoms outside the declared whole gaps")
    return axis


def prepare_payload(request, source, proposal):
    axis = validate_inputs(request, source, proposal)
    residue_order = [r for target in request.target_sequences for r in axis[target.chain]]
    indices = {(r.chain, r.residue_number, r.insertion_code): i + 1 for i, r in enumerate(residue_order)}
    coordinates, gaps = {}, []
    for gap in request.gaps:
        residues = [gap.left_anchor, *gap.generated_residues, gap.right_anchor]
        def point(r, name):
            return np.array(proposal[(r.chain, r.residue_number, r.insertion_code, name)]["xyz"])
        torsions = []
        for previous, residue, following in zip(residues, residues[1:], residues[2:]):
            torsions.append((dihedral(point(previous, "C"), point(residue, "N"), point(residue, "CA"), point(residue, "C")),
                             dihedral(point(residue, "N"), point(residue, "CA"), point(residue, "C"), point(following, "N"))))
        right_phi = dihedral(point(gap.generated_residues[-1], "C"), point(gap.right_anchor, "N"),
                             point(gap.right_anchor, "CA"), point(gap.right_anchor, "C"))
        left = [point(gap.left_anchor, n) for n in ("CA", "C", "O")]
        right = [point(gap.right_anchor, n) for n in ("N", "CA", "C")]
        torsions, route_fit = fit_route(left, right, torsions, right_phi,
                                       [point(r, "CA") for r in gap.generated_residues])
        backbone = build_backbone(left, torsions)
        reverse = build_backbone_from_right([point(gap.right_anchor, n) for n in ("N", "CA", "C")], torsions, right_phi)
        split = len(backbone) // 2
        backbone[split:] = reverse[split:]
        for residue, points in zip(gap.generated_residues, backbone, strict=True):
            coordinates.update({AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, n): p for n, p in points.items()})
        ids = [indices[(r.chain, r.residue_number, r.insertion_code)] for r in residues]
        if ids != list(range(ids[0], ids[-1] + 1)):
            raise ValueError("gap anchors are not consecutive on the request axis")
        gaps.append({"start": ids[1], "stop": ids[-2], "left": ids[0], "right": ids[-1], "cut": ids[split],
                     "route_fit": route_fit})
    text = "".join(r["line"] for r in proposal.values()) + "END\n"
    text = rewrite_generated_coordinates(text, coordinates)
    generated = set(map(key, request.generated_atoms))
    kept = []
    # Explicit ordinal order for the Rosetta adapter; original records stay in
    # source/proposal artifacts and only generated coordinates are exported later.
    lines = {}
    for line in text.splitlines(keepends=True):
        if line.startswith("ATOM  "):
            a = line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()
            lines[a] = line
    for residue in residue_order:
        rk = residue.chain, residue.residue_number, residue.insertion_code
        for a in proposal:
            if a[:3] == rk and (a not in generated or a[3] in {"N", "CA", "C", "O"}):
                kept.append(lines[a])
        if residue == axis[residue.chain][-1]:
            kept.append("TER\n")
    return "".join(kept) + "END\n", {
        "gaps": gaps, "residues": [asdict(r) for r in residue_order],
        "fixed": [{"identity": list(a), "xyz": r["xyz"]} for a, r in source.items()],
        "generated": [list(key(a)) for a in request.generated_atoms],
        "reference_coordinates_read": False,
    }


def audit_native(payload, native):
    expected = {tuple(r["identity"]): r["xyz"] for r in payload["fixed"]}
    actual = {tuple(r["identity"]): r["xyz"] for r in native["atoms"]}
    if len(actual) != len(native["atoms"]) or set(actual) != set(expected) | set(map(tuple, payload["generated"])):
        raise ValueError("native axis changed")
    error = max((math.dist(xyz, actual[a]) for a, xyz in expected.items()), default=0.)
    if error > 1e-7:
        raise ValueError(f"native observed motion before copying: {error} Angstrom")
    if not np.isfinite(list(actual.values())).all():
        raise ValueError("nonfinite native output")
    return actual, error


def launch_worker(command, output, env, timeout):
    """Reuse the Mac pressure-v2 limits, without a false MPS-execution claim."""
    import psutil

    start = time.monotonic()
    peak = 0
    critical_since = low_since = None
    reason = None
    telemetry = []
    with (output / "worker.log").open("xb") as stream:
        process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, env=env, start_new_session=True)
        try:
            while process.poll() is None:
                now = time.monotonic()
                try:
                    child = psutil.Process(process.pid)
                    rss = child.memory_info().rss + sum(p.memory_info().rss for p in child.children(recursive=True))
                    peak = max(peak, rss)
                except (psutil.NoSuchProcess, psutil.ZombieProcess):
                    pass
                available = psutil.virtual_memory().available
                status = subprocess.run(["sysctl", "-n", "kern.memorystatus_vm_pressure_level"],
                                        capture_output=True, text=True, timeout=5)
                if status.returncode or not status.stdout.strip().isdigit():
                    reason = "host-pressure-telemetry-unavailable"
                    pressure = -1
                else:
                    pressure = int(status.stdout.strip())
                    critical_since = (critical_since if critical_since is not None else now) if pressure >= 4 else None
                    low_since = (low_since if low_since is not None else now) if available < 1024**3 else None
                    reason = resource_stop_reason(elapsed=now-start, timeout=timeout, peak_rss=peak,
                        pressure=pressure, critical_seconds=0 if critical_since is None else now-critical_since,
                        low_headroom_seconds=0 if low_since is None else now-low_since)
                telemetry.append({"elapsed_seconds": now-start, "rss_peak_bytes": peak,
                                  "available_bytes": available, "pressure_level": pressure})
                if reason:
                    stop_owned_process(process)
                    break
                time.sleep(.5)
        except BaseException:
            stop_owned_process(process)
            raise
        returncode = process.wait()
    report = {"returncode": returncode, "guard_stop_reason": reason,
              "wall_seconds": time.monotonic()-start, "peak_aggregate_rss_bytes": peak, "telemetry": telemetry}
    save(output / "execution.json", report)
    return report


def run(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    args.owned_output = True
    paths = {"source": args.source.resolve(), "proposal": args.proposal.resolve(), "request": args.request.resolve()}
    request = DiffusionRequest.from_json(paths["request"].read_text())
    if digest(paths["source"]) != request.normalized_pdb.sha256:
        raise ValueError("source bytes differ from the frozen request")
    source, proposal = atoms(paths["source"]), atoms(paths["proposal"])
    validate_inputs(request, source, proposal)
    for name, path in paths.items():
        (output / (name + (".json" if name == "request" else ".pdb"))).write_bytes(path.read_bytes())
    worker = Path(__file__).with_name("hybrid_gap_rosetta.py").resolve()
    files = [*paths.values(), Path(__file__).resolve(), worker, Path(args.rosetta_python).resolve(),
             Path(__file__).with_name("audit_hybrid_gap_reconstruction.py").resolve(),
             Path("tests/test_hybrid_gap_reconstruction.py").resolve(),
             Path("scripts/mac_reconstruction_comparison.py").resolve(), Path("scripts/retained_evidence.py").resolve(),
             Path("scripts/audit_anchor_calibration_grid.py").resolve(),
             Path("scripts/audit_retained_full_context_controls.py").resolve(),
             *Path("src/dvbfixer/model/diffusion").glob("*.py")]
    files = [p.resolve() for p in files]
    save(output / "plan.json", {
        "protocol": "private-saved-diffusion-genkic-packing-v3", "stage": "route-fit-closure-carbonyls-and-packing-no-relaxation",
        "source_code_and_input_sha256": {str(p): digest(p) for p in files},
        "new_diffusion_calls": 0, "MODELLER_calls": 0, "observed_anchor_release": False,
        "reference_coordinates_read": False, "candidates": args.candidates, "closure_attempts_per_gap": args.attempts,
        "seed": args.seed, "timeout_seconds": args.timeout, "public_qualification": False,
        "ranking": "independent geometry pass first, then Rosetta full-atom score; no reference score",
        "diagnostic_geometry_categories": sorted(GEOMETRY_CATEGORIES),
        "other_diagnostic_errors": "retained and reported separately; no downstream-preparation readiness claimed",
        "perturbation_degrees": [0, 10, 20, 40], "selection_after_packing": True,
        "stop_after_first_successful_closure_batch": True,
        "route_fit": "CA weight1; endpoint weight5; periodic torsion prior0.02; analytic Jacobian; max_nfev200",
        "carbonyl_materialization": "generated O only, 1.23A/120.8deg in the closed peptide plane before joint packing",
        "preselector_limit": "native per-root packing precedes the explicit post-closure carbonyl rebuild",
        "resources": "serial CPU; threads1; aggregate RSS8GiB; pressure-v2 host guard; bounded wall time; no solvent",
    })
    freeze_code(output / "plan.json")
    prepare_start = time.monotonic()
    backbone, payload = prepare_payload(request, source, proposal)
    (output / "input-backbone.pdb").write_text(backbone)
    payload.update(candidates=args.candidates, attempts=args.attempts, seed=args.seed)
    save(output / "payload.json", payload)
    prepared = {str(p): digest(p) for p in (output / "input-backbone.pdb", output / "payload.json")}
    save(output / "prepared-bindings.json", {"plan_sha256": digest(output / "plan.json"),
         "files": prepared, "wall_seconds": time.monotonic()-prepare_start})
    assess = assessment_adapter(request, paths["source"].read_text())
    baseline = assess(paths["proposal"].read_text())
    save(output / "baseline-assessment.json", {**asdict(baseline), "geometry_passed": baseline.geometry_passed})
    baseline_counts = diagnostic_counts(diagnose(paths["proposal"], output / "baseline-diagnose.json"))
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "PYTHONNOUSERSITE": "1"}
    command = [str(Path(args.rosetta_python).resolve()), str(worker), str(output)]
    save(output / "command.json", command)
    execution = launch_worker(command, output, env, args.timeout)
    returncode = execution["returncode"]
    timed_out = execution["guard_stop_reason"] == "wall-time-limit"
    verify_bindings(json.loads((output / "plan.json").read_text())["source_code_and_input_sha256"])
    verify_bindings(prepared)
    results = []
    for path in sorted(output.glob("candidate-*/native.json")):
        native = json.loads(path.read_text())
        points, error = audit_native(payload, native)
        updates = {a: np.array(points[key(a)]) for a in request.generated_atoms}
        text = rewrite_generated_coordinates(paths["proposal"].read_text(), updates)
        record_audit(paths["proposal"].read_text(), text, set(map(key, request.generated_atoms)))
        candidate = path.parent / "candidate.pdb"
        candidate.write_text(text)
        assessment = assess(text)
        save(path.parent / "assessment.json", asdict(assessment))
        findings = diagnose(candidate, path.parent / "diagnose.json")
        movement = np.array([points[key(a)] for a in request.generated_atoms]) - np.array([proposal[key(a)]["xyz"] for a in request.generated_atoms])
        results.append({"candidate": str(candidate.relative_to(output)), "sha256": digest(candidate),
                        "native_fixed_max_A": error, "geometry_passed": assessment.geometry_passed,
                        "admission_failures": list(assessment.admission_failures),
                        "diagnose": diagnostic_counts(findings),
                        "generated_heavy_rms_from_diffusion_A": float(np.sqrt(np.mean(np.sum(movement**2, axis=1)))),
                        "score": native["score"], "closure": native["closure"]})
    passing = [r for r in results if r["geometry_passed"] and not r["admission_failures"] and not r["diagnose"]["geometry_errors"]]
    completed = (returncode == 0 and execution["guard_stop_reason"] is None
                 and (output / "worker-completed.json").is_file())
    winner = min(passing, key=lambda r: r["score"]) if passing and completed else None
    save(output / "summary.json", {"returncode": returncode, "timed_out": timed_out,
        "guard_stop_reason": execution["guard_stop_reason"],
        "baseline_geometry_passed": baseline.geometry_passed, "baseline_diagnose": baseline_counts,
        "candidates": results, "winner": winner, "public_output_published": False,
        "downstream_preparation_tested": False, "topology_qualification": False,
        "status": "private-geometry-pass" if winner else "no-accepted-candidate"})
    print(json.dumps({"output": str(output), "completed": len(results), "passing": len(passing), "returncode": returncode}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "proposal", "request", "output", "rosetta-python"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--candidates", type=int, default=4, choices=range(1, 9))
    parser.add_argument("--attempts", type=int, default=32, choices=range(1, 129))
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--timeout", type=int, default=900, choices=range(1, 1801))
    args = parser.parse_args()
    try:
        run(args)
    except Exception as exc:
        import traceback
        if getattr(args, "owned_output", False) and not (args.output / "failure.json").exists():
            save(args.output / "failure.json", {"error": str(exc), "traceback": traceback.format_exc()})
        raise


if __name__ == "__main__":
    main()
