"""Research-only GFN-FF context experiment with exactly fixed surroundings.

Uses xTB's C single-point API and SciPy L-BFGS-B on selected coordinates. No
production adapter changes, peptide cuts, solvent boxes or accepted PDB outputs.
All energies are Hartree; optimizer coordinates/gradients use Angstrom.
"""
from __future__ import annotations

import argparse
import ctypes as ct
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from audit_xtb_controls import geometry_snapshot, pocket_snapshot, read_xyz
from audit_xtb_geometry_deltas import geometry_deltas

ROOT = Path(__file__).resolve().parents[1]
BOHR = 0.529177210903
KCAL_PER_HARTREE = 627.5094740631


class GFNFF:
    """Minimal typed C ABI binding; each process owns one calculator."""

    def __init__(self, library, numbers, coordinates, log):
        self.lib = ct.CDLL(str(library))
        pointer, double, integer = ct.c_void_p, ct.POINTER(ct.c_double), ct.POINTER(ct.c_int)
        signatures = {
            "xtb_newEnvironment": ([], pointer), "xtb_newCalculator": ([], pointer),
            "xtb_newResults": ([], pointer), "xtb_checkEnvironment": ([pointer], ct.c_int),
            "xtb_getError": ([pointer, ct.c_char_p, integer], None),
            "xtb_setOutput": ([pointer, ct.c_char_p], None),
            "xtb_setVerbosity": ([pointer, ct.c_int], None),
            "xtb_newMolecule": ([pointer, integer, integer, double, double, integer, double, pointer], pointer),
            "xtb_updateMolecule": ([pointer, pointer, double, double], None),
            "xtb_loadGFNFF": ([pointer, pointer, pointer, ct.c_char_p], None),
            "xtb_singlepoint": ([pointer, pointer, pointer, pointer], None),
            "xtb_getEnergy": ([pointer, pointer, double], None),
            "xtb_getGradient": ([pointer, pointer, double], None),
            **{f"xtb_del{x}": ([ct.POINTER(pointer)], None)
               for x in ["Environment", "Calculator", "Results", "Molecule"]},
        }
        for name, (args, result) in signatures.items():
            function = getattr(self.lib, name)
            function.argtypes, function.restype = args, result
        self.env = pointer(self.lib.xtb_newEnvironment())
        self.lib.xtb_setOutput(self.env, os.fsencode(log))
        self.lib.xtb_setVerbosity(self.env, 1)
        self.numbers = np.ascontiguousarray(numbers, dtype=np.int32)
        self.count = ct.c_int(len(numbers))
        self.coords = np.ascontiguousarray(coordinates / BOHR, dtype=np.float64)
        self.mol = pointer(self.lib.xtb_newMolecule(
            self.env, ct.byref(self.count), self.numbers.ctypes.data_as(integer),
            self.coords.ctypes.data_as(double), None, None, None, None))
        self.check()
        self.calc = pointer(self.lib.xtb_newCalculator())
        self.result = pointer(self.lib.xtb_newResults())
        self.lib.xtb_loadGFNFF(self.env, self.mol, self.calc, None)
        self.check()
        self.numeric_libraries = []
        if sys.platform == "linux":
            paths = {line.split()[-1] for line in Path("/proc/self/maps").read_text().splitlines()
                     if any(name in line for name in ["blas", "lapack", "gfortran", "mkl", "blis"])}
            self.numeric_libraries = [{"path": path, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
                                      for path in sorted(paths) if Path(path).is_file()]
        self.resolved_blas_symbols = {}
        if sys.platform == "linux":
            class DlInfo(ct.Structure):
                _fields_ = [("filename", ct.c_char_p), ("base", ct.c_void_p),
                            ("symbol", ct.c_char_p), ("address", ct.c_void_p)]

            dladdr = ct.CDLL(None).dladdr
            dladdr.argtypes, dladdr.restype = [ct.c_void_p, ct.POINTER(DlInfo)], ct.c_int
            for name in ["dgemm_", "dsyev_", "dsytrf_", "dgetrf_", "dpotrf_"]:
                info = DlInfo()
                function = getattr(self.lib, name, None)
                if function is not None and dladdr(ct.cast(function, ct.c_void_p), ct.byref(info)):
                    self.resolved_blas_symbols[name] = os.fsdecode(info.filename)

    def check(self):
        if self.lib.xtb_checkEnvironment(self.env):
            buffer = ct.create_string_buffer(8192)
            count = ct.c_int(len(buffer))
            self.lib.xtb_getError(self.env, buffer, ct.byref(count))
            raise RuntimeError(buffer.value.decode(errors="replace"))

    def evaluate(self, coordinates):
        double = ct.POINTER(ct.c_double)
        self.coords[:] = coordinates / BOHR
        self.lib.xtb_updateMolecule(self.env, self.mol, self.coords.ctypes.data_as(double), None)
        self.check()
        self.lib.xtb_singlepoint(self.env, self.mol, self.calc, self.result)
        self.check()
        energy, gradient = ct.c_double(), np.empty_like(self.coords)
        self.lib.xtb_getEnergy(self.env, self.result, ct.byref(energy))
        self.lib.xtb_getGradient(self.env, self.result, gradient.ctypes.data_as(double))
        self.check()
        if not np.isfinite(energy.value) or not np.isfinite(gradient).all():
            raise ValueError("Nonfinite xTB energy/gradient")
        return energy.value, gradient / BOHR

    def close(self):
        for name, handle in [("Results", self.result), ("Calculator", self.calc),
                             ("Molecule", self.mol), ("Environment", self.env)]:
            getattr(self.lib, f"xtb_del{name}")(ct.byref(handle))


def restrained_energy(energy, gradient, coordinates, reference, heavy, k):
    delta = coordinates[heavy] - reference[heavy]
    gradient = gradient.copy()
    gradient[heavy] += k * delta / KCAL_PER_HARTREE
    return energy + .5 * k * float(np.sum(delta**2)) / KCAL_PER_HARTREE, gradient


def validate_api(library, binary, output):
    import re

    from openmm.app import Element

    source = ROOT / ".artifacts/probe-xtb/linux-fresh/fixed-controls/ethanol-default/input.xyz"
    elements, coordinates = read_xyz(source)
    (output / "input.xyz").write_bytes(source.read_bytes())
    calc = GFNFF(library, [Element.getBySymbol(e).atomic_number for e in elements], coordinates,
                 output / "api.log")
    try:
        energy, gradient = calc.evaluate(coordinates)
        comparisons = []
        for atom, axis in [(0, 0), (0, 2), (2, 1), (3, 0), (5, 1), (8, 2)]:
            plus, minus = coordinates.copy(), coordinates.copy()
            plus[atom, axis] += 1e-5
            minus[atom, axis] -= 1e-5
            finite_difference = (calc.evaluate(plus)[0] - calc.evaluate(minus)[0]) / 2e-5
            comparisons.append({"atom": atom, "axis": axis, "analytic": float(gradient[atom, axis]),
                                "finite_difference": finite_difference,
                                "absolute_error": abs(finite_difference - gradient[atom, axis])})
    finally:
        calc.close()
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    run = subprocess.run([str(binary), "input.xyz", "--gfnff", "--sp", "--norestart"],
                         cwd=output, env=env, capture_output=True, timeout=60)
    (output / "native.stdout").write_bytes(run.stdout)
    (output / "native.stderr").write_bytes(run.stderr)
    matches = re.findall(r"TOTAL ENERGY\s+([-+\d.Ee]+)", run.stdout.decode(errors="replace"))
    cli_energy = float(matches[-1]) if matches else None
    passed = bool(run.returncode == 0 and cli_energy is not None and abs(cli_energy - energy) < 1e-7
                  and max(r["absolute_error"] for r in comparisons) < 1e-5)
    report = {"passed": passed, "api_energy_hartree": energy, "native_energy_hartree": cli_energy,
              "returncode": run.returncode, "gradient_checks": comparisons,
              "numeric_libraries": calc.numeric_libraries,
              "resolved_blas_symbols": calc.resolved_blas_symbols,
              "library_sha256": hashlib.sha256(library.read_bytes()).hexdigest()}
    (output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    if not passed:
        raise RuntimeError("API energy/gradient validation failed")


def load_case(path, glycam):
    from openmm.app import PDBFile
    from openmm.unit import nanometer

    from dvbfixer.ffutils import build_glycam_system
    from dvbfixer.minimize.refine import _extract_heterogen_subsystem

    pdb = PDBFile(str(path))
    if glycam:
        build_glycam_system(["amber14-all.xml", "amber14/GLYCAM_06j-1.xml", "amber14/tip3pfb.xml"],
                            pdb.topology, pdb.positions)
    subtop, subpos, mapping, frozen = _extract_heterogen_subsystem(pdb.topology, pdb.positions)
    before = np.asarray(pdb.positions.value_in_unit(nanometer)) * 10
    subcoords = np.asarray(subpos.value_in_unit(nanometer)) * 10
    mobile = sorted(full for full, sub in mapping.items() if sub not in frozen)
    return pdb, subtop, before, subcoords, mapping, mobile


def finding_sets(topology, coordinates, mobile):
    from openmm.unit import nanometer

    from dvbfixer.diagnose import chemistry, steric

    def key(finding):
        atoms = (finding.extra or {}).get("atoms", [])
        return (finding.category, str(finding.severity), tuple(sorted(
            (str(a.get("chain", "")), str(a.get("resid", "")), str(a.get("icode", "")),
             str(a.get("atom", ""))) for a in atoms)),
            finding.chain, finding.resid, finding.atom)

    positions = coordinates / 10 * nanometer
    clashes = steric.clashes_python(topology, positions, clash_warn_a=.4, clash_error_a=.5)
    chemical = chemistry.run_all(topology, positions)
    identities = {(a.residue.chain.id, a.residue.id, a.residue.insertionCode.strip(), a.name)
                  for a in topology.atoms() if a.index in mobile}
    errors = set()
    for finding in chemical:
        atoms = (finding.extra or {}).get("atoms", [])
        touches = any((str(a.get("chain", "")), str(a.get("resid", "")),
                       str(a.get("icode", "")), str(a.get("atom", ""))) in identities for a in atoms)
        if str(finding.severity) == "ERROR" and touches:
            errors.add(key(finding))
    return {key(f) for f in clashes}, errors


def worker(args):
    from scipy.optimize import minimize

    worker_source = Path(__file__).read_bytes()
    (args.output / "worker-script.py").write_bytes(worker_source)
    pdb, subtop, before, subcoords, mapping, mobile = load_case(args.baseline, args.glycam)
    full_atoms = list(pdb.topology.atoms())
    if args.scope == "isolated":
        reference, atoms = subcoords, list(subtop.atoms())
        moving = np.array(sorted(mapping[i] for i in mobile), dtype=int)
    else:
        reference, atoms = before, full_atoms
        moving = np.array(mobile, dtype=int)
    fixed = np.array(sorted(set(range(len(reference))) - set(moving)), dtype=int)
    heavy = np.array([i for i in moving if atoms[i].element.symbol != "H"], dtype=int)
    report = {"baseline_sha256": hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
              "threads": {key: os.environ.get(key) for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]},
              "scope": args.scope, "restraint_k_kcal_mol_angstrom2": args.k,
              "atoms": len(atoms), "mobile_atoms": len(moving), "fixed_atoms": len(fixed),
              "input_atom_identities": [(a.residue.chain.id, a.residue.id, a.residue.insertionCode.strip(),
                                         a.residue.name, a.name, a.element.symbol) for a in atoms],
              "full_to_sub_indices": mapping, "mobile_full_indices": mobile,
              "optimizer": "SciPy L-BFGS-B on mobile coordinates only", "gtol_hartree_per_angstrom": 1e-4,
              "total_charge_parameter": "xTB API default zero; no formal-charge reconstruction",
              "maxiter": args.maxiter, "history": []}
    work = args.output
    np.save(work / "reference.npy", reference)
    report_path = work / "report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    calc = GFNFF(args.library, [a.element.atomic_number for a in atoms], reference, work / "native.log")
    report["numeric_libraries"] = calc.numeric_libraries
    report["resolved_blas_symbols"] = getattr(calc, "resolved_blas_symbols", {})
    report["script_sha256"] = hashlib.sha256(worker_source).hexdigest()
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    first_energy = None

    def evaluate(vector):
        nonlocal first_energy
        coordinates = reference.copy()
        coordinates[moving] = vector.reshape(-1, 3)
        energy, gradient = calc.evaluate(coordinates)
        if first_energy is None:
            first_energy = energy
        value, gradient = restrained_energy(energy, gradient, coordinates, reference, heavy, args.k)
        row = {"evaluation": len(report["history"]), "native_energy_hartree": energy,
               "objective_hartree": value, "mobile_gradient_max_hartree_angstrom": float(np.abs(gradient[moving]).max())}
        report["history"].append(row)
        np.save(work / "last-evaluation.npy", coordinates)
        with (work / "evaluations.jsonl").open("a") as handle:
            handle.write(json.dumps(row) + "\n")
        return value, gradient[moving].reshape(-1)

    try:
        start = time.monotonic()
        def checkpoint(vector):
            coordinates = reference.copy()
            coordinates[moving] = vector.reshape(-1, 3)
            np.save(work / "last-optimizer-iterate.npy", coordinates)

        result = minimize(evaluate, reference[moving].reshape(-1), jac=True, method="L-BFGS-B",
                          callback=checkpoint,
                          options={"maxiter": args.maxiter, "maxfun": 3000, "gtol": 1e-4, "ftol": 1e-12, "maxls": 30})
        after = reference.copy()
        after[moving] = result.x.reshape(-1, 3)
        final_energy, final_gradient = calc.evaluate(after)
        _, projected = restrained_energy(final_energy, final_gradient, after, reference, heavy, args.k)
        np.save(work / "raw-final.npy", after)
        full_after = before.copy()
        if args.scope == "isolated":
            for full, sub in mapping.items():
                full_after[full] = after[sub]
        else:
            full_after = after
        assert not len(fixed) or np.array_equal(after[fixed], reference[fixed])
        np.save(work / "full-final.npy", full_after)
        before_pairs, before_errors = finding_sets(pdb.topology, before, set(mobile))
        after_pairs, after_errors = finding_sets(pdb.topology, full_after, set(mobile))
        delta = geometry_deltas(pdb.topology, before, full_after)
        pockets = pocket_snapshot(pdb.topology, before, full_after)
        gradient_max = float(np.abs(projected[moving]).max())
        converged = bool(result.success and gradient_max <= 1e-4)
        gates = {"converged": converged, "frozen_drift": True,
                 "bond_delta": not delta["bonds_changed_over_0_1_angstrom"],
                 "angle_delta": not delta["angles_changed_over_10_degrees"],
                 "tetrahedral_signs": not delta["potential_stereovolume_flips"],
                 "no_new_mobile_chemical_errors": not after_errors - before_errors,
                 "clash_count": len(after_pairs) <= len(before_pairs),
                 "pocket_contacts": all(len(p["retained_contacts"]) >= .8 * len(p["original_contacts"]) for p in pockets)}
        report.update({"seconds": time.monotonic() - start, "scipy_success": bool(result.success),
                       "optimizer_message": str(result.message), "iterations": int(result.nit),
                       "converged": converged, "projected_gradient_max_hartree_angstrom": gradient_max,
                       "initial_native_energy_hartree": first_energy, "final_native_energy_hartree": final_energy,
                       "raw_frozen_drift_max_angstrom": 0., "gates": gates, "screens_pass": all(gates.values()),
                       "geometry_before": geometry_snapshot(pdb.topology, before),
                       "geometry_after": geometry_snapshot(pdb.topology, full_after), "geometry_deltas": delta,
                       "clash_pairs_before": len(before_pairs), "clash_pairs_after": len(after_pairs),
                       "new_clash_pairs": sorted(after_pairs - before_pairs), "lost_clash_pairs": sorted(before_pairs - after_pairs),
                       "new_mobile_chemical_errors": sorted(after_errors - before_errors), "pockets": pockets})
        report_path.write_text(json.dumps(report, indent=2) + "\n")
    finally:
        calc.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--scope", choices=["isolated", "context"], default="context")
    parser.add_argument("--k", type=float, default=0.)
    parser.add_argument("--glycam", action="store_true")
    parser.add_argument("--maxiter", type=int, default=500)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--threads", type=int, default=1, help="Labelled numerical-thread sensitivity; default one")
    parser.add_argument("--single-arm", action="store_true", help="Bounded supplemental replay of --scope / --k only")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--validate-api", action="store_true")
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    if min(args.maxiter, args.timeout, args.threads) <= 0 or not np.isfinite(args.k) or args.k < 0:
        parser.error("Iterations, timeout and threads must be positive; k must be finite and nonnegative")
    args.library, args.output = args.library.resolve(), args.output.resolve()
    if args.worker:
        if args.validate_api:
            validate_api(args.library, args.binary.resolve(), args.output)
        else:
            worker(args)
        return
    args.output.mkdir(parents=True, exist_ok=False)
    if args.validate_api:
        if args.binary is None:
            parser.error("--validate-api requires --binary")
        command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--validate-api",
                   "--library", str(args.library), "--binary", str(args.binary.resolve()), "--output", str(args.output)]
        run = subprocess.run(command, cwd=args.output, capture_output=True, timeout=120)
        (args.output / "api.stdout").write_bytes(run.stdout)
        (args.output / "api.stderr").write_bytes(run.stderr)
        if run.returncode:
            raise RuntimeError("API validation worker failed; inspect api.stderr")
        return
    (args.output / "audit-script.py").write_bytes(Path(__file__).read_bytes())
    report = {"platform": platform.platform(), "python": sys.version,
              "library_sha256": hashlib.sha256(args.library.read_bytes()).hexdigest(),
              "invocation": sys.argv, "jobs": [], "acceptance_note": "Comparative screens, not production qualification"}
    cases = [("hinge", ROOT / ".artifacts/probe-xtb/linux-fresh/baselines/hinge-openmm.pdb", True),
             ("plm", ROOT / ".artifacts/probe-xtb/linux-fresh/baselines/7x35-openmm.pdb", False),
             ("1vcu", ROOT / ".artifacts/probe-xtb/linux-1vcu/openmm-fresh-cache.pdb", False)]
    if args.baseline:
        cases = [(args.baseline.stem, args.baseline.resolve(), args.glycam)]

    def run_job(job):
        label, path, glycam, scope, k = job
        name = f"{label}-{scope}-k{k:g}"
        work = args.output / name
        work.mkdir()
        command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--library", str(args.library),
                   "--output", str(work), "--baseline", str(path), "--scope", scope, "--k", str(k),
                   "--maxiter", str(args.maxiter), *(["--glycam"] if glycam else [])]
        row = {"name": name, "command": command}
        env = {**os.environ, "OMP_NUM_THREADS": str(args.threads), "OPENBLAS_NUM_THREADS": str(args.threads),
               "MKL_NUM_THREADS": str(args.threads)}
        start = time.monotonic()
        try:
            run = subprocess.run(command, cwd=work, env=env, capture_output=True, timeout=args.timeout)
            row["returncode"] = run.returncode
            (work / "stdout.raw").write_bytes(run.stdout)
            (work / "stderr.raw").write_bytes(run.stderr)
        except subprocess.TimeoutExpired as exc:
            row["timeout"] = True
            (work / "stdout.raw").write_bytes(exc.stdout or b"")
            (work / "stderr.raw").write_bytes(exc.stderr or b"")
        row["seconds"] = time.monotonic() - start
        print(name, row, flush=True)
        return row

    arms = [(args.scope, args.k)] if args.single_arm else [("isolated", 0.), ("context", 0.), ("context", 5.)]
    jobs = [(label, path, glycam, scope, k) for label, path, glycam in cases for scope, k in arms]
    with ThreadPoolExecutor(max_workers=2) as pool:
        for row in pool.map(run_job, jobs):
            report["jobs"].append(row)
            (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest = {str(p.relative_to(args.output)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(args.output.rglob("*")) if p.is_file()}
    (args.output / "manifest.sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
