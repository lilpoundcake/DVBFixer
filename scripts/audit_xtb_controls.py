"""Run bounded xTB controls, retaining raw output before any coordinate merge.

Not a protein refinement benchmark. Requires RDKit and the existing xTB binary.
Output directories must not exist, so evidence cannot be silently overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text().splitlines()
    count = int(lines[0])
    rows = [line.split() for line in lines[2:]]
    if len(rows) != count or any(len(row) != 4 for row in rows):
        raise ValueError("Invalid XYZ row count/shape")
    coords = np.array([[float(x) for x in row[1:]] for row in rows])
    if not np.isfinite(coords).all():
        raise ValueError("Nonfinite XYZ coordinates")
    return [row[0] for row in rows], coords


def drift_metrics(before: np.ndarray, after: np.ndarray) -> dict:
    delta = after - before
    # Separate raw drift, translation and residual deformation after rigid fit.
    p, q = before - before.mean(0), after - after.mean(0)
    u, _, vt = np.linalg.svd(p.T @ q)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    residual = q - p @ rotation
    return {
        "max_angstrom": float(np.linalg.norm(delta, axis=1).max()),
        "rms_angstrom": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
        "centroid_shift_angstrom": float(np.linalg.norm(delta.mean(0))),
        "rigid_fit_rms_angstrom": float(np.sqrt(np.mean(np.sum(residual**2, axis=1)))),
    }


def last_complete_xyz_frame(path: Path) -> str | None:
    """Recover only complete trajectory frames; never imply final convergence."""
    lines = path.read_text().splitlines()
    index, last = 0, None
    while index < len(lines):
        try:
            count = int(lines[index])
        except ValueError:
            break
        end = index + count + 2
        if count <= 0 or end > len(lines):
            break
        last = "\n".join(lines[index:end]) + "\n"
        index = end
    return last


def geometry_snapshot(topology, coordinates_angstrom):
    from collections import Counter

    from openmm.unit import Quantity, nanometer

    from dvbfixer.diagnose import chemistry, steric

    positions = Quantity(coordinates_angstrom / 10, nanometer)
    findings = chemistry.run_all(topology, positions)
    findings += steric.clashes_python(topology, positions, clash_warn_a=0.4, clash_error_a=0.5)
    return {"findings_by_category": dict(Counter(f.category for f in findings)),
            "findings_by_severity": dict(Counter(str(f.severity) for f in findings))}


def pocket_snapshot(topology, before, after):
    """Measure geometric contacts, not affinity or interaction quality."""
    from scipy.spatial import cKDTree

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    protein = [a.index for a in topology.atoms()
               if a.residue.name in PROTEIN_RESIDUES and a.element.symbol != "H"]
    if not protein:
        return []
    results = []
    for residue in topology.residues():
        if residue.name in PROTEIN_RESIDUES | SOLVENT_IONS:
            continue
        indices = [a.index for a in residue.atoms() if a.element.symbol != "H"]
        if not indices:
            continue
        def contacts(coords, ligand_indices):
            nearby = cKDTree(coords[ligand_indices]).query_ball_tree(cKDTree(coords[protein]), r=4.0)
            return {(ligand_indices[i], protein[j]) for i, matches in enumerate(nearby) for j in matches}

        original, current = contacts(before, indices), contacts(after, indices)
        displacement = np.linalg.norm(after[indices] - before[indices], axis=1)
        results.append({"chain": residue.chain.id, "resid": residue.id,
                        "icode": residue.insertionCode.strip(), "resname": residue.name,
                        "heavy_displacement_rms_angstrom": float(np.sqrt(np.mean(displacement**2))),
                        "heavy_displacement_max_angstrom": float(displacement.max()),
                        "contact_cutoff_angstrom": 4.0,
                        "original_contacts": sorted(original), "retained_contacts": sorted(original & current),
                        "lost_contacts": sorted(original - current), "new_contacts": sorted(current - original)})
    return results


def method_arguments(method, charge, uhf, solvent):
    """Keep the production GFN-FF route separate from explicit QM research."""
    if method == "gfn2" and (charge is None or uhf is None):
        raise ValueError("GFN2 research requires explicit --charge and --uhf (including zero)")
    if uhf is not None and uhf < 0:
        raise ValueError("--uhf must be a nonnegative unpaired-electron count")
    flags = ["--gfn", "2"] if method == "gfn2" else ["--gfnff"]
    if charge is not None:
        flags.extend(["--chrg", str(charge)])
    if uhf is not None:
        flags.extend(["--uhf", str(uhf)])
    if solvent != "none":
        flags.extend(["--alpb", solvent])
    return flags


def save_stream(path, data):
    """Keep original bytes even when a failing native executable emits invalid UTF-8."""
    data = data or b""
    path.write_bytes(data)
    text = data.decode("utf-8", errors="replace")
    path.with_suffix(".txt").write_text(text)
    return text


def frozen_atom_spec(indices):
    """Use compact 1-based ranges, matching the production xcontrol writer."""
    numbers = sorted({i + 1 for i in indices})
    ranges = []
    if numbers:
        start = previous = numbers[0]
        for number in numbers[1:]:
            if number != previous + 1:
                ranges.append((start, previous))
                start = number
            previous = number
        ranges.append((start, previous))
    return ",".join(str(a) if a == b else f"{a}-{b}" for a, b in ranges)


def subsystem_components(topology):
    """Connected components in stable atom order; no distance-derived bonds."""
    adjacent = {a.index: set() for a in topology.atoms()}
    for a, b in topology.bonds():
        adjacent[a.index].add(b.index)
        adjacent[b.index].add(a.index)
    remaining, components = set(adjacent), []
    while remaining:
        pending, component = [min(remaining)], set()
        while pending:
            index = pending.pop()
            if index not in component:
                component.add(index)
                pending.extend(adjacent[index] - component)
        remaining -= component
        components.append(sorted(component))
    return components


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--xtb", default=shutil.which("xtb"))
    parser.add_argument("--cycles", type=int, default=200)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--pdb", type=Path, action="append", default=[],
                        help="Already hydrogenated input; extract the production heterogen subsystem")
    parser.add_argument("--fragmented-control", action="store_true",
                        help="Add 60 spatially separated ethanols (540 atoms) to exercise the >500-atom path")
    parser.add_argument("--freeze-control-heavy", action="store_true",
                        help="Mechanism control only: freeze all 180 carbon atoms of the synthetic alkane")
    parser.add_argument("--engines", nargs="+", choices=["default", "rf", "lbfgs", "unfrozen"],
                        default=["default", "rf", "unfrozen"])
    parser.add_argument("--glycam-bonds", action="store_true",
                        help="Use the canonical GLYCAM setup wrapper to populate missing topology bonds")
    parser.add_argument("--case", action="append", default=[], help="Run only exact matching case names")
    parser.add_argument("--component", type=int, help="Research only: optimize one zero-based bonded subsystem component")
    parser.add_argument("--charge", type=int,
                        help="Explicit research charge arm; omitted matches production (no --chrg argument)")
    parser.add_argument("--method", choices=["gfnff", "gfn2"], default="gfnff")
    parser.add_argument("--uhf", type=int, help="Explicit unpaired-electron count; required for GFN2")
    parser.add_argument("--solvent", choices=["none", "water"], default="none",
                        help="ALPB implicit-water sensitivity arm, not an explicit solvent box")
    args = parser.parse_args()
    try:
        method_flags = method_arguments(args.method, args.charge, args.uhf, args.solvent)
    except ValueError as exc:
        parser.error(str(exc))
    if not args.xtb:
        parser.error("xTB executable is required")
    args.xtb = str(Path(args.xtb).resolve())
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.output / "audit-script.py")
    version = subprocess.run([args.xtb, "--version"], capture_output=True, text=True, timeout=30)
    (args.output / "version.txt").write_text(version.stdout + version.stderr)
    report = {
        "platform": platform.platform(), "python": sys.version,
        "binary": args.xtb,
        "binary_sha256": hashlib.sha256(Path(args.xtb).read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "invocation": sys.argv,
        "charge_requested": args.charge,
        "method_requested": args.method, "uhf_requested": args.uhf,
        "solvent_requested": args.solvent,
        "scope": "mechanism controls and optional complex snapshots; not general refinement acceptance",
        "frozen_drift_tolerance_angstrom": 1e-5,
        "threads": 1, "seed": 20261027, "runs": [],
    }
    cases = []
    pdb_contexts = {}
    for case, smiles in [("ethanol", "CCO"), ("two_ethanols", "CCO.CCO")]:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if AllChem.EmbedMolecule(mol, randomSeed=20261027) != 0:
            raise RuntimeError("RDKit embedding failed")
        conf = mol.GetConformer()
        # RDKit embeds fragments independently: separate them explicitly.
        if case == "two_ethanols":
            for idx in Chem.GetMolFrags(mol)[1]:
                point = conf.GetAtomPosition(idx)
                conf.SetAtomPosition(idx, (point.x + 8, point.y, point.z))
        xyz = Chem.MolToXYZBlock(mol)
        frozen = [0, 1, 2]  # three noncollinear heavy atoms of first ethanol
        bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
        cases.append((case, xyz, frozen, bonds))
    if args.fragmented_control:
        template_rows = cases[0][1].splitlines()[2:]
        rows, bonds = [], []
        for fragment in range(60):
            shift = np.array([fragment % 5, (fragment // 5) % 4, fragment // 20]) * 12.0
            for row in template_rows:
                fields = row.split()
                position = np.array([float(x) for x in fields[1:]]) + shift
                rows.append(f"{fields[0]} {position[0]:.8f} {position[1]:.8f} {position[2]:.8f}")
            bonds.extend((a + fragment * 9, b + fragment * 9) for a, b in cases[0][3])
        cases.append(("sixty_ethanols_540_atoms", "540\n60 separated ethanol controls\n" +
                      "\n".join(rows) + "\n", [0, 1, 2], bonds))
        alkane = Chem.MolFromSmiles("C" * 180)
        conf = Chem.Conformer(180)
        for i in range(180):
            conf.SetAtomPosition(i, (i * 1.28, (i % 2) * 0.86, 0.08 * np.sin(i * 0.3)))
        alkane.AddConformer(conf)
        alkane = Chem.AddHs(alkane, addCoords=True)
        alkane_frozen = list(range(180)) if args.freeze_control_heavy else [0, 1, 2]
        cases.append(("connected_alkane_542_atoms", Chem.MolToXYZBlock(alkane), alkane_frozen,
                      [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in alkane.GetBonds()]))
    for index, pdb_path in enumerate(args.pdb):
        from openmm.app import PDBFile

        from dvbfixer.minimize.refine import _extract_heterogen_subsystem, _write_xyz

        pdb = PDBFile(str(pdb_path))
        if args.glycam_bonds:
            from dvbfixer.ffutils import build_glycam_system

            build_glycam_system(["amber14-all.xml", "amber14/GLYCAM_06j-1.xml",
                                 "amber14/tip3pfb.xml"], pdb.topology, pdb.positions, verbose=True)
        top, pos, index_map, frozen = _extract_heterogen_subsystem(pdb.topology, pdb.positions)
        components = subsystem_components(top)
        if args.component is not None:
            from openmm.app import Modeller

            if not 0 <= args.component < len(components):
                raise ValueError("Component index outside extracted subsystem")
            selected = set(components[args.component])
            remap = {old: new for new, old in enumerate(sorted(selected))}
            model = Modeller(top, pos)
            model.delete([a for a in top.atoms() if a.index not in selected])
            top, pos = model.topology, model.positions
            index_map = {full: remap[sub] for full, sub in index_map.items() if sub in selected}
            frozen = {remap[sub] for sub in frozen if sub in selected}
        if top.getNumAtoms() == 0:
            raise ValueError("PDB control requires a nonempty heterogen subsystem")
        case = f"pdb-{index}-{pdb_path.stem}"
        xyz_path = args.output / f"{case}-subsystem.xyz"
        _write_xyz(xyz_path, top, pos)
        metadata = {
            "input": str(pdb_path.resolve()),
            "input_sha256": hashlib.sha256(pdb_path.read_bytes()).hexdigest(),
            "subsystem_atoms": top.getNumAtoms(),
            "frozen_atoms": len(frozen), "full_to_sub_indices": index_map,
            "joint_component_sizes": [len(c) for c in components], "selected_component": args.component,
            "glycam_bonds_requested": args.glycam_bonds,
            "subsystem_atom_identities": [
                {"chain": a.residue.chain.id, "resid": a.residue.id,
                 "icode": a.residue.insertionCode.strip(), "resname": a.residue.name,
                 "atom": a.name, "element": a.element.symbol} for a in top.atoms()
            ],
            "cut_bonds_full_indices": [(a.index, b.index) for a, b in pdb.topology.bonds()
                                       if (a.index in index_map) != (b.index in index_map)],
            "limitation": "Production extraction, no pocket or added caps; not a validated chemical model",
        }
        (args.output / f"{case}-extraction.json").write_text(json.dumps(metadata, indent=2) + "\n")
        cases.append((case, xyz_path.read_text(), sorted(frozen),
                      [(a.index, b.index) for a, b in top.bonds()]))
        from openmm.unit import nanometer

        full_coords = np.asarray(pdb.positions.value_in_unit(nanometer)) * 10
        pdb_contexts[case] = (pdb, index_map, full_coords, geometry_snapshot(pdb.topology, full_coords))
    for case, xyz, frozen, bonds in cases:
        if args.case and case not in args.case:
            continue
        for engine in args.engines:
            work = args.output / f"{case}-{engine}"
            work.mkdir()
            (work / "input.xyz").write_text(xyz)
            atom_list = frozen_atom_spec(frozen)
            control = f"$fix\n  atoms: {atom_list}\n$end\n" if frozen and engine != "unfrozen" else ""
            if engine in ("rf", "lbfgs"):
                control += f"$opt\n  engine={engine}\n$end\n"
            (work / "xcontrol").write_text(control)
            command = [args.xtb, "input.xyz", "--opt", *method_flags, "--cycles",
                       str(args.cycles), "--norestart", "--verbose"]
            if control:
                command.extend(["--input", "xcontrol"])
            env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                   "MKL_NUM_THREADS": "1"}
            record = {"case": case, "engine": engine, "command": command,
                      "frozen_indices_zero_based": frozen if control else [],
                      "tracked_indices_zero_based": frozen}
            start = time.monotonic()
            try:
                run = subprocess.run(command, cwd=work, env=env, capture_output=True,
                                     timeout=args.timeout)
                stdout = save_stream(work / "stdout.raw", run.stdout)
                save_stream(work / "stderr.raw", run.stderr)
                record["returncode"] = run.returncode
                record["converged"] = "GEOMETRY OPTIMIZATION CONVERGED" in stdout
                cycles = re.findall(r"CYCLE\s+(\d+)", stdout)
                record["last_cycle_banner"] = max(map(int, cycles), default=0)
                record["fixed_atoms_acknowledged"] = f"{len(frozen)} fixed atom positions" in stdout
                record["fragmented_diagonalization_observed"] = "fragmented diagonalization" in stdout
                record["optimizer_banner_lines"] = [
                    line.strip() for line in stdout.splitlines()
                    if "Rational Function Optimizer" in line or "L-ANC" in line or "L-BFGS" in line
                ]
                record["fragmentation_log_lines"] = [
                    line for line in stdout.splitlines()
                    if "fragment" in line.lower() or "hessian" in line.lower()
                ]
                coordinate_path = work / "xtbopt.xyz"
                record["final_xyz_present"] = coordinate_path.exists()
                record["coordinate_source"] = "xtbopt.xyz" if coordinate_path.exists() else None
                if not coordinate_path.exists() and (work / "xtbopt.log").exists():
                    partial = last_complete_xyz_frame(work / "xtbopt.log")
                    if partial:
                        coordinate_path = work / "last-trajectory-frame.xyz"
                        coordinate_path.write_text(partial)
                        record["coordinate_source"] = "partial trajectory; not an optimized result"
                if coordinate_path.exists():
                    elements, before = read_xyz(work / "input.xyz")
                    new_elements, after = read_xyz(coordinate_path)
                    if elements != new_elements:
                        raise ValueError("Atom count/element order changed")
                    if frozen:
                        record["raw_frozen_drift"] = drift_metrics(before[frozen], after[frozen])
                        record["frozen_drift_pass"] = record["raw_frozen_drift"]["max_angstrom"] <= 1e-5
                    else:
                        record["fixation_status"] = "not applicable: no covalent protein anchors"
                    record["all_atom_displacement"] = drift_metrics(before, after)
                    merged = after.copy()
                    if control:
                        merged[frozen] = before[frozen]
                    frozen_set = set(frozen)
                    record["frozen_mobile_bonds"] = [
                        {"indices": [a, b], "before_angstrom": float(np.linalg.norm(before[a] - before[b])),
                         "raw_angstrom": float(np.linalg.norm(after[a] - after[b])),
                         "restored_angstrom": float(np.linalg.norm(merged[a] - merged[b]))}
                        for a, b in bonds if (a in frozen_set) != (b in frozen_set)
                    ]
                    if case in pdb_contexts:
                        pdb, index_map, full_coords, baseline_geometry = pdb_contexts[case]
                        raw_full, merged_full = full_coords.copy(), full_coords.copy()
                        for full_index, sub_index in index_map.items():
                            raw_full[full_index] = after[sub_index]
                            merged_full[full_index] = merged[sub_index]
                        record["whole_structure_geometry"] = {
                            "openmm_only": baseline_geometry,
                            "raw_reinserted": geometry_snapshot(pdb.topology, raw_full),
                            "anchors_restored": geometry_snapshot(pdb.topology, merged_full),
                            "note": "Diagnostic snapshots only; unconverged/failing outputs are not accepted structures",
                        }
                        record["heterogen_pocket_contacts"] = pocket_snapshot(pdb.topology, full_coords, merged_full)
            except subprocess.TimeoutExpired as exc:
                record["timeout"] = True
                for name, content in [("stdout", exc.stdout), ("stderr", exc.stderr)]:
                    save_stream(work / f"{name}.raw", content)
            except ValueError as exc:
                record["validation_error"] = str(exc)
            record["seconds"] = time.monotonic() - start
            report["runs"].append(record)
            (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
            print(case, engine, "converged", record.get("converged"),
                  "max drift", record.get("raw_frozen_drift", {}).get("max_angstrom"), flush=True)
    manifest = {str(path.relative_to(args.output)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(args.output.rglob("*")) if path.is_file()}
    (args.output / "manifest.sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
