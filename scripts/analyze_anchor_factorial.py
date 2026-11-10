#!/usr/bin/env python3
"""Declared, read-only factorial and label-symmetry analysis of retained pilots."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from dataclasses import asdict
from io import StringIO
from pathlib import Path


def main():
    import numpy as np
    from openmm.app import PDBFile
    from research_anchor_calibration import parse_model

    from dvbfixer.diagnose.chemistry import check_bond_lengths, check_ca_chirality
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.model.diffusion.anchor_calibration import compare_anchor
    from dvbfixer.model.diffusion.anchor_symmetry import measure_label_sensitivity
    from dvbfixer.model.diffusion.contract import ResidueIdentity
    from dvbfixer.model.diffusion.geometry import weighted_kabsch
    from dvbfixer.model.diffusion.portfolio_evidence import json_bytes, strict_json
    from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--declaration", type=Path, default=Path("docs/research/anchor-calibration-factorial-v1.json"))
    parser.add_argument("--inputs-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    declaration_bytes = args.declaration.read_bytes()
    declaration = strict_json(declaration_bytes)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "declaration.json").write_bytes(declaration_bytes)
    (args.output / "run-script.py").write_bytes(Path(__file__).read_bytes())
    import dvbfixer.model.diffusion.anchor_symmetry as symmetry_module

    (args.output / "anchor_symmetry.py").write_bytes(Path(symmetry_module.__file__).read_bytes())
    grids = list(itertools.product(declaration["backbone_rmsd_grid_angstrom"],
        declaration["ca_displacement_grid_angstrom"], declaration["heavy_maximum_grid_angstrom"]))
    if len(grids) != 36 or any(not np.isfinite(x) or x <= 0 for triplet in grids for x in triplet):
        raise ValueError("invalid predeclared factorial grid")
    summaries, custody = [], {}
    for group in declaration["inputs"]:
        root = args.inputs_root / group["directory"]
        if root.is_symlink() or not root.is_dir():
            raise ValueError("calibration inputs must be regular retained directories")
        root = root.resolve(strict=True)

        def read(name, directory=root):
            data = read_private_file_bytes(directory, name, max_bytes=16_000_000)
            custody[str(directory / name)] = hashlib.sha256(data).hexdigest()
            return data

        bound = strict_json(read("input-custody.json"))
        old_declaration = strict_json(read("declaration.json"))
        for name, expected in bound.items():
            data = read(name)
            if len(data) != expected["bytes"] or hashlib.sha256(data).hexdigest() != expected["sha256"]:
                raise ValueError("retained experimental input differs from its original custody")
        for case_index in group["cases"]:
            original = strict_json(read(f"case-{case_index}.json"))
            case = original["case"]
            if case != old_declaration["cases"][case_index]:
                raise ValueError("retained comparison differs from original declaration")
            reference_text = read(case["reference"] + ".pdb").decode("ascii")
            mobile_text = read(case["mobile"] + ".pdb").decode("ascii")
            ref, names, bad_ref, order, ref_model = parse_model(reference_text, case["reference_model"], old_declaration["chain"])
            mobile, mobile_names, bad_mobile, mobile_order, mobile_model = parse_model(mobile_text, case["mobile_model"], old_declaration["chain"])
            if order != mobile_order or names != mobile_names:
                raise ValueError("retained pilot no longer has an exact residue mapping")
            valid = set(order) - bad_ref.keys() - bad_mobile.keys()

            def residue_of(atom):
                return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)

            valid = {r for r in valid if {a for a in ref if residue_of(a) == r} ==
                     {a for a in mobile if residue_of(a) == r}}
            shared = tuple(sorted(a for a in ref.keys() & mobile.keys() if residue_of(a) in valid))
            rows = []
            for old in original["measurements"]:
                residue = ResidueIdentity(**old["residue"])
                position = order.index(residue)
                fit_residues = valid - {residue}
                if old["frame"] == "local-window-4-leave-anchor-out":
                    fit_residues &= set(order[max(0, position - 4):position + 5])
                elif old["frame"] != "global-leave-anchor-out":
                    raise ValueError("unknown retained fitting frame")
                fit = tuple(a for a in shared if residue_of(a) in fit_residues and a.atom_name in {"N", "CA", "C", "O"})
                contacts = tuple(a for a in shared if residue_of(a) not in set(order[max(0, position - 1):position + 2]))
                exact = compare_anchor(ref, mobile, residue, fit, contact_atoms=contacts,
                    contact_cutoff_angstrom=old_declaration["contact_cutoff_angstrom"])
                if strict_json(json_bytes(asdict(exact))) != {k: v for k, v in old.items() if k != "frame"}:
                    raise ValueError("exact-name recomputation differs from retained measurement")
                transform = weighted_kabsch(np.array([mobile[a] for a in fit]), np.array([ref[a] for a in fit]))
                atoms = tuple(a for a in shared if residue_of(a) == residue)
                reference_axis = {a.atom_name: ref[a] for a in atoms}
                transformed = transform.apply(np.array([mobile[a] for a in atoms]))
                mobile_axis = {a.atom_name: point for a, point in zip(atoms, transformed, strict=True)}
                row = dict(old)
                try:
                    measurement = measure_label_sensitivity(names[residue], reference_axis, mobile_axis)
                except ValueError as exc:
                    row["symmetry"] = {"status": "unavailable", "reason": str(exc)}
                else:
                    if not np.isclose(measurement.exact_maximum_angstrom, exact.maximum_heavy_displacement_angstrom, atol=1e-10, rtol=0):
                        raise ValueError("symmetry frame differs from exact-name measurement")
                    row["symmetry"] = {"status": "available", **asdict(measurement)}
                rows.append(row)
            qc = []
            for role, text in (("reference", ref_model), ("mobile", mobile_model)):
                structure = PDBFile(StringIO(text))
                findings = [*check_ca_chirality(structure.topology, structure.positions),
                            *check_bond_lengths(structure.topology, structure.positions)]
                qc.append({"role": role, "error_count": sum(f.severity is Severity.ERROR for f in findings),
                    "findings": [asdict(f) for f in findings]})
            factorial = []
            for frame in declaration["frames"]:
                frame_rows = [r for r in rows if r["frame"] == frame]
                available = [r for r in frame_rows if r["symmetry"]["status"] == "available"]
                for backbone, ca, heavy in grids:
                    def passed(row, symmetry=False, bb_limit=backbone, ca_limit=ca, heavy_limit=heavy):
                        maximum = (row["symmetry"]["symmetry_minimax_angstrom"] if symmetry
                                   else row["maximum_heavy_displacement_angstrom"])
                        return row["backbone_rmsd_angstrom"] <= bb_limit and row["ca_displacement_angstrom"] <= ca_limit and maximum <= heavy_limit

                    factorial.append({"frame": frame, "backbone_max": backbone, "ca_max": ca, "heavy_max": heavy,
                        "exact_total": len(frame_rows), "exact_passed": sum(passed(r) for r in frame_rows),
                        "paired_total": len(available), "exact_paired_passed": sum(passed(r) for r in available),
                        "symmetry_paired_passed": sum(passed(r, True) for r in available)})
            label = group["directory"] + f"-case-{case_index}"
            report = {"case": case, "measurements": rows, "factorial": factorial, "source_qc": qc,
                "source_defect_stratum": any(r["error_count"] for r in qc), "production_limits_changed": False}
            (args.output / (label + ".json")).write_bytes(json_bytes(report))
            baseline = [r for r in factorial if r["backbone_max"] == 0.5 and r["ca_max"] == 0.5 and r["heavy_max"] == 1.25]
            summary = {"case": case, "source_defect_stratum": report["source_defect_stratum"], "baseline_grid": baseline}
            summaries.append(summary)
            print(json.dumps(summary), flush=True)
    (args.output / "custody.json").write_bytes(json_bytes(custody))
    (args.output / "summary.json").write_bytes(json_bytes({"scope": "retrospective label/factorial sensitivity; no adopted ceiling",
        "declaration_sha256": hashlib.sha256(declaration_bytes).hexdigest(), "cases": summaries,
        "production_limits_changed": False, "coordinate_or_label_rewriting": False}))


if __name__ == "__main__":
    main()
