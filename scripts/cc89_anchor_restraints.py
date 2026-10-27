"""Research-only paired positional-restraint ablation; never publish candidates."""

from __future__ import annotations

import argparse
import shutil
import time
from collections import defaultdict
from dataclasses import asdict
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest, ResidueIdentity

BACKBONE = frozenset({"N", "CA", "C", "O"})
POLICIES = (("fixed", 0), ("free", 0), ("uniform", 1000), ("bfactor", 1000),
            ("uniform", 10000), ("bfactor", 10000))
SOURCE_SHA = "6bcc9638c3d299ea12151c3894b7e6c3482ebd2fdfe1d8e54dcc1f6f6316a281"
RAW_SHA = {"reinjection-off": "d0dae11433d3da6318faec026d9062c5d04d5be6602ce55f00bacf6df11067bf",
           "reinjection-on": "7453570841f46986a782d74501c730cdb01c56f34a5ebb190aa0e00745fd95a5"}


def residue(atom):
    return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)


def bfactor_weights(source_text, anchor_atoms):
    """Mean-one atom-weighted weights from original residue backbone B values."""
    values = defaultdict(dict)
    for line in source_text.splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() in BACKBONE:
            key = ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())
            name = line[12:16].strip()
            if name in values[key]:
                raise ValueError("ambiguous backbone atom/alternate state")
            values[key][name] = float(line[60:66])
    if not anchor_atoms or len(set(anchor_atoms)) != len(anchor_atoms):
        raise ValueError("anchor atoms must be nonempty and unique")
    selected = {residue(a) for a in anchor_atoms}
    means = {}
    for key in selected:
        if set(values[key]) != BACKBONE:
            raise ValueError("anchor needs complete original backbone B values")
        mean = float(np.mean(list(values[key].values())))
        if not np.isfinite(mean) or mean <= 0 or any(not np.isfinite(b) or b <= 0 for b in values[key].values()):
            raise ValueError("anchor B values must be finite and positive")
        means[key] = mean
    median = float(np.median(list(means.values())))
    raw = np.array([np.clip(median / means[residue(a)], 0.5, 2.0) for a in anchor_atoms])
    return dict(zip(anchor_atoms, raw / raw.mean(), strict=True)), means


def add_anchor_force(system, topology, source_coordinates, coefficients):
    """Source coordinates are Å; OpenMM per-particle centres are nm, k is kJ/nm²."""
    from openmm import CustomExternalForce

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    for name in ("k", "x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    seen = set()
    for atom in topology.atoms():
        identity = AtomIdentity(atom.residue.chain.id, str(atom.residue.id),
                                atom.residue.insertionCode.strip(), atom.name)
        if identity in coefficients:
            k = coefficients[identity]
            xyz = np.asarray(source_coordinates[identity], dtype=float)
            if k <= 0 or not np.isfinite(k) or not np.isfinite(xyz).all():
                raise ValueError("invalid restraint parameters")
            if identity in seen:
                raise ValueError("duplicate topology atom")
            seen.add(identity)
            force.addParticle(atom.index, [k, *(xyz / 10)])
    if seen != set(coefficients):
        raise ValueError("restraint identities missing from final topology")
    system.addForce(force)
    return force


def assess(path, request, source_text, mobile_atoms):
    from openmm.app import PDBFile
    from scripts.run_cc89_lrd_ablation import atom_records

    from dvbfixer.diagnose.chemistry import (
        check_backbone_bond_angles,
        check_bond_lengths,
        check_peptide_omegas,
    )
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.diagnose.steric import clashes_python
    from dvbfixer.model.diffusion.geometry import weighted_kabsch
    from dvbfixer.model.diffusion.quality import check_ramachandran, check_sidechain_chi12
    from dvbfixer.model.diffusion.validate import (
        ValidationThresholds,
        _parse_pdb_text,
        _Structure,
        _validate_candidate,
    )

    def load(text, filename):
        atoms, residues, links = _parse_pdb_text(text)
        return _Structure(filename, PDBFile(StringIO(text)), atoms, residues, frozenset(links), text)

    source = load(source_text, Path("source.pdb"))
    candidate = load(path.read_text(), path)
    summary, _ = _validate_candidate(request, source, candidate, ValidationThresholds())
    geometry_gates = [g for g in summary.hard_gate_failures
                      if g not in {"fixed-heavy-atom-rmsd", "fixed-heavy-atom-max-displacement"}]
    # Audit the SAME generated+anchor neighborhood for every policy, plus one
    # target-position residue at each edge. Author-number arithmetic is not used.
    anchors = {r for gap in request.gaps for r in (gap.left_anchor, gap.right_anchor) if r is not None}
    impacted = set()
    for placement in request.sequence_placements:
        positions = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
        for gap in request.gaps:
            if gap.chain == placement.chain:
                positions.update(zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True))
        active = {residue(a) for a in request.generated_atoms} | anchors
        indices = {i for i, r in positions.items() if r in active}
        expanded = {j for i in indices for j in (i - 1, i, i + 1)}
        impacted.update((r.chain, r.residue_number + r.insertion_code) for i, r in positions.items() if i in expanded)
    findings = []
    for check in (check_bond_lengths, check_backbone_bond_angles, check_peptide_omegas, check_ramachandran, check_sidechain_chi12):
        for f in check(candidate.pdb.topology, candidate.pdb.positions):
            if (f.chain, f.resid) in impacted and (f.severity is Severity.ERROR or f.category == "non_planar_amide"):
                findings.append(asdict(f))
    for f in clashes_python(candidate.pdb.topology, candidate.pdb.positions, clash_warn_a=0.9, clash_error_a=0.9):
        partner = f.extra.get("clash_partner", {})
        if (f.chain, f.resid) in impacted or (partner.get("chain"), partner.get("resid")) in impacted:
            findings.append(asdict(f))
    records, original = atom_records(path.read_text()), atom_records(source_text)
    fit = tuple(a for a in request.fixed_atoms if a.chain == "W" and a.atom_name in BACKBONE)
    if len(fit) != 844:
        raise ValueError("unexpected W-chain fitting population")
    transform = weighted_kabsch(np.array([records[a][0] for a in fit]), np.array([original[a][0] for a in fit]))
    anchor_metrics = {}
    all_sq, all_raw_sq = [], []
    for r in sorted(anchors):
        heavy = [a for a in request.fixed_atoms if residue(a) == r]
        bb = [a for a in heavy if a.atom_name in BACKBONE]
        xyz, ref = np.array([records[a][0] for a in bb]), np.array([original[a][0] for a in bb])
        sq = np.sum((transform.apply(xyz) - ref)**2, axis=1)
        raw_sq = np.sum((xyz - ref)**2, axis=1)
        all_sq.extend(sq)
        all_raw_sq.extend(raw_sq)
        anchor_metrics[f"{r.chain}{r.residue_number}{r.insertion_code}"] = {
            "backbone_rmsd_W_aligned_A": float(np.sqrt(np.mean(sq))),
            "backbone_rmsd_original_frame_A": float(np.sqrt(np.mean(raw_sq))),
            "max_heavy_displacement_original_frame_A": float(max(np.linalg.norm(records[a][0] - original[a][0]) for a in heavy)),
            "max_heavy_displacement_W_aligned_A": float(max(np.linalg.norm(transform.apply(np.array([records[a][0]]))[0] - original[a][0]) for a in heavy)),
        }
    outside = [a for a in request.fixed_atoms if a not in mobile_atoms]
    if any(records[a][2] != original[a][2] for a in outside):
        raise ValueError("nonmobile source atom record changed")
    return {"strict_validation": asdict(summary), "geometry_gates_excluding_source_drift": geometry_gates,
            "expanded_geometry_findings": findings, "geometry_only_pass_not_acceptance": not geometry_gates and not findings,
            "anchors": anchor_metrics, "combined_anchor_backbone_rmsd_W_aligned_A": float(np.sqrt(np.mean(all_sq))),
            "combined_anchor_backbone_rmsd_original_frame_A": float(np.sqrt(np.mean(all_raw_sq))),
            "nonmobile_atom_records_exact": True, "fit_atom_count": len(fit),
            "rotation": transform.rotation.tolist(), "translation": transform.translation.tolist()}


def main():
    from scripts.run_cc89_lrd_ablation import (
        atom_records,
        digest,
        rewrite_generated_coordinates,
        write_json,
    )

    from dvbfixer.model.diffusion import boundary_refinement as br

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    base, out = args.baseline.resolve(), args.output.resolve()
    request = DiffusionRequest.from_json((base / "request.json").read_text())
    source_path = base / request.normalized_pdb.path
    if digest(source_path) != SOURCE_SHA or request.seeds != (7,) or request.candidate_count != 1:
        raise ValueError("unexpected source/request")
    for arm, sha in RAW_SHA.items():
        if digest(base / "initial" / arm / "candidate.pdb") != sha:
            raise ValueError("raw candidate digest changed")
    out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(source_path, out / "source.pdb")
    shutil.copy2(base / "request.json", out / "request.json")
    source_text = source_path.read_text()
    original = atom_records(source_text)
    anchors = {r for gap in request.gaps for r in (gap.left_anchor, gap.right_anchor) if r is not None}
    anchor_atoms = tuple(a for a in request.fixed_atoms if residue(a) in anchors)
    if len(anchor_atoms) != 24 or len(request.generated_atoms) != 30:
        raise ValueError("unexpected mobility mask")
    weights, b_means = bfactor_weights(source_text, anchor_atoms)
    # Snapshot the exact research code, production refiner and validator before runs.
    snapshots = [Path(__file__), Path(br.__file__), Path("scripts/run_cc89_lrd_ablation.py"),
                 Path("docs/research/9cqd-bfactor-anchor-restraints.md")]
    snapshot_dir = out / "code"
    snapshot_dir.mkdir()
    for p in snapshots:
        shutil.copy2(p, snapshot_dir / p.name)
    import importlib.metadata
    manifest = {"source_sha256": SOURCE_SHA, "request_sha256": digest(base / "request.json"),
        "raw_sha256": RAW_SHA, "code_sha256": {p.name: digest(p) for p in snapshots},
        "environment": {name: importlib.metadata.version(name) for name in ("openmm", "pdbfixer", "numpy")},
        "platform": "CPU", "threads": 1, "solvent": False, "seed": 7,
        "max_iterations": 500, "restart_count": 8, "perturbation_A": 0.25,
        "policies": POLICIES, "production_revision": br.BOUNDARY_REFINEMENT_REVISION,
        "restraints": [{"atom": asdict(a), "relative_weight": float(weights[a]),
                         "source_backbone_B": b_means[residue(a)]} for a in anchor_atoms]}
    write_json(out / "manifest.json", manifest)
    report = []
    original_freeze = br._freeze_observed_heavy_atoms
    for arm in RAW_SHA:
        input_path = base / "initial" / arm / "candidate.pdb"
        for policy, strength in POLICIES:
            directory = out / arm / f"{policy}-{strength}"
            directory.mkdir(parents=True)
            shutil.copy2(input_path, directory / "input.pdb")
            mobile = tuple(request.generated_atoms) if policy == "fixed" else tuple((*request.generated_atoms, *anchor_atoms))
            mobile_residues = tuple(dict.fromkeys(residue(a) for a in mobile))
            coefficients = {a: strength * (weights[a] if policy == "bfactor" else 1.0) for a in anchor_atoms} if strength else {}
            def freeze_and_restrain(system, topology, coordinates, generated, coefficients=coefficients):
                original_freeze(system, topology, coordinates, generated)
                if coefficients:
                    add_anchor_force(system, topology, {a: original[a][0] for a in anchor_atoms}, coefficients)
            row = {"arm": arm, "policy": policy, "mean_k_kJ_mol_nm2": strength,
                   "input_sha256": digest(input_path), "candidate": str((directory / "candidate.pdb").relative_to(out))}
            started = time.perf_counter()
            try:
                with patch.object(br, "_freeze_observed_heavy_atoms", freeze_and_restrain):
                    result = br.refine_generated_region(input_path.read_text(), generated_residues=mobile_residues,
                        generated_atoms=mobile, max_iterations=500, restart_count=8,
                        perturbation_angstrom=0.25, platform_name="CPU", random_seed=7)
                text = rewrite_generated_coordinates(input_path.read_text(), result.coordinates_angstrom)
                before, after = atom_records(input_path.read_text()), atom_records(text)
                if any(before[a][2] != after[a][2] for a in before if a not in set(mobile)):
                    raise ValueError("changed input atom outside mobility mask")
                path = directory / "candidate.pdb"
                path.write_text(text)
                row.update({"candidate_sha256": digest(path), "initial_energy_kJ_mol": result.initial_energy_kj_mol,
                            "final_energy_kJ_mol": result.final_energy_kj_mol, **assess(path, request, source_text, set(mobile))})
            except Exception as exc:
                row["error"] = f"{type(exc).__name__}: {exc}"
            row["wall_seconds"] = time.perf_counter() - started
            write_json(directory / "result.json", row)
            report.append(row)
            write_json(out / "comparison.json", report)
            print(arm, policy, strength, "geometry", row.get("geometry_only_pass_not_acceptance"),
                  "RMSD", row.get("combined_anchor_backbone_rmsd_W_aligned_A"), row.get("error", ""), flush=True)
    write_json(out / "completed.json", {"rows": len(report), "errors": sum("error" in row for row in report)})


if __name__ == "__main__":
    main()
