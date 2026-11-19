"""Frozen offline CCD ideal-roundtrip controls, NOT native-accuracy qualification."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any

from dvbfixer.chemical_reconstruction.ccd_authority import load_geometry_authority, verify_cip
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Observation,
    Operation,
    Ownership,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.plan_compiler import PlanCompilationError, compile_chain_plan
from dvbfixer.chemical_reconstruction.proposals import AuditPolicy
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures/chemical_reconstruction/ccd"
COMPONENTS = ("ALA", "THR", "ILE", "LYS", "PHE", "EOH", "ACT")
HEAVY_MASKS = {"ALA": ("CB",), "THR": ("CB", "OG1", "CG2"), "ILE": ("CB", "CG1", "CG2", "CD1"),
               "LYS": ("CB", "CG", "CD", "CE", "NZ"), "PHE": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ")}
POLICY = AuditPolicy(.02, 2., .7, .02)
PROTOCOL_SHA256 = "89f60a0b0ce8d2ed69a2db0a23a63f4ba52d78f48647f5c0d2e65e51abc57459"
SOURCE_MANIFEST_SHA256 = "29c4f83b28435a34b666b5c1489c4f3423c53beb276e75d16cd0d8c3c055dece"


def _json(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _json(asdict(value))
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(x) for x in value]
    if isinstance(value, frozenset):
        return sorted((_json(x) for x in value), key=lambda x: json.dumps(x, sort_keys=True))
    return value


def _write(path: Path, value: Any) -> None:
    with path.open("x") as handle:
        handle.write(json.dumps(_json(value), indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)+"\n")


def _transform(point: tuple[float, float, float], variant: str) -> tuple[float, float, float]:
    return point if variant == "identity" else (point[2]+7., point[0]-4., point[1]+2.)


def run_controls(output: Path) -> dict[str, Any]:
    """One complete declared run; never overwrite an earlier attempt's evidence."""
    manifest_bytes = (FIXTURES / "manifest.json").read_bytes()
    if hashlib.sha256(manifest_bytes).hexdigest() != SOURCE_MANIFEST_SHA256:
        raise ValueError("frozen source pins changed; declare a new protocol, do not silently replay")
    pins = json.loads(manifest_bytes)["components"]
    protocol = ROOT / "docs/plans/chemical-reconstruction-ccd-v1.md"
    if hashlib.sha256(protocol.read_bytes()).hexdigest() != PROTOCOL_SHA256:
        raise ValueError("frozen protocol changed")
    output.mkdir(parents=True, exist_ok=False)
    controls = []
    for name in COMPONENTS:
        try:
            authority = load_geometry_authority(FIXTURES / f"{name}.cif", name, expected_sha256=pins[name]["sha256"])
        except ValueError as exc:
            for arm in ("heavy", "hydrogens"):
                for variant in ("identity", "cyclic-translate"):
                    item = {"id": f"{name}-{arm}-{variant}", "status": "authority-refused", "reason": str(exc)}
                    _write(output / f"{item['id']}.json", item)
                    controls.append(item)
            continue
        graph = authority.state.graph
        identity = ComponentInstanceRef(1, "A", "82", "B", occurrence=2)
        for arm in ("heavy", "hydrogens"):
            if arm == "hydrogens":
                missing = {a.name for a in graph.atoms if a.element == "H"}
            elif name in HEAVY_MASKS:
                missing = set(HEAVY_MASKS[name])
            else:
                oxygens = sorted(a.name for a in graph.atoms if a.element == "O")
                missing = {oxygens[-1]}
            if not missing <= {a.name for a in graph.atoms}:
                raise ValueError(f"frozen mask absent from {name} source")
            for variant in ("identity", "cyclic-translate"):
                observations = tuple(Observation(ExactAtomRef(identity, a.name), a.element, _transform(a.ideal_position, variant),
                                                 Unit.ANGSTROM) for a in graph.atoms
                                     if a.name not in missing and a.element != "H")
                request = ReconstructionRequest(
                    digest_record(("ccd-derived-roundtrip-input-v1", authority.component.provenance.sha256, arm, variant)),
                    Operation.HYDROGENS if arm == "hydrogens" else Operation.HEAVY,
                    (Component(identity, ChemistryKind.ISOLATED, (authority.state,)),), observations,
                    frozenset({identity}), Ownership(frozenset(ExactAtomRef(identity, n) for n in missing)),
                    constraints=authority.constraints(identity, include_hydrogens=arm == "hydrogens",
                                                     minimum_volume=.05, planar_tolerance_angstrom=.02),
                )
                before = digest_record(observations)
                item: dict[str, Any] = {"id": f"{name}-{arm}-{variant}", "request_digest": request.digest,
                                        "authority_digest": authority.digest, "source_sha256": pins[name]["sha256"],
                                        "rdkit_version": authority.rdkit_version, "requested_atoms": sorted(missing),
                                        "expected_negative": arm == "heavy" and name in {"EOH", "ACT"}}
                evidence: dict[str, Any] = {"control": item, "authority": authority, "request": request}
                try:
                    plan = compile_chain_plan(request, POLICY)
                except PlanCompilationError as exc:
                    item.update(status="compiler-refused", reason=exc.code, blocked=_json(exc.blocked))
                else:
                    result = propose(request, plan)
                    item.update(status=result.status, reason=result.reason)
                    evidence.update(plan=plan, result=result)
                    errors = {}
                    for stage in ("raw", "rounded"):
                        candidate = getattr(result, stage)
                        if candidate is None:
                            continue
                        target = {a.name: _transform(a.ideal_position, variant) for a in graph.atoms}
                        deviations = [math.dist(o.position, target[o.identity.name]) for o in candidate.additions]
                        errors[stage] = max(deviations) if deviations else None
                        if arm == "hydrogens":
                            positions = {o.identity.name: o.position for o in observations}
                            positions.update({o.identity.name: o.position for o in candidate.additions})
                            try:
                                labels = verify_cip(graph, positions)
                            except ValueError as exc:
                                item[f"{stage}_cip"] = {"passed": False, "error": str(exc)}
                            else:
                                item[f"{stage}_cip"] = {"passed": True, "labels": labels}
                    item["authority_roundtrip_max_error_angstrom"] = errors
                item["observations_unchanged"] = digest_record(request.observations) == before
                _write(output / f"{item['id']}.json", evidence)
                controls.append(item)
    summary = {
        "schema": "dvbfixer.ccd-geometry-roundtrip-controls.v1", "protocol_commit_before_download": "67c99e8",
        "protocol_sha256": hashlib.sha256(protocol.read_bytes()).hexdigest(), "controls": controls,
        "declared": 28, "executed_entries": len(controls),
        "geometry_passed": sum(c["status"] == "declared-geometry-passed" for c in controls),
        "compiler_refused": sum(c["status"] == "compiler-refused" for c in controls),
        "held_out_or_native_accuracy_qualification": False, "polymer_state_or_h_qualification": False,
        "production_enabled": False,
    }
    _write(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new private evidence directory; never overwritten")
    args = parser.parse_args()
    result = run_controls(args.output)
    print(json.dumps({k: v for k, v in result.items() if k != "controls"}, indent=2))
