"""Independent source-geometry checks before diffusion sampling (no mutation)."""

from __future__ import annotations

import hashlib
import json
import tempfile
from dataclasses import asdict
from io import StringIO
from pathlib import Path

from openmm.app import PDBFile, Topology

from dvbfixer.diagnose.chemistry import (
    check_backbone_bond_angles,
    check_bond_lengths,
    check_ca_chirality,
    check_peptide_omegas,
)
from dvbfixer.diagnose.report import Finding, Severity
from dvbfixer.diagnose.steric import clashes_python
from dvbfixer.diagnose.structural import check_coincident_atoms
from dvbfixer.model.diffusion.contract import DiffusionRequest, ResidueIdentity
from dvbfixer.pdbutils.sequence_evidence import observed_peptide_continuity


class InputGeometryError(ValueError):
    """The source cannot safely serve as fixed diffusion conditioning."""


def assess_input_geometry(request: DiffusionRequest, source: bytes) -> tuple[Finding, ...]:
    """Check the measured input, excluding peptide bonds across known gaps.

    A target-sequence gap is not a stretched peptide bond. Conversely, adjacent
    observed sequence positions must not hide a broken peptide behind that rule.
    All original chains participate, not just the cropped neural motif.
    """
    if hashlib.sha256(source).hexdigest() != request.normalized_pdb.sha256:
        raise InputGeometryError("input geometry source digest mismatch")
    # Require the steric engine's dependencies instead of accepting its optional
    # no-dependency fallback (an empty findings list).
    import scipy.spatial  # noqa: F401

    pdb = PDBFile(StringIO(source.decode("utf-8")))
    indices = {residue: (placement.chain, index)
               for placement in request.sequence_placements
               for residue, index in zip(placement.observed_residues, placement.observed_target_indices, strict=True)}
    topology = Topology()
    atoms = {}
    for chain in pdb.topology.chains():
        copied_chain = topology.addChain(chain.id)
        for residue in chain.residues():
            copied_residue = topology.addResidue(residue.name, copied_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                atoms[atom] = topology.addAtom(atom.name, atom.element, copied_residue, atom.id)
    for a, b in pdb.topology.bonds():
        if a.residue is not b.residue and {a.name, b.name} == {"C", "N"}:
            left = indices.get(ResidueIdentity(a.residue.chain.id, a.residue.id, a.residue.insertionCode.strip()))
            right = indices.get(ResidueIdentity(b.residue.chain.id, b.residue.id, b.residue.insertionCode.strip()))
            if left and right and (left[0] != right[0] or abs(left[1] - right[1]) != 1):
                continue
        topology.addBond(atoms[a], atoms[b])
    findings = []
    for check in (check_coincident_atoms, check_bond_lengths, check_backbone_bond_angles,
                  check_peptide_omegas, check_ca_chirality):
        findings.extend(check(topology, pdb.positions))
    findings.extend(clashes_python(topology, pdb.positions))
    lines = source.decode("utf-8").splitlines()
    for placement in request.sequence_placements:
        residues = placement.observed_residues
        continuity = observed_peptide_continuity(lines, placement.chain,
                                                 [(r.residue_number, r.insertion_code) for r in residues])
        for i, connected in enumerate(continuity):
            if connected is False and placement.observed_target_indices[i + 1] == placement.observed_target_indices[i] + 1:
                left_residue, right_residue = residues[i:i + 2]
                findings.append(Finding(Severity.ERROR, "chain_break", left_residue.chain,
                                        left_residue.residue_number + left_residue.insertion_code, "",
                                        f"broken observed peptide before {right_residue.residue_number}{right_residue.insertion_code}; no target-sequence gap accounts for it"))
    return tuple(findings)


def require_input_geometry(request: DiffusionRequest, source: bytes, work_parent: Path) -> None:
    """Retain findings and stop on errors; never silently minimize fixed atoms."""
    findings = assess_input_geometry(request, source)
    report = record_input_geometry(request, source, findings, work_parent)
    errors = [finding for finding in findings if finding.severity is Severity.ERROR]
    if errors:
        first = errors[0]
        raise InputGeometryError(f"diffusion input has {len(errors)} severe geometry finding(s); "
                                 f"{first.chain}:{first.resid}: {first.message}. "
                                 f"Sampling was not started. Retained report: {report}. "
                                 "Repair the input before using it as fixed conditioning.")


def record_input_geometry(request, source, findings, work_parent, *, repair_planned=False):
    """Keep exact source/findings and notify without granting repair acceptance."""
    if not findings:
        return None
    directory = Path(tempfile.mkdtemp(prefix="dvbfixer-input-geometry-", dir=work_parent))
    (directory / "source.pdb").write_bytes(source)
    report = directory / "report.json"
    errors = [finding for finding in findings if finding.severity is Severity.ERROR]
    report.write_text(json.dumps({"schema_version": 1, "source_sha256": request.normalized_pdb.sha256,
                                 "passed": not errors, "coordinates_modified": False,
                                 "findings": [asdict(finding) for finding in findings]},
                                ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    warnings = [finding for finding in findings if finding.severity is Severity.WARNING]
    level = "ERROR" if errors and not repair_planned else "WARNING" if errors or warnings else "INFO"
    print(f"{level}: Diffusion input geometry: {len(errors)} error(s), {len(warnings)} warning(s). "
          "Input coordinates were not modified.", flush=True)
    # Keep the console/GUI job log readable, without hiding the full findings.
    ordered = sorted(findings, key=lambda finding: finding.sort_key())
    for finding in ordered[:5]:
        severity = "WARNING" if repair_planned and finding.severity is Severity.ERROR else finding.severity
        print(f"{severity}: Input {finding.chain or '_'}:{finding.resid} "
              f"{finding.resname} [{finding.category}]: {finding.message}", flush=True)
    if len(ordered) > 5:
        print(f"{level}: {len(ordered) - 5} additional input geometry findings are in the report.", flush=True)
    print(f"{level}: Full input geometry report: {report}", flush=True)
    return report
