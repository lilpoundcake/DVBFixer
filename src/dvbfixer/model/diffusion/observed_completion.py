"""Seeded heavy-atom preparation of existing residues, never gap generation.

Observations remain byte-identical. New atoms become fixed conditioning before
context selection/freezing, with independently checked publication evidence.
"""

from __future__ import annotations

import hashlib
import json
import math
import tempfile
from dataclasses import asdict, replace
from io import StringIO
from pathlib import Path

from dvbfixer.diagnose.report import Severity
from dvbfixer.model.diffusion.anchor_policy import atom_residue
from dvbfixer.model.diffusion.contract import ArtifactReference, AtomIdentity, BackendOption
from dvbfixer.model.diffusion.input_geometry import InputGeometryError, assess_input_geometry
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS

REVISION = "seeded-observed-heavy-completion-v1"


def missing_observed_atoms(request, source):
    """Inventory only canonical ATOM residues already present in the input."""
    from dvbfixer.model.diffusion.input_preparation import _records

    if hashlib.sha256(source).hexdigest() != request.normalized_pdb.sha256:
        raise InputGeometryError("heavy-atom preparation source digest mismatch")
    records = _records(source)
    names = {atom_residue(atom): line[17:20].strip() for atom, line in records.items()
             if line.startswith("ATOM  ")}
    missing = []
    for residue, name in names.items():
        for atom_name in CANONICAL_HEAVY_ATOMS.get(name, ()):
            atom = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, atom_name)
            if atom not in records:
                missing.append(atom)
    return tuple(missing)


def _eligible(request, source, missing):
    from dvbfixer.model.diffusion.input_preparation import _records
    from dvbfixer.model.diffusion.scope import _parse_coordinate_records

    if request.protein_context is not None or request.gap_contexts is not None:
        raise InputGeometryError("heavy-atom preparation must precede context selection; frozen requests are never rewritten")
    if any(o.name.startswith("input_preparation.") for o in request.backend_options):
        raise InputGeometryError("an already prepared request cannot be silently prepared again")
    records = _records(source)
    parsed, _, _ = _parse_coordinate_records(source.decode().splitlines())
    if {r.identity for r in parsed if r.element != "H"} != set(request.fixed_atoms):
        raise InputGeometryError("heavy-atom preparation requires every source heavy atom in the fixed mask")
    selected = {atom_residue(a) for a in missing}
    if any(atom_residue(link.atom1) in selected or atom_residue(link.atom2) in selected
           for link in request.retained_explicit_links):
        raise InputGeometryError("heavy-atom preparation refuses explicit links involving incomplete residues")
    placed = {r for p in request.sequence_placements for r in p.observed_residues}
    if not selected <= placed:
        raise InputGeometryError("incomplete residues are outside observed sequence placement")
    for atom, line in records.items():
        name = line[17:20].strip()
        if (line.startswith("HETATM") and name not in {"HOH", "WAT", "DOD"}
                or line.startswith("ATOM  ") and name not in CANONICAL_HEAVY_ATOMS):
            raise InputGeometryError("automatic heavy-atom preparation currently requires canonical protein without nonwater heterogens")
        if line[16].strip():
            raise InputGeometryError("select alternate locations before heavy-atom preparation")
        if atom_residue(atom) in selected:
            if float(line[54:60]) != 1.0:
                raise InputGeometryError("heavy-atom preparation refuses partial occupancy in incomplete residues")
            element = line[76:78].strip().upper()
            if element in {"H", "D"} or (not element and atom.atom_name.lstrip("0123456789").startswith(("H", "D"))):
                raise InputGeometryError("incomplete residues contain observed hydrogens; provide heavy-atom input first")


def _insert(source, added_records):
    from dvbfixer.model.diffusion.atom_completion_materialize import _source_layout
    from dvbfixer.model.diffusion.pdb_materialize import _coordinate_identity

    lines, _serials, last = _source_layout(source.decode())
    additions = {}
    for line in added_records:
        if not isinstance(line, str) or not line.startswith("ATOM  ") or len(line.splitlines()) != 1:
            raise InputGeometryError("invalid heavy-atom preparation record")
        identity = _coordinate_identity(line)
        if identity is None or atom_residue(identity) not in last:
            raise InputGeometryError("heavy-atom preparation cannot insert whole residues")
        additions.setdefault(last[atom_residue(identity)], []).append(line)
    result = []
    for index, line in enumerate(lines):
        if line[:6].strip() == "MASTER":
            continue  # Counts describe the incomplete input, not the prepared structure.
        result.append(line)
        if index in additions:
            if not line.endswith("\n"):
                raise InputGeometryError("cannot append atoms after an unterminated source record")
            newline = "\r\n" if line.endswith("\r\n") else "\n"
            result.extend(item.rstrip("\r\n") + newline for item in additions[index])
    return "".join(result).encode()


def _rebuild(request, source, missing):
    from openmm.unit import angstrom
    from pdbfixer import PDBFixer

    from dvbfixer.ffutils.geometry import assert_all_l, rebuild_missing_atoms_with_retry
    from dvbfixer.model.diffusion.atom_completion_materialize import _source_layout
    from dvbfixer.model.diffusion.input_preparation import _fragmented_forcefield_input, _records
    from dvbfixer.model.diffusion.pdb_materialize import _generated_lines

    # Temporary fragment termini avoid fictitious peptide springs across gaps.
    text = _fragmented_forcefield_input(request, source)
    text = "".join(line for line in text.splitlines(keepends=True) if not line.startswith("HETATM"))
    fixer = PDBFixer(pdbfile=StringIO(text))
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # Whole gaps belong to diffusion, never to this stage.
    fixer.findMissingAtoms()
    rebuild_missing_atoms_with_retry(fixer, verbose=True, log_prefix="pre-diffusion heavy atoms: ")
    assert_all_l(fixer.topology, fixer.positions)
    positions = fixer.positions.value_in_unit(angstrom)
    rebuilt = {}
    for atom in fixer.topology.atoms():
        residue = atom.residue
        key = AtomIdentity(residue.chain.id, residue.id, residue.insertionCode.strip(), atom.name)
        if key in rebuilt:
            raise InputGeometryError("PDBFixer produced ambiguous atom identities")
        rebuilt[key] = (residue.name, atom.element.symbol, tuple(float(v) for v in positions[atom.index]))
    for atom, line in _records(source).items():
        if not line.startswith("ATOM  "):
            continue
        original = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
        if atom not in rebuilt or math.dist(original, rebuilt[atom][2]) > 1e-4:
            raise InputGeometryError("PDBFixer moved or lost an observed atom during heavy-atom preparation")
    if any(atom not in rebuilt for atom in missing):
        raise InputGeometryError("PDBFixer left requested heavy atoms missing")
    _, serials, _ = _source_layout(source.decode())
    available = (serial for serial in range(1, 100_000) if serial not in serials)
    return _generated_lines(missing, rebuilt, available)


def _refine_additions(request, source, raw, missing):
    """If necessary, relax only new atoms; no observed-sidechain reflection."""
    from openmm.app import PDBFile

    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.atom_completion_materialize import _source_layout
    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
    from dvbfixer.model.diffusion.input_preparation import _fragmented_forcefield_input, _records
    from dvbfixer.model.diffusion.pdb_materialize import _generated_lines

    pdb = PDBFile(StringIO(raw.decode()))
    assert_all_l(pdb.topology, pdb.positions)  # Reject raw D, never reflect observations.
    text = _fragmented_forcefield_input(request, raw)
    text = "".join(line for line in text.splitlines(keepends=True) if not line.startswith("HETATM"))
    result = refine_generated_region(text, generated_residues=tuple(dict.fromkeys(atom_residue(a) for a in missing)),
                                     generated_atoms=missing, max_iterations=250, restart_count=8,
                                     platform_name="CPU", random_seed=1)
    if set(result.coordinates_angstrom) != set(missing):
        raise InputGeometryError("heavy-atom refinement returned coordinates outside the addition mask")
    names = {atom_residue(a): line[17:20].strip() for a, line in _records(source).items()}
    represented = {a: (names[atom_residue(a)], a.atom_name[0], tuple(float(v) for v in result.coordinates_angstrom[a]))
                   for a in missing}
    _, serials, _ = _source_layout(source.decode())
    records = _generated_lines(missing, represented, (s for s in range(1, 100_000) if s not in serials))
    return records, {"platform": result.platform, "max_iterations": 250, "restart_count": 8, "seed": 1,
                     "initial_energy_kj_mol": result.initial_energy_kj_mol,
                     "final_energy_kj_mol": result.final_energy_kj_mol, "observed_mobility": False}


def audit_completion_preparation(request, original, prepared, report):
    """Audit additions, immutable observations and geometry without rerunning PDBFixer."""
    from dvbfixer.model.diffusion.input_preparation import _records, preparation_options

    if report.get("revision") != REVISION or report.get("status") != "accepted":
        raise InputGeometryError("unsupported or unsuccessful heavy-atom preparation evidence")
    def digest(data):
        return hashlib.sha256(data).hexdigest()

    if (digest(original) != report["original_sha256"] or digest(prepared) != report["prepared_sha256"]
            or request.normalized_pdb.sha256 != digest(prepared)):
        raise InputGeometryError("heavy-atom preparation source binding mismatch")
    options = preparation_options(request)
    if options is not None and (options["original_sha256"] != digest(original) or options["revision"] != REVISION):
        raise InputGeometryError("heavy-atom preparation original binding mismatch")
    original_records = _records(original)
    added = set(_records(prepared)) - set(original_records)
    before_request = replace(request, normalized_pdb=ArtifactReference("original.pdb", digest(original)),
                             fixed_atoms=tuple(a for a in request.fixed_atoms if a not in added),
                             protein_context=None, gap_contexts=None,
                             backend_options=tuple(o for o in request.backend_options if not o.name.startswith("input_preparation.")))
    missing = missing_observed_atoms(before_request, original)
    _eligible(before_request, original, missing)
    if not missing or added != set(missing) or report["added_atoms"] != [asdict(a) for a in missing]:
        raise InputGeometryError("heavy-atom preparation added mask differs from independent inventory")
    if not added <= set(request.fixed_atoms) or added & set(request.generated_atoms):
        raise InputGeometryError("template-built heavy atoms must be fixed conditioning, not diffusion-generated")
    from dvbfixer.model.diffusion.pdb_materialize import _generated_lines

    names = {atom_residue(a): line[17:20].strip() for a, line in original_records.items()}
    represented = {}
    for atom, line in _records(prepared).items():
        if atom not in added:
            continue
        xyz = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
        if not all(math.isfinite(v) and len(f"{v:8.3f}") == 8 for v in xyz):
            raise InputGeometryError("heavy-atom preparation contains invalid coordinates")
        represented[atom] = (names[atom_residue(atom)], atom.atom_name[0], xyz)
    from dvbfixer.model.diffusion.atom_completion_materialize import _source_layout

    _, serials, _ = _source_layout(original.decode())
    canonical = _generated_lines(missing, represented, (s for s in range(1, 100_000) if s not in serials))
    if report["added_records"] != canonical:
        raise InputGeometryError("heavy-atom preparation records differ from canonical missing-atom identities")
    if prepared != _insert(original, report["added_records"]):
        raise InputGeometryError("heavy-atom preparation changed observations or added unreported records")
    # Reparse complete layout to reject duplicate/reused serials and record ambiguity.
    _source_layout(prepared.decode())
    if missing_observed_atoms(request, prepared):
        raise InputGeometryError("prepared input still has missing canonical heavy atoms")
    if any(report.get(key) is not False for key in ("diffusion_generated", "whole_gaps_built", "observed_atoms_moved")):
        raise InputGeometryError("heavy-atom preparation ownership policy differs from the audited operation")
    raw = _insert(original, report["raw_added_records"])
    raw_request = replace(request, normalized_pdb=ArtifactReference("raw-template.pdb", digest(raw)),
                          protein_context=None, gap_contexts=None)
    raw_findings = assess_input_geometry(raw_request, raw)
    if report["raw_findings"] != [asdict(f) for f in raw_findings]:
        raise InputGeometryError("raw template findings differ from independent geometry assessment")
    from openmm.app import PDBFile

    from dvbfixer.ffutils.geometry import assert_all_l

    pdb = PDBFile(StringIO(prepared.decode()))
    assert_all_l(pdb.topology, pdb.positions)
    findings = assess_input_geometry(request, prepared)
    if any(f.severity is Severity.ERROR for f in findings):
        raise InputGeometryError("heavy-atom preparation left or introduced severe geometry defects")
    return (), findings


def prepare_observed_completion(request, source, missing, work_parent: Path):
    from dvbfixer.model.diffusion.input_preparation import EVIDENCE_FILES, PREFIX

    directory = Path(tempfile.mkdtemp(prefix="dvbfixer-input-completion-", dir=work_parent))
    (directory / EVIDENCE_FILES["original"]).write_bytes(source)
    (directory / "request-before.json").write_text(request.to_json())
    original_digest = hashlib.sha256(source).hexdigest()
    report = {"revision": REVISION, "status": "failed", "original_sha256": original_digest,
              "added_atoms": [asdict(a) for a in missing], "method": "seeded-pdbfixer-heavy-atoms",
              "diffusion_generated": False, "whole_gaps_built": False, "observed_atoms_moved": False}
    try:
        _eligible(request, source, missing)
        print(f"WARNING: Completing {len(missing)} missing heavy atoms in observed residues with seeded PDBFixer "
              f"before diffusion; these are template-built fixed conditioning. Evidence: {directory}", flush=True)
        raw_records = _rebuild(request, source, missing)
        report["added_records"] = raw_records
        prepared = _insert(source, raw_records)
        report["raw_added_records"] = list(raw_records)
        (directory / "raw-template.pdb").write_bytes(prepared)
        raw_request = replace(request, normalized_pdb=ArtifactReference("raw-template.pdb", hashlib.sha256(prepared).hexdigest()),
                              fixed_atoms=(*request.fixed_atoms, *missing))
        raw_findings = assess_input_geometry(raw_request, prepared)
        report["raw_findings"] = [asdict(f) for f in raw_findings]
        report["refinement"] = None
        if any(f.severity is Severity.ERROR for f in raw_findings):
            report["added_records"], report["refinement"] = _refine_additions(request, source, prepared, missing)
            prepared = _insert(source, report["added_records"])
        path = directory / EVIDENCE_FILES["prepared"]
        path.write_bytes(prepared)
        prepared_digest = hashlib.sha256(prepared).hexdigest()
        report.update(status="accepted", prepared_sha256=prepared_digest)
        new_request = replace(request, normalized_pdb=ArtifactReference(path.name, prepared_digest),
                              fixed_atoms=(*request.fixed_atoms, *missing))
        movements, findings = audit_completion_preparation(new_request, source, prepared, report)
        report["movements"] = list(movements)
        report["after_findings"] = [asdict(f) for f in findings]
        report_bytes = (json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode()
        (directory / EVIDENCE_FILES["report"]).write_bytes(report_bytes)
        options = (*request.backend_options, BackendOption(PREFIX + "revision", REVISION),
                   BackendOption(PREFIX + "original_sha256", original_digest),
                   BackendOption(PREFIX + "report_sha256", hashlib.sha256(report_bytes).hexdigest()))
        new_request = replace(new_request, backend_options=options)
        (directory / "request-prepared.json").write_text(new_request.to_json())
        return new_request, path
    except Exception as exc:
        report.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        (directory / EVIDENCE_FILES["report"]).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
        raise InputGeometryError(f"pre-diffusion heavy-atom preparation failed: {exc}. Sampling was not started. "
                                 f"Original and attempt evidence retained: {directory}") from exc
