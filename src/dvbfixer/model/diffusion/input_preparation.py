"""Observed-heavy completion or bounded repair, before selection and freezing.

No diffusion-generated masks are changed; template completions join fixed atoms.
Original bytes, failed attempts and
the independent before/after audit are retained even if sampling never starts.
"""

from __future__ import annotations

import hashlib
import json
import math
import tempfile
from dataclasses import asdict, replace
from pathlib import Path

from dvbfixer.diagnose.report import Severity
from dvbfixer.model.diffusion.anchor_policy import (
    AnchorMovementLimits,
    atom_residue,
    build_anchor_restraints,
    check_anchor_movements,
)
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    BackendOption,
    DiffusionRequest,
    ResidueIdentity,
)
from dvbfixer.model.diffusion.input_geometry import (
    InputGeometryError,
    assess_input_geometry,
    record_input_geometry,
)

LEGACY_REVISION = "bounded-observed-input-repair-v1"
REVISION = "bounded-observed-input-repair-v2"
COMPOSED_REVISION = "completed-and-bounded-input-repair-v1"
PREFIX = "input_preparation."
EVIDENCE_FILES = {"original": "source-original.pdb", "prepared": "prepared.pdb", "report": "input-repair.json"}
MAX_SELECTED_RESIDUES = 64
MAX_SELECTED_FRACTION = 0.25
RESTRAINT_K = 10000.0
MAX_ITERATIONS = 1000
PUBLICATION_ROLES = frozenset({"input-original", "input-prepared", "input-repair-report", "input-repair-request"})


def _digest(data):
    return hashlib.sha256(data).hexdigest()


def preparation_options(request):
    from dvbfixer.model.diffusion.observed_completion import REVISION as COMPLETION_REVISION

    options = [option for option in request.backend_options if option.name.startswith(PREFIX)]
    if not options:
        return None
    values = {option.name.removeprefix(PREFIX): option.value for option in options}
    if (len(values) != len(options) or set(values) != {"revision", "original_sha256", "report_sha256"}
            or values["revision"] not in {LEGACY_REVISION, REVISION, COMPLETION_REVISION, COMPOSED_REVISION}):
        raise InputGeometryError("incomplete or unsupported input preparation binding")
    for key in ("original_sha256", "report_sha256"):
        ArtifactReference("input-evidence", values[key])
    return values


def _records(source):
    records = {}
    for line in source.decode("utf-8").splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom in records:
            raise InputGeometryError("automatic input repair requires unique observed atom identities")
        records[atom] = line
    return records


def _coordinates(records):
    return {atom: tuple(float(line[i:i + 8]) for i in (30, 38, 46)) for atom, line in records.items()}


def _check_selected_links(request, source, selected, *, allow_disulfides):
    """Only proven, unambiguous SG-SG bonds use the canonical protein FF path.

    Arbitrary covalent adducts are still unsupported. An annotated disulfide
    must actually exist in the parsed force-field topology, not merely be near
    the selected residue, and each sulfur must have exactly one partner.
    """
    from io import StringIO

    from openmm.app import PDBFile

    from dvbfixer.model.diffusion.boundary_refinement import _disulfide_bonds_by_identity

    links = [link for link in request.retained_explicit_links
             if atom_residue(link.atom1) in selected or atom_residue(link.atom2) in selected]
    if not links:
        return
    if not allow_disulfides:
        raise InputGeometryError("automatic input repair refuses explicit covalent links involving selected residues")
    records = _records(source)
    pdb = PDBFile(StringIO(source.decode()))
    supported = _disulfide_bonds_by_identity(pdb.topology)
    for link in links:
        pair = frozenset((link.atom1, link.atom2))
        if (len(pair) != 2 or any(atom.atom_name != "SG" or records[atom][17:20].strip() != "CYS" for atom in pair)
                or pair not in supported or any(sum(atom in bond for bond in supported) != 1 for atom in pair)):
            raise InputGeometryError("automatic input repair cannot represent a selected explicit covalent link "
                                     "as an unambiguous force-field disulfide")


def _select_residues(request, source, findings, *, allow_disulfides=True):
    records = _records(source)
    known = {atom_residue(atom) for atom in request.fixed_atoms}
    chosen = set()
    for finding in findings:
        if finding.severity is not Severity.ERROR:
            continue
        if finding.category not in {"clash", "bond_length", "bond_angle", "non_planar_amide"}:
            raise InputGeometryError(f"automatic input repair cannot safely handle {finding.category}; manual preparation is required")
        atoms = finding.extra.get("atoms")
        if not isinstance(atoms, list) or not atoms:
            raise InputGeometryError("input finding lacks exact affected-atom evidence; refusing to guess a repair mask")
        for atom in atoms:
            residue = ResidueIdentity(atom["chain"], atom["resid"], atom["icode"])
            if residue not in known:
                raise InputGeometryError("input finding references a residue outside measured protein atoms")
            chosen.add(residue)
    if not chosen or len(chosen) > MAX_SELECTED_RESIDUES or len(chosen) > MAX_SELECTED_FRACTION * len(known):
        raise InputGeometryError("input damage exceeds the bounded local repair mask; no global relaxation is authorized")
    if any(line.startswith("HETATM") and line[17:20].strip() not in {"HOH", "WAT", "DOD"}
           for line in source.decode().splitlines()):
        raise InputGeometryError("automatic input repair currently requires canonical protein without nonwater heterogens")
    _check_selected_links(request, source, chosen, allow_disulfides=allow_disulfides)
    for atom, line in records.items():
        if atom_residue(atom) not in chosen:
            continue
        element = line[76:78].strip().upper()
        if element in {"H", "D"} or (not element and atom.atom_name.lstrip("0123456789").startswith(("H", "D"))):
            raise InputGeometryError("selected input residues contain measured hydrogens; prepare a heavy-atom input first")
    for residue in chosen:
        if not {AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, name)
                for name in ("N", "CA", "C", "O")} <= records.keys():
            raise InputGeometryError("automatic input repair requires complete observed backbones")
    return tuple(sorted(chosen))


def _fragmented_forcefield_input(request, source):
    """Break inferred peptides across known missing sequence, in FF copy only.

    The temporary force field uses standard fragment termini, not a fictitious
    C-N spring across the missing region. No new termini enter the prepared PDB.
    """
    indices = {residue: index for p in request.sequence_placements
               for residue, index in zip(p.observed_residues, p.observed_target_indices, strict=True)}
    previous = None
    output = []
    for line in source.decode().splitlines(keepends=True):
        if line.startswith("TER"):
            previous = None
        if line.startswith("ATOM  "):
            residue = ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())
            if residue != previous:
                if (previous is not None and previous.chain == residue.chain and
                        previous in indices and residue in indices and indices[residue] != indices[previous] + 1):
                    output.append("TER\n")
                previous = residue
        output.append(line)
    return "".join(output)


def audit_prepared_input(request, original, prepared, report):
    """Recompute mask, movement, source identity and geometry, not report claims."""
    from dvbfixer.model.diffusion.observed_completion import REVISION as COMPLETION_REVISION
    from dvbfixer.model.diffusion.observed_completion import audit_completion_preparation

    if report.get("revision") == COMPOSED_REVISION:
        return audit_composed_preparation(request, original, prepared, report)
    if report.get("revision") == COMPLETION_REVISION:
        return audit_completion_preparation(request, original, prepared, report)
    if report.get("revision") not in {LEGACY_REVISION, REVISION} or report.get("status") != "accepted":
        raise InputGeometryError("unsupported or unsuccessful input repair evidence")
    if (_digest(original) != report["original_sha256"] or _digest(prepared) != report["prepared_sha256"]
            or request.normalized_pdb.sha256 != _digest(prepared)):
        raise InputGeometryError("input repair evidence source binding mismatch")
    options = preparation_options(request)
    if options is not None and (options["original_sha256"] != _digest(original)
                               or options["revision"] != report["revision"]):
        raise InputGeometryError("input repair original differs from request binding")
    before_request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=_digest(original)),
                             protein_context=None, gap_contexts=None)
    findings = assess_input_geometry(before_request, original)
    selected = _select_residues(before_request, original, findings,
                                allow_disulfides=report["revision"] == REVISION)
    if [asdict(r) for r in selected] != report["selected_residues"]:
        raise InputGeometryError("input repair mask differs from independently detected damage")
    before, after = _records(original), _records(prepared)
    if before.keys() != after.keys():
        raise InputGeometryError("input repair changed atom identities")
    selected_set = set(selected)
    # Check every physical record, including headers, CONECT, ANISOU and HETATM.
    old_lines, new_lines = original.splitlines(keepends=True), prepared.splitlines(keepends=True)
    if len(old_lines) != len(new_lines):
        raise InputGeometryError("input repair changed physical record count")
    for old, new in zip(old_lines, new_lines, strict=True):
        if old.startswith((b"ATOM  ", b"HETATM")):
            residue = ResidueIdentity(chr(old[21]), old[22:26].decode().strip(), old[26:27].decode().strip())
            if residue in selected_set:
                if old[:30] != new[:30] or old[54:] != new[54:]:
                    raise InputGeometryError("input repair changed observed atom metadata")
                continue
        if old != new:
            raise InputGeometryError("input repair changed a record outside the authorized coordinate mask")
    positions = _coordinates(after)
    if any(not all(math.isfinite(x) for x in xyz) for xyz in positions.values()):
        raise InputGeometryError("input repair produced nonfinite coordinates")
    restraints = build_anchor_restraints(original.decode(), selected, observed_atoms=request.fixed_atoms,
                                         base_k_kj_mol_nm2=RESTRAINT_K, bfactor_source="disabled")
    movements = check_anchor_movements(restraints, positions, AnchorMovementLimits())
    if not all(movement.within_limits for movement in movements):
        raise InputGeometryError("input repair exceeded observed-residue displacement bounds")
    _check_selected_links(request, prepared, selected, allow_disulfides=report["revision"] == REVISION)
    after_findings = assess_input_geometry(request, prepared)
    if any(f.severity is Severity.ERROR for f in after_findings):
        raise InputGeometryError("input repair left or introduced severe geometry defects")
    return movements, after_findings


def _without_preparation(request):
    return replace(request, backend_options=tuple(
        option for option in request.backend_options if not option.name.startswith(PREFIX)))


def audit_composed_preparation(request, original, prepared, report):
    """Independently audit original -> completed -> locally repaired observations.

    The completion-only intermediate is NOT sampling permission. Its records,
    added-atom inventory and L chirality are checked first; the second audit then
    requires clean final geometry and independently bounded observed movement.
    """
    from dvbfixer.model.diffusion.observed_completion import _insert, audit_completion_preparation

    if report.get("revision") != COMPOSED_REVISION or report.get("status") != "accepted":
        raise InputGeometryError("unsupported or unsuccessful composed preparation evidence")
    if (_digest(original) != report["original_sha256"] or _digest(prepared) != report["prepared_sha256"]
            or request.normalized_pdb.sha256 != _digest(prepared)):
        raise InputGeometryError("composed preparation source binding mismatch")
    options = preparation_options(request)
    if options is not None and (options["revision"] != COMPOSED_REVISION
                               or options["original_sha256"] != _digest(original)):
        raise InputGeometryError("composed preparation original binding mismatch")
    completion = report["completion"]
    completed = _insert(original, completion["added_records"])
    intermediate_request = replace(_without_preparation(request),
        normalized_pdb=ArtifactReference("completed.pdb", _digest(completed)), protein_context=None, gap_contexts=None)
    _, completion_findings = audit_completion_preparation(
        intermediate_request, original, completed, completion, require_final_geometry=False)
    if completion.get("after_findings") != [asdict(f) for f in completion_findings]:
        raise InputGeometryError("completion intermediate findings differ from independent assessment")
    movements, findings = audit_prepared_input(
        _without_preparation(request), completed, prepared, report["observed_repair"])
    if (report["observed_repair"].get("movements") != [asdict(m) for m in movements]
            or report["observed_repair"].get("after_findings") != [asdict(f) for f in findings]):
        raise InputGeometryError("observed repair evidence differs from independent assessment")
    # The movement audit above does not replace a global chirality assertion.
    from io import StringIO

    from openmm.app import PDBFile

    from dvbfixer.ffutils.geometry import assert_all_l

    pdb = PDBFile(StringIO(prepared.decode()))
    assert_all_l(pdb.topology, pdb.positions)
    return movements, findings


def complete_then_repair(request, original, completed, completion, directory):
    """Continue incomplete-source preparation instead of rejecting old clashes."""
    from dvbfixer.model.diffusion.observed_completion import audit_completion_preparation

    completed_path = directory / "completed.pdb"
    completed_path.write_bytes(completed)
    intermediate_request = replace(request,
        normalized_pdb=ArtifactReference(completed_path.name, _digest(completed)))
    _, completion_findings = audit_completion_preparation(
        intermediate_request, original, completed, completion, require_final_geometry=False)
    completion["after_findings"] = [asdict(f) for f in completion_findings]
    completion["movements"] = []
    (directory / "completion-evidence.json").write_text(
        json.dumps(completion, ensure_ascii=False, indent=2) + "\n")
    result = _prepare_observed_geometry(intermediate_request, completed, directory)
    if result is None:
        raise InputGeometryError("heavy-atom preparation left or introduced severe geometry defects")
    repaired_request, repaired_path = result
    repaired = repaired_path.read_bytes()
    report = {"revision": COMPOSED_REVISION, "status": "accepted",
              "original_sha256": _digest(original), "prepared_sha256": _digest(repaired),
              "completion": completion,
              "observed_repair": json.loads((repaired_path.parent / EVIDENCE_FILES["report"]).read_text())}
    path = directory / EVIDENCE_FILES["prepared"]
    path.write_bytes(repaired)
    new_request = replace(_without_preparation(repaired_request),
        normalized_pdb=ArtifactReference(path.name, _digest(repaired)))
    movements, findings = audit_composed_preparation(new_request, original, repaired, report)
    report["movements"] = [asdict(m) for m in movements]
    report["after_findings"] = [asdict(f) for f in findings]
    report_bytes = (json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode()
    (directory / EVIDENCE_FILES["report"]).write_bytes(report_bytes)
    new_request = replace(new_request, backend_options=(*new_request.backend_options,
        BackendOption(PREFIX + "revision", COMPOSED_REVISION),
        BackendOption(PREFIX + "original_sha256", _digest(original)),
        BackendOption(PREFIX + "report_sha256", _digest(report_bytes))))
    (directory / "request-prepared.json").write_text(new_request.to_json())
    return new_request, path


def audit_preparation_payloads(request, payloads):
    if set(payloads) != PUBLICATION_ROLES:
        raise InputGeometryError("input preparation publication is incomplete")
    options = preparation_options(request)
    if options is None or _digest(payloads["input-repair-report"]) != options["report_sha256"]:
        raise InputGeometryError("input repair report differs from request binding")
    if payloads["input-repair-request"] != request.to_json().encode():
        raise InputGeometryError("input repair snapshot differs from the scientific request")
    report = json.loads(payloads["input-repair-report"])
    movements, findings = audit_prepared_input(request, payloads["input-original"], payloads["input-prepared"], report)
    if report.get("movements") != [asdict(m) for m in movements] or report.get("after_findings") != [asdict(f) for f in findings]:
        raise InputGeometryError("recorded input repair audit differs from independent measurements")


def prepare_input(request: DiffusionRequest, source: bytes, work_parent: Path):
    """Return (prepared request, prepared path), or None if no repair is needed."""
    from dvbfixer.model.diffusion.observed_completion import (
        missing_observed_atoms,
        prepare_observed_completion,
    )

    missing = missing_observed_atoms(request, source)
    if missing:
        return prepare_observed_completion(request, source, missing, work_parent)
    return _prepare_observed_geometry(request, source, work_parent)


def _prepare_observed_geometry(request: DiffusionRequest, source: bytes, work_parent: Path):
    """Bounded repair of complete observed residues, also after atom completion."""
    findings = assess_input_geometry(request, source)
    if not any(f.severity is Severity.ERROR for f in findings):
        return None
    if request.protein_context is not None or request.gap_contexts is not None:
        raise InputGeometryError("automatic input repair must precede motif selection and freezing; frozen requests are never rewritten")
    directory = Path(tempfile.mkdtemp(prefix="dvbfixer-input-repair-", dir=work_parent))
    (directory / EVIDENCE_FILES["original"]).write_bytes(source)
    (directory / "request-before.json").write_text(request.to_json())
    (directory / "before-findings.json").write_text(json.dumps([asdict(f) for f in findings], ensure_ascii=False, indent=2))
    record_input_geometry(request, source, findings, directory, repair_planned=True)
    report = {"revision": REVISION, "status": "failed", "original_sha256": _digest(source),
              "policy": {"restraint_k_kj_mol_nm2": RESTRAINT_K, "max_iterations": MAX_ITERATIONS,
                         "movement_limits": asdict(AnchorMovementLimits()), "restart_count": 0,
                         "platform": "CPU", "solvent": False,
                         "temporary_gap_chemistry": "standard fragment termini; no gap-spanning peptide springs"}}
    try:
        selected = _select_residues(request, source, findings)
        report["selected_residues"] = [asdict(r) for r in selected]
        restraints = build_anchor_restraints(source.decode(), selected, observed_atoms=request.fixed_atoms,
                                             base_k_kj_mol_nm2=RESTRAINT_K, bfactor_source="disabled")
        ff_input = _fragmented_forcefield_input(request, source)
        (directory / "forcefield-input.pdb").write_text(ff_input)
        report["status"] = "running"
        (directory / EVIDENCE_FILES["report"]).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
        print(f"WARNING: Attempting bounded local input repair for {len(selected)} observed residues on CPU; "
              f"all other source records stay fixed. Evidence: {directory}", flush=True)
        from dvbfixer.model.diffusion.boundary_refinement import (
            BOUNDARY_REFINEMENT_FORCEFIELD,
            refine_observed_region,
        )
        from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates

        result = refine_observed_region(ff_input, residues=selected, atoms=tuple(r.atom for r in restraints),
                                        restraints=restraints, max_iterations=MAX_ITERATIONS)
        import openmm

        report["openmm_version"] = openmm.__version__
        report["forcefield"] = list(BOUNDARY_REFINEMENT_FORCEFIELD)
        prepared = rewrite_generated_coordinates(source.decode(), result.coordinates_angstrom).encode()
        path = directory / EVIDENCE_FILES["prepared"]
        path.write_bytes(prepared)
        new_request = replace(request, normalized_pdb=ArtifactReference(path.name, _digest(prepared)))
        report.update(status="accepted", prepared_sha256=_digest(prepared),
                      initial_energy_kj_mol=result.initial_energy_kj_mol, final_energy_kj_mol=result.final_energy_kj_mol,
                      chirality_repairs=[asdict(r) for r in result.chirality_repairs])
        movements, after = audit_prepared_input(new_request, source, prepared, report)
        report["movements"] = [asdict(m) for m in movements]
        report["after_findings"] = [asdict(f) for f in after]
        report_bytes = (json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode()
        (directory / EVIDENCE_FILES["report"]).write_bytes(report_bytes)
        options = (*request.backend_options, BackendOption(PREFIX + "revision", REVISION),
                   BackendOption(PREFIX + "original_sha256", _digest(source)),
                   BackendOption(PREFIX + "report_sha256", _digest(report_bytes)))
        new_request = replace(new_request, backend_options=options)
        (directory / "request-prepared.json").write_text(new_request.to_json())
        print(f"WARNING: Input repair accepted after independent geometry/movement checks; original preserved. "
              f"Diffusion will use {path}. Report: {directory / EVIDENCE_FILES['report']}", flush=True)
        return new_request, path
    except Exception as exc:
        report.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        (directory / EVIDENCE_FILES["report"]).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
        raise InputGeometryError(f"automatic input repair failed: {exc}. Sampling was not started. "
                                 f"Original and attempt evidence retained: {directory}") from exc
