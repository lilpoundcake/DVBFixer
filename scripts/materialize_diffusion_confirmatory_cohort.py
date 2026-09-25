#!/usr/bin/env python3
"""Screen and materialize the locked diffusion confirmatory cohort.

The stable output layout is::

    <root>/downloads/<PDB>.cif
    <root>/cases/<pdb-id>/workspace/{input/normalized.pdb,reference.pdb,
                                    target.fasta,request.json,source.json}
    <root>/cases/<pdb-id>/{accepted.json|rejected.json}
    <root>/cases/<pdb-id>/.complete.json
    <root>/report.json

Only an accepted case with an atomically published ``.complete.json`` marker is
runnable.  A batch runner can consume ``report.json`` and use the relative
``workspace`` value for every case whose status is ``accepted``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.scope import assess_diffusion_scope
from dvbfixer.structure_input import PDB_CHAIN_IDS, StructureInputError, normalize_structure

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_COHORT = REPO_ROOT / "docs/research/diffusion-confirmatory-cohort.json"
DEFAULT_OUTPUT = REPO_ROOT / ".artifacts/diffusion-confirmatory-cohort"
REPORT_SCHEMA_VERSION = 1
LAYOUT_VERSION = "dvbfixer-diffusion-confirmatory-workspace-v4"
MASK_SELECTION_SEED = "dvbfixer-diffusion-confirmatory-mask-v1"
CONTEXT_RADIUS_ANGSTROM = 5.0
CANONICAL = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
BACKBONE = frozenset({"N", "CA", "C", "O"})
_CASE_ID = re.compile(r"^[A-Za-z0-9]+$")


class ScreeningRejection(RuntimeError):
    """A stable scientific rejection, as opposed to a transient run failure."""

    def __init__(self, code: str, detail: str, **evidence: Any) -> None:
        super().__init__(detail)
        self.code = code
        self.detail = detail
        self.evidence = evidence


@dataclass(frozen=True)
class ObservedResidue:
    entity_index: int
    residue: Any


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode()


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _download_https(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": "dvbfixer-confirmatory-materializer/1"})
    with urlopen(request, timeout=120) as response:
        final_url = response.geturl()
        if urlsplit(final_url).scheme.lower() != "https":
            raise RuntimeError(f"download redirected to a non-HTTPS URL: {final_url}")
        return response.read()


def _source_cif(
    case: dict[str, Any],
    output_root: Path,
    fetch: Callable[[str], bytes],
) -> tuple[Path, str]:
    url = case.get("mmcif_url")
    if not isinstance(url, str) or urlsplit(url).scheme.lower() != "https":
        raise ScreeningRejection("non-https-source", "mmCIF source URL must use HTTPS")
    pdb_id = case.get("pdb_id")
    if not isinstance(pdb_id, str) or not _CASE_ID.fullmatch(pdb_id):
        raise ScreeningRejection("invalid-case-metadata", "pdb_id is not path-safe")
    parsed = urlsplit(url)
    expected_path = f"/download/{pdb_id.upper()}.cif"
    if (
        parsed.hostname != "files.rcsb.org"
        or parsed.path != expected_path
        or parsed.query
        or parsed.fragment
    ):
        raise ScreeningRejection(
            "unpinned-rcsb-source",
            f"mmCIF URL must be exactly https://files.rcsb.org{expected_path}",
        )
    destination = output_root / "downloads" / f"{pdb_id.upper()}.cif"
    metadata_path = destination.with_suffix(".source.json")
    if not destination.is_file() or not metadata_path.is_file():
        payload = fetch(url)
        if not payload:
            raise RuntimeError(f"empty response for {url}")
        _atomic_write(destination, payload)
        _atomic_write(
            metadata_path,
            _json_bytes({"url": url, "sha256": _sha256(payload), "bytes": len(payload)}),
        )
    payload = destination.read_bytes()
    digest = _sha256(payload)
    if not metadata_path.is_file():
        raise RuntimeError(f"download cache metadata is missing: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("url") != url or metadata.get("sha256") != digest:
        raise RuntimeError(f"download cache metadata does not verify {destination}")
    return destination, digest


def _required_case(case: dict[str, Any], key: str, expected: type[Any]) -> Any:
    value = case.get(key)
    if not isinstance(value, expected):
        raise ScreeningRejection("invalid-case-metadata", f"{key} has the wrong type")
    return value


def _load_structure(path: Path) -> Any:
    import gemmi

    try:
        structure = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Mmcif)
    except Exception as exc:
        raise ScreeningRejection("invalid-mmcif", f"Gemmi could not parse the mmCIF: {exc}") from exc
    if len(structure) != 1:
        raise ScreeningRejection(
            "multiple-models",
            f"expected exactly one coordinate model, found {len(structure)}",
            model_count=len(structure),
        )
    return structure


def _select_instance(structure: Any, case: dict[str, Any]) -> tuple[Any, list[ObservedResidue]]:
    label_chain = _required_case(case, "label_asym_id", str)
    author_chain = _required_case(case, "auth_asym_id", str)
    entity_id = _required_case(case, "entity_id", str)
    sequence = _required_case(case, "sequence", str)
    if len(author_chain) != 1 or author_chain not in PDB_CHAIN_IDS:
        raise ScreeningRejection(
            "unsupported-chain-identity",
            f"author chain {author_chain!r} is not a one-character PDB chain ID",
        )
    matches: list[Any] = []
    conflicting_entities: set[str] = set()
    for chain in structure[0]:
        for residue in chain:
            if chain.name == author_chain and residue.subchain == label_chain:
                if residue.entity_id != entity_id:
                    conflicting_entities.add(str(residue.entity_id))
                    continue
                matches.append(residue)
    if conflicting_entities:
        raise ScreeningRejection(
            "chain-entity-mismatch",
            "declared label/auth chain pair does not uniquely identify the declared entity",
            observed_entity_ids=sorted(conflicting_entities),
        )
    if not matches:
        raise ScreeningRejection(
            "chain-not-found",
            "declared case-sensitive label_asym_id/auth_asym_id pair was not found",
        )

    observed: list[ObservedResidue] = []
    seen_indices: set[int] = set()
    previous_index = -1
    for residue in matches:
        if residue.label_seq is None:
            continue
        if residue.name not in CANONICAL:
            raise ScreeningRejection(
                "noncanonical-residue",
                f"selected polymer contains noncanonical residue {residue.name!r}",
            )
        index = int(residue.label_seq) - 1
        if index < 0 or index >= len(sequence):
            raise ScreeningRejection(
                "sequence-placement-mismatch",
                f"label_seq_id {residue.label_seq} is outside the declared sequence",
            )
        if index in seen_indices:
            raise ScreeningRejection(
                "ambiguous-sequence-placement",
                f"multiple residues map to declared sequence index {index}",
            )
        if index <= previous_index:
            raise ScreeningRejection(
                "ambiguous-sequence-placement",
                "selected residues are not in strictly increasing label_seq_id order",
                previous_label_seq_id=previous_index + 1,
                label_seq_id=index + 1,
            )
        if CANONICAL[residue.name] != sequence[index]:
            raise ScreeningRejection(
                "sequence-placement-mismatch",
                f"coordinate residue {residue.name} disagrees with declared sequence index {index}",
            )
        seen_indices.add(index)
        previous_index = index
        observed.append(ObservedResidue(index, residue))
    if not observed:
        raise ScreeningRejection(
            "insufficient-contiguous-observed-residues",
            "selected entity has no canonical observed-coordinate residues",
        )
    return matches, observed


def _has_altloc(residue: Any) -> bool:
    return any(atom.altloc not in {"\x00", " ", "."} for atom in residue)


def _distance_squared(first: Any, second: Any) -> float:
    return (
        (first.x - second.x) ** 2
        + (first.y - second.y) ** 2
        + (first.z - second.z) ** 2
    )


def _connection_touches(connection: Any, author_chain: str, identities: set[tuple[int, str]]) -> bool:
    for partner in (connection.partner1, connection.partner2):
        seqid = partner.res_id.seqid
        if partner.chain_name == author_chain and (seqid.num, seqid.icode.strip()) in identities:
            return True
    return False


def _window_reasons(
    structure: Any,
    observed: list[ObservedResidue],
    start: int,
    stop: int,
    author_chain: str,
) -> tuple[str, ...]:
    local = observed[start - 1 : stop + 1]
    local_positions = [atom.pos for item in local for atom in item.residue]
    if any(_has_altloc(item.residue) for item in local):
        return ("alternate-location-ambiguity-near-mask",)

    radius_squared = CONTEXT_RADIUS_ANGSTROM ** 2
    reasons: list[str] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.het_flag == "A":
                    continue
                if any(
                    _distance_squared(atom.pos, local_position) <= radius_squared
                    for atom in residue
                    for local_position in local_positions
                ):
                    reasons.append("heteroatom-context-near-mask")
                    break
            if reasons:
                break
        if reasons:
            break
    identities = {
        (item.residue.seqid.num, item.residue.seqid.icode.strip())
        for item in local
    }
    if any(
        _connection_touches(connection, author_chain, identities)
        for connection in structure.connections
    ):
        reasons.append("covalent-link-context-near-mask")
    return tuple(dict.fromkeys(reasons))


def _select_mask(
    structure: Any,
    observed: list[ObservedResidue],
    case: dict[str, Any],
) -> tuple[int, int, dict[str, int]]:
    length = _required_case(case, "planned_gap_length", int)
    if length not in {5, 10}:
        raise ScreeningRejection(
            "unsupported-gap-length", f"planned gap length {length} is not 5 or 10"
        )
    candidates: list[tuple[str, int, int]] = []
    for start in range(1, len(observed) - length):
        stop = start + length
        local = observed[start - 1 : stop + 1]
        if len(local) != length + 2:
            continue
        # A natural entity-sequence gap may exist elsewhere, but the withheld
        # segment and its anchors must itself be a genuinely contiguous run.
        entity_indices = [item.entity_index for item in local]
        if entity_indices != list(range(entity_indices[0], entity_indices[0] + length + 2)):
            continue
        if any(not BACKBONE <= {atom.name for atom in item.residue} for item in local):
            continue
        token = (
            f"{MASK_SELECTION_SEED}:{case['pdb_id']}:{case['label_asym_id']}:"
            f"{case['auth_asym_id']}:{length}:{start}:{observed[start].entity_index}"
        )
        candidates.append((hashlib.sha256(token.encode()).hexdigest(), start, stop))
    if not candidates:
        raise ScreeningRejection(
            "insufficient-contiguous-observed-residues",
            f"no internal run supports a {length}-residue mask with two complete anchors",
        )
    rejected_reasons: Counter[str] = Counter()
    for _rank, start, stop in sorted(candidates):
        reasons = _window_reasons(structure, observed, start, stop, case["auth_asym_id"])
        if not reasons:
            return start, stop, dict(sorted(rejected_reasons.items()))
        rejected_reasons.update(reasons)
    raise ScreeningRejection(
        "no-clean-internal-mask",
        "every deterministic mask candidate had unsupported local context",
        candidate_count=len(candidates),
        local_rejection_counts=dict(sorted(rejected_reasons.items())),
    )


def _selected_structure(structure: Any, observed: list[ObservedResidue], author_chain: str) -> Any:
    import gemmi

    selected = gemmi.Structure()
    selected.name = structure.name or "dvbfixer_confirmatory_case"
    model = gemmi.Model(1)
    chain = gemmi.Chain(author_chain)
    for item in observed:
        chain.add_residue(item.residue.clone())
    model.add_chain(chain)
    selected.add_model(model)
    return selected


def _pdb_residue(line: str) -> ResidueIdentity:
    return ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())


def _pdb_atom(line: str) -> AtomIdentity:
    residue = _pdb_residue(line)
    return AtomIdentity(
        residue.chain,
        residue.residue_number,
        residue.insertion_code,
        line[12:16].strip(),
    )


def _is_heavy(line: str) -> bool:
    element = line[76:78].strip() if len(line) >= 78 else ""
    return (element or next((char for char in line[12:16] if char.isalpha()), "")).upper() != "H"


def _build_workspace(
    staging: Path,
    structure: Any,
    observed: list[ObservedResidue],
    case: dict[str, Any],
    start: int,
    stop: int,
    source: dict[str, Any],
    *,
    seed: int,
) -> dict[str, Any]:
    author_chain = case["auth_asym_id"]
    selected_cif = staging / "selected.cif"
    _selected_structure(structure, observed, author_chain).make_mmcif_document().write_file(
        str(selected_cif)
    )
    converted = staging / "selected.pdb"
    try:
        normalized = normalize_structure(selected_cif, converted)
    except StructureInputError as exc:
        message = str(exc)
        code = (
            "pdb-fixed-column-overflow"
            if any(
                phrase in message
                for phrase in (
                    "outside the PDB",
                    "PDB atom-serial limit",
                    "PDB supports at most",
                    "exceeds 3 PDB",
                    "exceeds 4 PDB",
                    "not representable",
                )
            )
            else "cif-normalization-failed"
        )
        raise ScreeningRejection(code, message) from exc
    if normalized.chain_map != {author_chain: author_chain}:
        raise ScreeningRejection(
            "unsupported-chain-mapping",
            "selected author chain was not preserved exactly by the CIF boundary",
            chain_map=normalized.chain_map,
        )

    atom_lines = [
        line
        for line in converted.read_text(encoding="utf-8").splitlines(keepends=True)
        if line.startswith("ATOM  ") and line[21] == author_chain and _is_heavy(line)
    ]
    if not atom_lines:
        raise ScreeningRejection("cif-normalization-failed", "selected chain produced no ATOM records")
    reference_bytes = ("".join(atom_lines) + "TER\nEND\n").encode()
    ordered: list[tuple[ResidueIdentity, str]] = []
    seen: set[ResidueIdentity] = set()
    for line in atom_lines:
        identity = _pdb_residue(line)
        if identity not in seen:
            seen.add(identity)
            ordered.append((identity, line[17:20].strip()))
    if len(ordered) != len(observed):
        raise ScreeningRejection(
            "identity-preservation-failed",
            "CIF boundary changed the selected residue count",
        )
    sequence = "".join(CANONICAL[name] for _identity, name in ordered)
    expected_observed_sequence = "".join(
        CANONICAL[item.residue.name] for item in observed
    )
    if sequence != expected_observed_sequence:
        raise ScreeningRejection(
            "identity-preservation-failed",
            "PDB materialization changed the observed-coordinate residue sequence",
        )

    generated_residues = tuple(identity for identity, _name in ordered[start:stop])
    generated_set = set(generated_residues)
    source_bytes = "".join(
        line for line in atom_lines if _pdb_residue(line) not in generated_set
    ).encode() + b"TER\nEND\n"
    generated_atoms = tuple(
        _pdb_atom(line)
        for line in atom_lines
        if _pdb_residue(line) in generated_set and _is_heavy(line)
    )
    fixed_atoms = tuple(
        _pdb_atom(line)
        for line in source_bytes.decode().splitlines()
        if line.startswith("ATOM  ") and _is_heavy(line)
    )
    observed_indices = tuple(index for index in range(len(ordered)) if not start <= index < stop)
    request = DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("input/normalized.pdb", _sha256(source_bytes)),
        target_sequences=(TargetSequence(author_chain, sequence),),
        sequence_placements=(
            SequencePlacement(
                chain=author_chain,
                target_length=len(sequence),
                observed_target_indices=observed_indices,
                observed_residues=tuple(ordered[index][0] for index in observed_indices),
            ),
        ),
        gaps=(
            GapRegion(
                chain=author_chain,
                target_interval=TargetInterval(start, stop),
                left_anchor=ordered[start - 1][0],
                right_anchor=ordered[stop][0],
                generated_residues=generated_residues,
                movable_junction_residues=(
                    ordered[start - 1][0],
                    *generated_residues,
                    ordered[stop][0],
                ),
            ),
        ),
        fixed_atoms=fixed_atoms,
        generated_atoms=generated_atoms,
        retained_explicit_links=(),
        candidate_count=1,
        seeds=(seed,),
    )
    admission = assess_diffusion_scope(request, source_bytes)
    if not admission.supported:
        raise ScreeningRejection(
            "unsupported-diffusion-scope",
            "materialized request failed the existing diffusion scope gate",
            scope_reasons=list(admission.reasons),
        )

    workspace = staging / "workspace"
    (workspace / "input").mkdir(parents=True)
    artifacts = {
        "input/normalized.pdb": source_bytes,
        "reference.pdb": reference_bytes,
        "target.fasta": f">chain_{author_chain}\n{sequence}\n".encode(),
        "request.json": request.to_json().encode(),
        "source.json": _json_bytes(source),
    }
    for relative, payload in artifacts.items():
        destination = workspace / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
    selected_cif.unlink()
    converted.unlink()
    return {
        "workspace": "workspace",
        "artifacts": {
            relative: {"sha256": _sha256(payload), "bytes": len(payload)}
            for relative, payload in sorted(artifacts.items())
        },
        "mask": {
            "target_interval_zero_based_half_open": [start, stop],
            "entity_interval_zero_based_half_open": [
                observed[start].entity_index,
                observed[stop - 1].entity_index + 1,
            ],
            "length": stop - start,
            "sequence": sequence[start:stop],
            "left_anchor": {
                "chain": request.gaps[0].left_anchor.chain,
                "residue_number": request.gaps[0].left_anchor.residue_number,
                "insertion_code": request.gaps[0].left_anchor.insertion_code,
            },
            "right_anchor": {
                "chain": request.gaps[0].right_anchor.chain,
                "residue_number": request.gaps[0].right_anchor.residue_number,
                "insertion_code": request.gaps[0].right_anchor.insertion_code,
            },
            "generated_residues": [
                {
                    "chain": item.chain,
                    "residue_number": item.residue_number,
                    "insertion_code": item.insertion_code,
                }
                for item in generated_residues
            ],
        },
    }


def _case_id(case: dict[str, Any]) -> str:
    pdb_id = _required_case(case, "pdb_id", str)
    if not _CASE_ID.fullmatch(pdb_id):
        raise ScreeningRejection("invalid-case-metadata", "pdb_id is not path-safe")
    return pdb_id.lower()


def _valid_completion(case_dir: Path, cohort_sha256: str) -> dict[str, Any] | None:
    marker = case_dir / ".complete.json"
    if not marker.is_file():
        return None
    try:
        value = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"invalid completion marker {marker}: {exc}") from exc
    if value.get("cohort_sha256") != cohort_sha256 or value.get("layout_version") != LAYOUT_VERSION:
        raise RuntimeError(f"stale completion marker requires a new output root: {marker}")
    record_name = "accepted.json" if value.get("status") == "accepted" else "rejected.json"
    record_path = case_dir / record_name
    if not record_path.is_file() or _sha256(record_path.read_bytes()) != value.get("record_sha256"):
        raise RuntimeError(f"completion marker does not verify its case record: {marker}")
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if record.get("status") == "accepted":
        if not isinstance(record.get("selected_for_inference"), bool):
            raise RuntimeError(f"accepted record lacks an inference-selection decision: {record_path}")
        for relative, metadata in record.get("artifacts", {}).items():
            artifact = case_dir / "workspace" / relative
            if not artifact.is_file() or _sha256(artifact.read_bytes()) != metadata.get("sha256"):
                raise RuntimeError(f"completed artifact failed digest verification: {artifact}")
    return record


def _publish_case(case_dir: Path, staging: Path, record: dict[str, Any], cohort_sha256: str) -> None:
    record_name = "accepted.json" if record["status"] == "accepted" else "rejected.json"
    record_bytes = _json_bytes(record)
    (staging / record_name).write_bytes(record_bytes)
    marker = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "layout_version": LAYOUT_VERSION,
        "cohort_sha256": cohort_sha256,
        "case_id": record["case_id"],
        "status": record["status"],
        "record": record_name,
        "record_sha256": _sha256(record_bytes),
    }
    (staging / ".complete.json").write_bytes(_json_bytes(marker))
    if case_dir.exists():
        shutil.rmtree(case_dir)
    os.replace(staging, case_dir)


def _rewrite_accepted_record(
    case_dir: Path,
    record: dict[str, Any],
    cohort_sha256: str,
) -> None:
    """Publish an inference-selection decision, with the marker written last."""
    record_bytes = _json_bytes(record)
    _atomic_write(case_dir / "accepted.json", record_bytes)
    marker = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "layout_version": LAYOUT_VERSION,
        "cohort_sha256": cohort_sha256,
        "case_id": record["case_id"],
        "status": "accepted",
        "record": "accepted.json",
        "record_sha256": _sha256(record_bytes),
    }
    _atomic_write(case_dir / ".complete.json", _json_bytes(marker))


def _finalize_inference_selection(
    records: list[dict[str, Any]],
    output_root: Path,
    cohort_sha256: str,
    *,
    expected_count: int,
    final_maximum: int,
) -> list[dict[str, Any]]:
    counts = Counter(record["status"] for record in records)
    if (
        len(records) != expected_count
        or counts["pending"]
        or counts["operational-error"]
    ):
        return records

    selected_count = 0
    finalized: list[dict[str, Any]] = []
    for record in records:
        if record["status"] != "accepted":
            finalized.append(record)
            continue
        selected = selected_count < final_maximum
        if selected:
            selected_count += 1
        selection_status = "selected-for-inference" if selected else "reserve-not-selected"
        updated = {
            **record,
            "selected_for_inference": selected,
            "selection_status": selection_status,
        }
        if updated != record:
            _rewrite_accepted_record(
                output_root / "cases" / record["case_id"],
                updated,
                cohort_sha256,
            )
        finalized.append(updated)
    return finalized


def _screen_case(
    case: dict[str, Any],
    output_root: Path,
    cohort_sha256: str,
    fetch: Callable[[str], bytes],
    *,
    seed: int,
    leakage_resolved: bool,
) -> dict[str, Any]:
    case_id = _case_id(case)
    case_dir = output_root / "cases" / case_id
    completed = _valid_completion(case_dir, cohort_sha256)
    if completed is not None:
        return completed
    case_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{case_id}.", dir=case_dir.parent))
    source: dict[str, Any] | None = None
    try:
        cif_path, cif_sha256 = _source_cif(case, output_root, fetch)
        source = {
            "url": case["mmcif_url"],
            "download": str(cif_path.relative_to(output_root)),
            "sha256": cif_sha256,
            "entity_sequence_sha256": _sha256(case["sequence"].encode()),
            "entity_sequence_length": len(case["sequence"]),
            "entity_sequence_role": "rcsb-clustering-and-provenance-only",
            "request_sequence_basis": "observed-coordinate-sequence",
            "pdb_id": case["pdb_id"],
            "entity_id": case["entity_id"],
            "label_asym_id": case["label_asym_id"],
            "auth_asym_id": case["auth_asym_id"],
            "independence_group": case["independence_group"],
        }
        structure = _load_structure(cif_path)
        _matches, observed = _select_instance(structure, case)
        observed_sequence = "".join(CANONICAL[item.residue.name] for item in observed)
        source.update(
            {
                "observed_coordinate_sequence_sha256": _sha256(observed_sequence.encode()),
                "observed_coordinate_residue_count": len(observed),
                "observed_entity_indices_zero_based": [
                    item.entity_index for item in observed
                ],
            }
        )
        if any(_has_altloc(item.residue) for item in observed):
            raise ScreeningRejection(
                "alternate-location-ambiguity-in-target-chain",
                "target-chain altlocs cannot be represented without choosing a conformer",
            )
        start, stop, skipped_contexts = _select_mask(structure, observed, case)
        workspace = _build_workspace(
            staging,
            structure,
            observed,
            case,
            start,
            stop,
            source,
            seed=seed,
        )
        record = {
            "schema_version": REPORT_SCHEMA_VERSION,
            "layout_version": LAYOUT_VERSION,
            "case_id": case_id,
            "screening_index": case["screening_index"],
            "status": "accepted",
            "selected_for_inference": False,
            "selection_status": "awaiting-full-screen",
            "independence_group": case["independence_group"],
            "leakage_resolved": leakage_resolved,
            "source": source,
            "deterministic_candidates_skipped_by_reason": skipped_contexts,
            **workspace,
        }
        record["workspace"] = str(Path("cases") / case_id / "workspace")
    except ScreeningRejection as exc:
        record = {
            "schema_version": REPORT_SCHEMA_VERSION,
            "layout_version": LAYOUT_VERSION,
            "case_id": case_id,
            "screening_index": case["screening_index"],
            "status": "rejected",
            "source": source,
            "reason_codes": [exc.code],
            "detail": exc.detail,
            "evidence": exc.evidence,
        }
        workspace = staging / "workspace"
        if workspace.exists():
            shutil.rmtree(workspace)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    _publish_case(case_dir, staging, record, cohort_sha256)
    return record


def _report(
    cohort_path: Path,
    cohort_sha256: str,
    records: Iterable[dict[str, Any]],
    *,
    expected_count: int,
    final_minimum: int,
    final_maximum: int,
) -> dict[str, Any]:
    cases = list(records)
    counts = Counter(record["status"] for record in cases)
    screening_complete = (
        len(cases) == expected_count
        and counts["pending"] == 0
        and counts["operational-error"] == 0
    )
    selected_for_inference = sum(
        record.get("selected_for_inference") is True
        for record in cases
        if record["status"] == "accepted"
    )
    reserve_accepted = sum(
        record.get("selection_status") == "reserve-not-selected"
        for record in cases
        if record["status"] == "accepted"
    )
    accepted_count_sufficient = counts["accepted"] >= final_minimum
    final_count_valid = final_minimum <= selected_for_inference <= final_maximum
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "layout_version": LAYOUT_VERSION,
        "cohort": str(cohort_path),
        "cohort_sha256": cohort_sha256,
        "context_radius_angstrom": CONTEXT_RADIUS_ANGSTROM,
        "mask_selection_seed": MASK_SELECTION_SEED,
        "summary": {
            "total": len(cases),
            "accepted": counts["accepted"],
            "rejected": counts["rejected"],
            "operational_errors": counts["operational-error"],
            "pending": counts["pending"],
            "selected_for_inference": selected_for_inference,
            "reserve_accepted": reserve_accepted,
            "final_minimum": final_minimum,
            "final_maximum": final_maximum,
            "accepted_count_sufficient": accepted_count_sufficient,
            "screening_complete": screening_complete,
            "confirmatory_sample_valid": screening_complete and final_count_valid,
        },
        "cases": cases,
    }


def materialize_cohort(
    cohort_path: Path,
    output_root: Path,
    *,
    seed: int = 7,
    case_ids: set[str] | None = None,
    fetch: Callable[[str], bytes] = _download_https,
) -> dict[str, Any]:
    """Materialize selected cohort cases and atomically update ``report.json``."""
    cohort_bytes = cohort_path.read_bytes()
    cohort_sha256 = _sha256(cohort_bytes)
    cohort = json.loads(cohort_bytes)
    if cohort.get("schema_version") != 1 or not isinstance(cohort.get("cases"), list):
        raise ValueError("unsupported confirmatory cohort schema")
    cases = cohort["cases"]
    expected_count = cohort.get("screening_pool_count")
    final_minimum = cohort.get("final_minimum_independence_groups")
    final_maximum = cohort.get("final_maximum_independence_groups")
    if (
        cohort.get("request_sequence_basis") != "observed-coordinate-sequence"
        or cohort.get("screening_order_locked") is not True
        or not isinstance(expected_count, int)
        or expected_count != len(cases)
        or not isinstance(final_minimum, int)
        or not isinstance(final_maximum, int)
        or final_minimum <= 0
        or final_maximum < final_minimum
        or final_maximum > expected_count
    ):
        raise ValueError("confirmatory cohort does not declare the frozen screening scope")
    for index, case in enumerate(cases):
        if (
            case.get("screening_index") != index
            or case.get("request_sequence_basis") != "observed-coordinate-sequence"
        ):
            raise ValueError("confirmatory cohort screening order is not locked")
    leakage_resolved = "protenix-v1" in cohort.get("eligible_backends", [])
    output_root = output_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    seen_case_ids: set[str] = set()
    for case in cases:
        case_id = _case_id(case)
        if case_id in seen_case_ids:
            raise ValueError(f"cohort repeats case ID {case_id!r}")
        seen_case_ids.add(case_id)
        if case_ids is not None and case_id not in case_ids:
            completed = _valid_completion(output_root / "cases" / case_id, cohort_sha256)
            record = completed or {
                "schema_version": REPORT_SCHEMA_VERSION,
                "layout_version": LAYOUT_VERSION,
                "case_id": case_id,
                "screening_index": case["screening_index"],
                "status": "pending",
            }
        else:
            try:
                record = _screen_case(
                    case,
                    output_root,
                    cohort_sha256,
                    fetch,
                    seed=seed,
                    leakage_resolved=leakage_resolved,
                )
            except Exception as exc:
                record = {
                    "schema_version": REPORT_SCHEMA_VERSION,
                    "layout_version": LAYOUT_VERSION,
                    "case_id": case_id,
                    "screening_index": case["screening_index"],
                    "status": "operational-error",
                    "detail": f"{type(exc).__name__}: {exc}",
                }
        records.append(record)
        _atomic_write(
            output_root / "report.json",
            _json_bytes(
                _report(
                    cohort_path,
                    cohort_sha256,
                    records,
                    expected_count=expected_count,
                    final_minimum=final_minimum,
                    final_maximum=final_maximum,
                )
            ),
        )
    records = _finalize_inference_selection(
        records,
        output_root,
        cohort_sha256,
        expected_count=expected_count,
        final_maximum=final_maximum,
    )
    report = _report(
        cohort_path,
        cohort_sha256,
        records,
        expected_count=expected_count,
        final_minimum=final_minimum,
        final_maximum=final_maximum,
    )
    _atomic_write(output_root / "report.json", _json_bytes(report))
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument(
        "--case-id",
        action="append",
        default=[],
        help="materialize only this lowercase PDB case ID; repeatable",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = materialize_cohort(
        args.cohort,
        args.output_root,
        seed=args.seed,
        case_ids=set(args.case_id) or None,
    )
    print(args.output_root / "report.json")
    return 1 if report["summary"]["operational_errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
