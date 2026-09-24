#!/usr/bin/env python3
"""Build a contained withheld-coordinate diffusion benchmark workspace."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import tomllib
from pathlib import Path
from typing import Any

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

REPO_ROOT = Path(__file__).resolve().parent.parent
INVENTORY = REPO_ROOT / "docs/research/diffusion-gap-reconstruction-inventory.toml"
ONE_LETTER = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _residue(line: str) -> ResidueIdentity:
    return ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())


def _atom(line: str) -> AtomIdentity:
    residue = _residue(line)
    return AtomIdentity(
        residue.chain,
        residue.residue_number,
        residue.insertion_code,
        line[12:16].strip(),
    )


def _is_heavy(line: str) -> bool:
    element = line[76:78].strip() or line[12:16].strip()[:1]
    return element.upper() != "H"


def _identity_label(identity: ResidueIdentity, residue_name: str) -> str:
    insertion_code = identity.insertion_code or "_"
    return (
        f"{identity.chain}:{identity.residue_number}:"
        f"{insertion_code}:{residue_name}"
    )


def _fasta_sequence(path: Path, record_id: str) -> str:
    records: dict[str, list[str]] = {}
    current = ""
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split(maxsplit=1)[0]
            if current in records:
                raise ValueError(f"duplicate FASTA record {current!r}")
            records[current] = []
        elif current:
            records[current].append(line)
        else:
            raise ValueError("FASTA sequence appears before its header")
    if record_id not in records:
        raise ValueError(f"FASTA record {record_id!r} is absent")
    return "".join(records[record_id])


def _case(case_id: str) -> dict[str, Any]:
    with INVENTORY.open("rb") as handle:
        inventory = tomllib.load(handle)
    matches = [case for case in inventory["benchmark_cases"] if case["id"] == case_id]
    if len(matches) != 1:
        raise ValueError(f"benchmark case {case_id!r} does not exist uniquely")
    case = matches[0]
    if case["classification"] != "supported":
        raise ValueError(f"benchmark case {case_id!r} is not supported")
    if case["request_sequence_basis"] != "observed-coordinate-sequence":
        raise ValueError(f"benchmark case {case_id!r} has an unsupported sequence basis")
    return case


def build_workspace(case_id: str, destination: Path, *, seed: int) -> Path:
    """Build one new benchmark workspace and return its request manifest."""
    case = _case(case_id)
    fixture = REPO_ROOT / case["fixture"]
    fixture_bytes = fixture.read_bytes()
    if _sha256(fixture_bytes) != case["fixture_sha256"]:
        raise ValueError(f"fixture digest mismatch for {case_id}")
    sequence_fixture = REPO_ROOT / case["sequence_fixture"]
    sequence_bytes = sequence_fixture.read_bytes()
    if _sha256(sequence_bytes) != case["sequence_fixture_sha256"]:
        raise ValueError(f"sequence fixture digest mismatch for {case_id}")
    reference_sequence = _fasta_sequence(sequence_fixture, case["sequence_record"])
    target_start, target_stop = case["target_interval_zero_based_half_open"]
    if reference_sequence[target_start:target_stop] != case["expected_sequence"]:
        raise ValueError(f"reference sequence interval drifted for {case_id}")
    fixture_text = fixture_bytes.decode("utf-8")
    lines = fixture_text.splitlines(keepends=True)
    chain = case["target_chain"]
    structure_scope = case.get("request_structure_scope", "full-fixture")
    if structure_scope == "full-fixture":
        benchmark_lines = lines
        reference_bytes = fixture_bytes
    elif structure_scope == "target-protein-chain-only":
        benchmark_lines = [
            line
            for line in lines
            if line.startswith("ATOM  ") and line[21] == chain
        ]
        if not benchmark_lines:
            raise ValueError(f"target chain {chain!r} has no protein atoms")
        benchmark_lines.extend(("TER\n", "END\n"))
        reference_bytes = "".join(benchmark_lines).encode("utf-8")
    else:
        raise ValueError(f"unsupported request structure scope {structure_scope!r}")

    ordered: list[tuple[ResidueIdentity, str]] = []
    seen: set[ResidueIdentity] = set()
    for line in benchmark_lines:
        if not line.startswith("ATOM  ") or line[21] != chain:
            continue
        identity = _residue(line)
        if identity not in seen:
            seen.add(identity)
            ordered.append((identity, line[17:20].strip()))

    start, stop = case["fixture_interval_zero_based_half_open"]
    generated_rows = ordered[start:stop]
    if not generated_rows or start == 0 or stop >= len(ordered):
        raise ValueError(f"benchmark case {case_id!r} is not an internal mask")
    generated_residues = tuple(identity for identity, _name in generated_rows)
    generated_set = set(generated_residues)
    labels = [_identity_label(identity, name) for identity, name in generated_rows]
    if labels != case["generated_residues"]:
        raise ValueError(f"generated residue identities drifted for {case_id}")
    if _identity_label(*ordered[start - 1]) != case["left_anchor"]:
        raise ValueError(f"left anchor identity drifted for {case_id}")
    if _identity_label(*ordered[stop]) != case["right_anchor"]:
        raise ValueError(f"right anchor identity drifted for {case_id}")
    sequence = "".join(ONE_LETTER[name] for _identity, name in ordered)
    if sequence[start:stop] != case["expected_sequence"]:
        raise ValueError(f"generated sequence drifted for {case_id}")

    source_text = "".join(
        line
        for line in benchmark_lines
        if not (
            line.startswith(("ATOM  ", "HETATM"))
            and _residue(line) in generated_set
        )
    )
    source_bytes = source_text.encode("utf-8")
    generated_atoms = tuple(
        _atom(line)
        for line in benchmark_lines
        if line.startswith(("ATOM  ", "HETATM"))
        and _residue(line) in generated_set
        and _is_heavy(line)
    )
    fixed_atoms = tuple(
        _atom(line)
        for line in source_text.splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and _is_heavy(line)
    )
    observed_indices = tuple(
        index for index in range(len(ordered)) if not start <= index < stop
    )
    request = DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference(
            "input/normalized.pdb",
            _sha256(source_bytes),
        ),
        target_sequences=(TargetSequence(chain, sequence),),
        sequence_placements=(
            SequencePlacement(
                chain=chain,
                target_length=len(sequence),
                observed_target_indices=observed_indices,
                observed_residues=tuple(ordered[index][0] for index in observed_indices),
            ),
        ),
        gaps=(
            GapRegion(
                chain=chain,
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

    output = destination.expanduser().absolute()
    if output.exists():
        raise FileExistsError(f"destination already exists: {output}")
    (output / "input").mkdir(parents=True)
    (output / "input/normalized.pdb").write_bytes(source_bytes)
    (output / "request.json").write_text(request.to_json(), encoding="utf-8")
    (output / "target.fasta").write_text(
        f">chain_{chain}\n{sequence}\n",
        encoding="utf-8",
    )
    if structure_scope == "full-fixture":
        shutil.copy2(fixture, output / "reference.pdb")
    else:
        (output / "reference.pdb").write_bytes(reference_bytes)
    return output / "request.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_id")
    parser.add_argument("destination", type=Path)
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(build_workspace(args.case_id, args.destination, seed=args.seed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
