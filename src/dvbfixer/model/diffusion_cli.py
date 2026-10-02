"""Narrow public orchestration for ``model --backend diffusion``."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

from dvbfixer.model.cli import AA3TO1
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    BackendOption,
    DiffusionRequest,
    ExplicitLink,
    GapRegion,
    ResidueIdentity,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.masks import (
    DiffusionMaskError,
    ObservedResidue,
    build_sequence_placement,
)
from dvbfixer.model.diffusion.pipeline import run_diffusion_pipeline
from dvbfixer.model.diffusion.runner import RunnerLimits
from dvbfixer.model.diffusion.scope import (
    MAXIMUM_GAP_LENGTH,
    MINIMUM_GAP_LENGTH,
    _parse_coordinate_records,
    _parse_explicit_links,
    canonical_heavy_atom_identities,
)


class DiffusionCliError(ValueError):
    """Raised when a CLI input cannot form the initial diffusion request."""


def build_cli_diffusion_request(
    input_path: Path,
    target_sequences: dict[str, str],
    *,
    seeds: tuple[int, ...],
    profile: str,
) -> DiffusionRequest:
    """Build the initial one-chain, one-gap canonical-protein request."""
    source = input_path.read_bytes()
    try:
        text = source.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise DiffusionCliError("normalized PDB is not valid UTF-8") from exc
    lines = text.splitlines()
    atoms, residue_names, _altlocs = _parse_coordinate_records(lines)

    protein_residues: list[ObservedResidue] = []
    seen: set[ResidueIdentity] = set()
    atom_names: dict[ResidueIdentity, list[str]] = {}
    for record in atoms:
        if record.record_name != "ATOM":
            continue
        residue = ResidueIdentity(
            record.identity.chain,
            record.identity.residue_number,
            record.identity.insertion_code,
        )
        atom_names.setdefault(residue, []).append(record.identity.atom_name)
        if residue in seen:
            continue
        seen.add(residue)
        one_letter = AA3TO1.get(record.residue_name)
        if one_letter is None:
            raise DiffusionCliError(
                f"non-canonical protein residue is unsupported: "
                f"{residue.chain}/{record.residue_name}{residue.residue_number}{residue.insertion_code}"
            )
        protein_residues.append(ObservedResidue(residue, one_letter, ("CA",)))

    chains = tuple(dict.fromkeys(residue.identity.chain for residue in protein_residues))
    if len(chains) != 1:
        raise DiffusionCliError("diffusion currently requires exactly one protein chain")
    chain = chains[0]
    if set(target_sequences) != {chain}:
        raise DiffusionCliError(
            "diffusion target sequence must cover exactly the single input protein chain"
        )
    observed = tuple(
        ObservedResidue(
            residue.identity,
            residue.one_letter_code,
            tuple(dict.fromkeys(atom_names[residue.identity])),
        )
        for residue in protein_residues
    )
    target = target_sequences[chain].upper()
    try:
        placement = build_sequence_placement(chain, observed, target)
    except DiffusionMaskError as exc:
        raise DiffusionCliError(str(exc)) from exc

    gaps = [
        (left + 1, right)
        for left, right in zip(
            placement.observed_target_indices,
            placement.observed_target_indices[1:],
        )
        if right > left + 1
    ]
    if placement.observed_target_indices[0] != 0 or placement.observed_target_indices[-1] != len(target) - 1:
        raise DiffusionCliError("terminal gaps are unsupported by diffusion")
    if len(gaps) != 1:
        raise DiffusionCliError("diffusion currently requires exactly one internal gap")
    start, stop = gaps[0]
    gap_length = stop - start
    if not MINIMUM_GAP_LENGTH <= gap_length <= MAXIMUM_GAP_LENGTH:
        raise DiffusionCliError(
            f"gap length {gap_length} is outside supported range "
            f"{MINIMUM_GAP_LENGTH}-{MAXIMUM_GAP_LENGTH}"
        )

    by_target = dict(zip(placement.observed_target_indices, placement.observed_residues))
    left_anchor = by_target[start - 1]
    right_anchor = by_target[stop]
    try:
        left_number = int(left_anchor.residue_number)
        right_number = int(right_anchor.residue_number)
    except ValueError as exc:
        raise DiffusionCliError("diffusion requires integer PDB residue numbers") from exc
    if right_number - left_number - 1 < gap_length:
        raise DiffusionCliError(
            "input numbering does not reserve enough residue numbers for the gap"
        )
    generated_residues = tuple(
        ResidueIdentity(chain, str(left_number + offset))
        for offset in range(1, gap_length + 1)
    )
    if set(generated_residues) & set(placement.observed_residues):
        raise DiffusionCliError("generated residue numbering collides with observed residues")
    gap = GapRegion(
        chain=chain,
        target_interval=TargetInterval(start, stop),
        left_anchor=left_anchor,
        right_anchor=right_anchor,
        generated_residues=generated_residues,
        movable_junction_residues=(left_anchor, *generated_residues, right_anchor),
    )

    fixed_atoms = tuple(
        record.identity for record in atoms if record.element.upper() != "H"
    )
    generated_atoms = tuple(
        atom
        for target_index, residue in zip(range(start, stop), generated_residues)
        for atom in canonical_heavy_atom_identities(
            residue,
            _one_to_three(target[target_index]),
        )
    )
    links = tuple(
        ExplicitLink(first, second, "PDB")
        for first, second in sorted(_parse_explicit_links(lines, atoms))
    )
    return DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference(input_path.name, hashlib.sha256(source).hexdigest()),
        target_sequences=(TargetSequence(chain, target),),
        sequence_placements=(placement,),
        gaps=(gap,),
        fixed_atoms=fixed_atoms,
        generated_atoms=generated_atoms,
        retained_explicit_links=links,
        candidate_count=len(seeds),
        seeds=seeds,
        backend_options=(BackendOption("profile", profile),),
    )


def run_diffusion_model(args: argparse.Namespace) -> None:
    input_path = Path(args.input).resolve()
    if not input_path.is_file():
        raise DiffusionCliError(f"input file does not exist: {input_path}")
    output_value = args.output
    output_path = (
        Path(output_value).resolve()
        if output_value
        else input_path.with_name(input_path.stem + "_model_diffusion")
    )
    if output_path.exists() or output_path.is_symlink():
        raise DiffusionCliError(f"diffusion output directory already exists: {output_path}")

    checkpoint = Path(args.diffusion_checkpoint).expanduser().resolve()
    if not checkpoint.is_file():
        raise DiffusionCliError(f"diffusion checkpoint does not exist: {checkpoint}")
    expected_digest = args.diffusion_checkpoint_sha256
    if expected_digest:
        actual_digest = _sha256(checkpoint)
        if actual_digest.lower() != expected_digest.lower():
            raise DiffusionCliError("diffusion checkpoint SHA-256 does not match")

    target_sequences = _target_sequences(input_path, args.fasta)
    request = build_cli_diffusion_request(
        input_path,
        target_sequences,
        seeds=tuple(args.diffusion_seeds),
        profile=args.diffusion_profile,
    )
    work_parent_value = args.diffusion_work_parent
    work_parent = (
        Path(work_parent_value).expanduser().resolve()
        if work_parent_value
        else output_path.parent
    )
    if not work_parent.is_dir():
        raise DiffusionCliError(f"diffusion work parent is not a directory: {work_parent}")
    command = (
        os.fspath(Path(args.diffusion_runner).expanduser()),
        "--profile",
        args.diffusion_profile,
        "--checkpoint",
        os.fspath(checkpoint),
    )
    outcome = run_diffusion_pipeline(
        request,
        command,
        source_root=input_path.parent,
        work_parent=work_parent,
        destination_bundle=output_path,
        limits=RunnerLimits(timeout_seconds=args.diffusion_timeout),
    )
    if outcome.published_bundle is None:
        raise DiffusionCliError(outcome.message or f"diffusion ended with {outcome.status.value}")
    print(f"Wrote diffusion bundle {outcome.published_bundle}")


def _target_sequences(input_path: Path, fasta_path: str | None) -> dict[str, str]:
    from dvbfixer.model.pipeline import parse_fasta, parse_seqres

    if fasta_path:
        return {chain: sequence.upper() for chain, sequence in parse_fasta(fasta_path).items()}
    seqres = parse_seqres(input_path.read_text(encoding="utf-8").splitlines())
    if not seqres:
        raise DiffusionCliError("diffusion requires SEQRES records or --fasta")
    return {
        chain: "".join(AA3TO1.get(residue, "X") for residue in residues)
        for chain, residues in seqres.items()
    }


def _one_to_three(one_letter: str) -> str:
    matches = [name for name, letter in AA3TO1.items() if letter == one_letter and name in {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    }]
    if not matches:
        raise DiffusionCliError(f"non-canonical target residue is unsupported: {one_letter}")
    return matches[0]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()
