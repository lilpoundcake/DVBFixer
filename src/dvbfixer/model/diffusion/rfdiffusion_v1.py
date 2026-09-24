"""Internal RFdiffusion v1 backbone-inpainting benchmark adapter.

This module intentionally remains outside the public ``model`` CLI.  It maps one
admitted internal gap onto RFdiffusion's lossy chain/numbering model, restores
DVBFixer identities, and materializes canonical side-chain heavy atoms before
the generic diffusion pipeline performs independent validation.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import resource
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import NamedTuple

import numpy as np
from openmm import Platform
from openmm.app import PDBFile
from pdbfixer import PDBFixer

from dvbfixer.ffutils.geometry import rebuild_missing_atoms_with_retry
from dvbfixer.model.diffusion.boundary_refinement import (
    BOUNDARY_REFINEMENT_REVISION,
    refine_generated_region,
)
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    ResidueIdentity,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResourceMetrics,
    RunnerResult,
)
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.runner import (
    DIFFUSION_RUNNER_PROTOCOL_VERSION,
    REQUEST_MANIFEST,
    RESULT_MANIFEST,
)

RFDIFFUSION_REPOSITORY = "https://github.com/RosettaCommons/RFdiffusion"
RFDIFFUSION_REVISION = "bf42b54c20a99dd7350456c85985ed4d83b95d48"
RFDIFFUSION_CHECKPOINT_SHA256 = (
    "0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca"
)
RFDIFFUSION_ADAPTER_REVISION = "dvbfixer-rfdiffusion-v1-adapter-2"
RFDIFFUSION_ENVIRONMENT_SHA256 = (
    "e36300c79f489f420d069d9688953ee4830fe563b8edc77f5de458e77464fe23"
)

_BACKBONE_ATOMS = ("N", "CA", "C", "O")
PARTNER_CONTEXT_BACKBONE_RMSD_MAX_ANGSTROM = 0.5
PARTNER_CONTEXT_BACKBONE_DISPLACEMENT_MAX_ANGSTROM = 2.0
_ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


class RFdiffusionAdapterError(RuntimeError):
    """Raised when pinned RFdiffusion evidence or output is invalid."""


class _PDBAtom(NamedTuple):
    identity: AtomIdentity
    residue_name: str
    coordinate: np.ndarray
    element: str


@dataclass(frozen=True, slots=True)
class RFdiffusionV1Config:
    repository: Path
    python: Path
    checkpoint: Path
    checkpoint_sha256: str = RFDIFFUSION_CHECKPOINT_SHA256


@dataclass(frozen=True, slots=True)
class PartnerChainMapping:
    source_chain: str
    synthetic_chain: str
    residue_count: int
    output_start: int


@dataclass(frozen=True, slots=True)
class SyntheticInput:
    pdb_text: str
    contig: str
    target_chain: str
    target_length: int
    partner_chains: tuple[PartnerChainMapping, ...] = ()


def build_synthetic_input(request: DiffusionRequest, source_text: str) -> SyntheticInput:
    """Map the target and fixed protein partners into RFdiffusion identity space."""
    reason = unsupported_request_reason(request)
    if reason:
        raise RFdiffusionAdapterError(reason)

    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    gap = request.gaps[0]
    target_index_by_residue = dict(
        zip(placement.observed_residues, placement.observed_target_indices)
    )
    atom_lines = [line for line in source_text.splitlines() if line.startswith("ATOM  ")]
    seen_residues: set[ResidueIdentity] = set()
    output: list[str] = []
    serial = 1
    for line in atom_lines:
        identity = _residue_identity(line)
        target_index = target_index_by_residue.get(identity)
        if target_index is None:
            continue
        seen_residues.add(identity)
        remapped = _replace_pdb_identity(
            line,
            serial=serial,
            chain="A",
            residue_number=str(target_index + 1),
            insertion_code="",
        )
        output.append(remapped + "\n")
        serial += 1

    if seen_residues != set(placement.observed_residues):
        raise RFdiffusionAdapterError(
            "normalized PDB does not contain every residue in the sequence placement"
        )
    output.append(f"TER   {serial:5d}\n")
    serial += 1

    source_partner_chains = tuple(
        dict.fromkeys(line[21] for line in atom_lines if line[21] != target.chain)
    )
    if len(source_partner_chains) > 25:
        raise RFdiffusionAdapterError(
            "RFdiffusion v1 benchmark adapter supports at most 25 partner chains"
        )
    partner_chains: list[PartnerChainMapping] = []
    receptor_offset = 0
    for index, source_chain in enumerate(source_partner_chains, start=1):
        synthetic_chain = chr(ord("A") + index)
        partner_lines = [line for line in atom_lines if line[21] == source_chain]
        residue_ordinals: dict[ResidueIdentity, int] = {}
        for line in partner_lines:
            identity = _residue_identity(line)
            residue_number = residue_ordinals.setdefault(
                identity,
                len(residue_ordinals) + 1,
            )
            output.append(
                _replace_pdb_identity(
                    line,
                    serial=serial,
                    chain=synthetic_chain,
                    residue_number=str(residue_number),
                    insertion_code="",
                )
                + "\n"
            )
            serial += 1
        if not residue_ordinals:
            continue
        partner_chains.append(
            PartnerChainMapping(
                source_chain=source_chain,
                synthetic_chain=synthetic_chain,
                residue_count=len(residue_ordinals),
                output_start=len(target.sequence) + receptor_offset + 1,
            )
        )
        receptor_offset += len(residue_ordinals)
        output.append(f"TER   {serial:5d}\n")
        serial += 1

    gap_length = gap.target_interval.stop - gap.target_interval.start
    left_stop = gap.target_interval.start
    right_start = gap.target_interval.stop + 1
    contig = (
        f"A1-{left_stop}/{gap_length}-{gap_length}/"
        f"A{right_start}-{placement.target_length}"
    )
    if partner_chains:
        receptor_contigs = [
            f"{partner.synthetic_chain}1-{partner.residue_count}"
            for partner in partner_chains
        ]
        contig = "/0 ".join((contig, *receptor_contigs))
    output.append("END\n")
    return SyntheticInput(
        pdb_text="".join(output),
        contig=contig,
        target_chain=target.chain,
        target_length=len(target.sequence),
        partner_chains=tuple(partner_chains),
    )


def unsupported_request_reason(request: DiffusionRequest) -> str:
    """Return the first adapter-specific unsupported reason, if any."""
    if len(request.target_sequences) != 1 or len(request.sequence_placements) != 1:
        return "RFdiffusion v1 benchmark adapter supports exactly one target chain"
    if len(request.gaps) != 1:
        return "RFdiffusion v1 benchmark adapter supports exactly one internal gap"
    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    gap = request.gaps[0]
    if target.chain != placement.chain or target.chain != gap.chain:
        return "RFdiffusion target, placement, and gap chains must match"
    if gap.target_interval.start <= 0 or gap.target_interval.stop >= len(target.sequence):
        return "RFdiffusion v1 benchmark adapter requires two observed gap anchors"
    expected_observed = set(range(len(target.sequence))) - set(
        range(gap.target_interval.start, gap.target_interval.stop)
    )
    if set(placement.observed_target_indices) != expected_observed:
        return "RFdiffusion sequence placement must cover every non-generated target residue"
    if request.retained_explicit_links:
        return "RFdiffusion v1 benchmark adapter does not yet carry explicit links"
    try:
        generated_numbers = tuple(
            int(residue.residue_number) for residue in gap.generated_residues
        )
        if any(number < -999 or number > 9999 for number in generated_numbers):
            return "RFdiffusion generated residue numbers are outside the PDB range"
    except ValueError:
        return "RFdiffusion generated residue identities require integer PDB residue numbers"
    unknown_options = sorted(
        option.name
        for option in request.backend_options
        if option.name not in {"boundary_refinement", "final_step"}
    )
    if unknown_options:
        return "unsupported RFdiffusion backend options: " + ", ".join(unknown_options)
    return ""


def run_adapter(request: DiffusionRequest, config: RFdiffusionV1Config) -> RunnerResult:
    """Run one pinned RFdiffusion design per requested seed."""
    reason = unsupported_request_reason(request)
    provenance = _backend_provenance(config)
    if reason:
        return RunnerResult(
            schema_version=DIFFUSION_SCHEMA_VERSION,
            status=DiffusionStatus.UNSUPPORTED,
            candidates=(),
            runner_diagnostics=RunnerDiagnostics(exit_code=None, timed_out=False),
            backend_provenance=provenance,
            message=reason,
        )

    _verify_config(config)
    source_path = Path(request.normalized_pdb.path)
    source_text = source_path.read_text(encoding="utf-8")
    synthetic = build_synthetic_input(request, source_text)
    synthetic_path = Path("rfdiffusion") / "input.pdb"
    synthetic_path.parent.mkdir(parents=True, exist_ok=True)
    synthetic_path.write_text(synthetic.pdb_text, encoding="utf-8")

    started = time.monotonic()
    candidates: list[RunnerCandidate] = []
    stdout_parts: list[str] = []
    stderr_parts: list[str] = []
    peak_vram_bytes = 0
    final_step = _backend_option(request, "final_step", "1")
    if not final_step.isdigit() or int(final_step) < 1 or int(final_step) > 50:
        raise RFdiffusionAdapterError("final_step must be an integer between 1 and 50")
    boundary_refinement = _boolean_backend_option(
        request,
        "boundary_refinement",
        default=True,
    )

    for index, seed in enumerate(request.seeds, start=1):
        prefix = Path("rfdiffusion") / f"seed-{seed}" / "design"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        command = (
            str(config.python),
            str(config.repository / "run_inference.py"),
            f"inference.input_pdb={synthetic_path.resolve()}",
            f"inference.output_prefix={prefix.resolve()}",
            "inference.num_designs=1",
            f"inference.design_startnum={seed}",
            "inference.deterministic=True",
            "inference.write_trajectory=False",
            f"inference.final_step={final_step}",
            f"inference.ckpt_override_path={config.checkpoint.resolve()}",
            f"contigmap.contigs=[{synthetic.contig}]",
            f"hydra.run.dir={(prefix.parent / 'hydra').resolve()}",
            "hydra.output_subdir=null",
            "hydra.job.chdir=False",
        )
        returncode, stdout, stderr, invocation_peak_vram = _run_inference(
            command,
            cwd=config.repository,
        )
        peak_vram_bytes = max(peak_vram_bytes, invocation_peak_vram)
        stdout_parts.append(stdout)
        stderr_parts.append(stderr)
        if returncode != 0:
            raise RFdiffusionAdapterError(
                f"RFdiffusion failed for seed {seed} with exit code {returncode}: "
                f"{stderr[-2000:]}"
            )

        raw_pdb = prefix.with_name(f"{prefix.name}_{seed}.pdb").resolve()
        raw_trb = prefix.with_name(f"{prefix.name}_{seed}.trb").resolve()
        if not raw_pdb.is_file() or not raw_trb.is_file():
            raise RFdiffusionAdapterError(
                f"RFdiffusion did not produce the expected PDB/TRB pair for seed {seed}"
            )
        candidate_id = f"candidate-{index:04d}"
        candidate_path = Path("candidates") / f"{candidate_id}.pdb"
        candidate_path.parent.mkdir(parents=True, exist_ok=True)
        candidate_text = materialize_candidate(
            request,
            source_text=source_text,
            synthetic_text=synthetic.pdb_text,
            raw_text=raw_pdb.read_text(encoding="utf-8"),
            partner_chains=synthetic.partner_chains,
            boundary_refinement=boundary_refinement,
            refinement_seed=seed,
            refinement_platform="Reference",
        )
        candidate_path.write_text(candidate_text, encoding="utf-8")
        warnings = [
            "RFdiffusion v1 is backbone-only; DVBFixer materialized target-sequence "
            "side chains with PDBFixer before independent validation",
            "fixed source coordinates were reinserted after post-sampling Kabsch alignment",
        ]
        if synthetic.partner_chains:
            warnings.append(
                "fixed partner chains conditioned RFdiffusion as receptor context and "
                "passed post-sampling backbone-adherence checks"
            )
        warnings.append(
            "generated residues received seeded post-sampling OpenMM boundary "
            f"refinement ({BOUNDARY_REFINEMENT_REVISION}); this is not per-step "
            "diffusion reinjection"
            if boundary_refinement
            else "post-sampling OpenMM boundary refinement was disabled for this ablation"
        )
        candidates.append(
            RunnerCandidate(
                candidate_id=candidate_id,
                seed=seed,
                coordinate_artifact=ArtifactReference(
                    candidate_path.as_posix(), _sha256(candidate_path)
                ),
                generated_atoms=request.generated_atoms,
                generated_residues=request.gaps[0].generated_residues,
                raw_backend_score=None,
                score_provenance="rfdiffusion-v1:no-comparable-backend-score",
                warnings=tuple(warnings),
            )
        )

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_ram_bytes = int(usage.ru_maxrss * 1024) if sys.platform != "darwin" else int(usage.ru_maxrss)
    return RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=tuple(candidates),
        runner_diagnostics=RunnerDiagnostics(
            exit_code=0,
            timed_out=False,
            stdout="".join(stdout_parts),
            stderr="".join(stderr_parts),
        ),
        backend_provenance=provenance,
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=time.monotonic() - started,
            peak_ram_bytes=peak_ram_bytes,
            peak_vram_bytes=peak_vram_bytes or None,
        ),
    )


def materialize_candidate(
    request: DiffusionRequest,
    *,
    source_text: str,
    synthetic_text: str,
    raw_text: str,
    partner_chains: tuple[PartnerChainMapping, ...] = (),
    boundary_refinement: bool = False,
    refinement_seed: int = 1,
    refinement_platform: str = "Reference",
    refinement_restart_count: int | None = None,
) -> str:
    """Restore identities and complete generated canonical heavy atoms."""
    gap = request.gaps[0]
    target = request.target_sequences[0]
    sampled = _atoms_by_target_index(raw_text)
    authoritative = _atoms_by_target_index(synthetic_text)
    observed_indices = request.sequence_placements[0].observed_target_indices
    anchor_keys = tuple(
        (target_index, atom_name)
        for target_index in observed_indices
        for atom_name in _BACKBONE_ATOMS
        if atom_name in sampled.get(target_index, {})
        and atom_name in authoritative.get(target_index, {})
    )
    if len(anchor_keys) < 3:
        raise RFdiffusionAdapterError(
            "RFdiffusion output contains fewer than three mapped backbone anchors"
        )
    mobile_anchor_coordinates = [
        sampled[index][atom_name].coordinate for index, atom_name in anchor_keys
    ]
    target_anchor_coordinates = [
        authoritative[index][atom_name].coordinate for index, atom_name in anchor_keys
    ]
    partner_mobile, partner_target = _partner_backbone_anchors(
        synthetic_text,
        raw_text,
        partner_chains,
    )
    mobile_anchor_coordinates.extend(partner_mobile)
    target_anchor_coordinates.extend(partner_target)
    mobile_anchors = np.vstack(mobile_anchor_coordinates)
    target_anchors = np.vstack(target_anchor_coordinates)
    transform = weighted_kabsch(mobile_anchors, target_anchors)
    if partner_mobile:
        displacements = np.linalg.norm(
            transform.apply(np.vstack(partner_mobile)) - np.vstack(partner_target),
            axis=1,
        )
        rmsd = float(np.sqrt(np.mean(np.square(displacements))))
        maximum = float(np.max(displacements))
        if (
            rmsd > PARTNER_CONTEXT_BACKBONE_RMSD_MAX_ANGSTROM
            or maximum > PARTNER_CONTEXT_BACKBONE_DISPLACEMENT_MAX_ANGSTROM
        ):
            raise RFdiffusionAdapterError(
                "RFdiffusion partner context drifted after sampling: "
                f"backbone RMSD {rmsd:.3f} Angstrom, maximum displacement "
                f"{maximum:.3f} Angstrom"
            )

    generated_lines: list[str] = []
    serial = _maximum_serial(source_text) + 1
    for target_index, residue in zip(
        range(gap.target_interval.start, gap.target_interval.stop),
        gap.generated_residues,
    ):
        residue_atoms = sampled.get(target_index)
        if residue_atoms is None or set(_BACKBONE_ATOMS) - set(residue_atoms):
            raise RFdiffusionAdapterError(
                f"RFdiffusion output is missing backbone atoms at target index {target_index}"
            )
        residue_name = _ONE_TO_THREE[target.sequence[target_index]]
        coordinates = transform.apply(
            np.vstack([residue_atoms[name].coordinate for name in _BACKBONE_ATOMS])
        )
        for atom_name, coordinate in zip(_BACKBONE_ATOMS, coordinates):
            generated_lines.append(
                _format_atom(
                    serial,
                    atom_name,
                    residue_name,
                    residue,
                    coordinate,
                    atom_name[0],
                )
            )
            serial += 1

    backbone_candidate = _insert_generated_lines(source_text, gap.right_anchor, generated_lines)
    fixer = PDBFixer(pdbfile=StringIO(backbone_candidate))
    fixer.platform = Platform.getPlatformByName("CPU")
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    rebuild_missing_atoms_with_retry(fixer)
    rebuilt = StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, rebuilt, keepIds=True)
    rebuilt_atoms = {
        atom.identity: atom
        for atom in _parse_atoms(rebuilt.getvalue())
        if atom.identity in set(request.generated_atoms)
    }
    missing = set(request.generated_atoms) - set(rebuilt_atoms)
    if missing:
        labels = ", ".join(
            f"{atom.chain}/{atom.residue_number}{atom.insertion_code}/{atom.atom_name}"
            for atom in sorted(missing)
        )
        raise RFdiffusionAdapterError(
            f"side-chain materialization did not produce expected generated atoms: {labels}"
        )

    completed_lines: list[str] = []
    serial = _maximum_serial(source_text) + 1
    for identity in request.generated_atoms:
        atom = rebuilt_atoms[identity]
        completed_lines.append(
            _format_atom(
                serial,
                identity.atom_name,
                atom.residue_name,
                ResidueIdentity(identity.chain, identity.residue_number, identity.insertion_code),
                atom.coordinate,
                atom.element,
            )
        )
        serial += 1
    completed_candidate = _insert_generated_lines(
        source_text,
        gap.right_anchor,
        completed_lines,
    )
    if not boundary_refinement:
        return completed_candidate

    short_gap = len(gap.generated_residues) <= 5
    refinement_kwargs: dict[str, object] = {}
    if short_gap:
        refinement_kwargs["max_iterations"] = 500
        refinement_kwargs["perturbation_angstrom"] = 0.25
    if refinement_restart_count is not None:
        refinement_kwargs["restart_count"] = refinement_restart_count
    refinement = refine_generated_region(
        completed_candidate,
        generated_residues=gap.generated_residues,
        generated_atoms=request.generated_atoms,
        random_seed=refinement_seed,
        platform_name=refinement_platform,
        **refinement_kwargs,
    )
    refined_lines: list[str] = []
    serial = _maximum_serial(source_text) + 1
    for identity in request.generated_atoms:
        atom = rebuilt_atoms[identity]
        refined_lines.append(
            _format_atom(
                serial,
                identity.atom_name,
                atom.residue_name,
                ResidueIdentity(
                    identity.chain,
                    identity.residue_number,
                    identity.insertion_code,
                ),
                refinement.coordinates_angstrom[identity],
                atom.element,
            )
        )
        serial += 1
    return _insert_generated_lines(source_text, gap.right_anchor, refined_lines)


def _atoms_by_target_index(text: str) -> dict[int, dict[str, _PDBAtom]]:
    atoms: dict[int, dict[str, _PDBAtom]] = {}
    for atom in _parse_atoms(text):
        if atom.identity.chain != "A":
            continue
        try:
            target_index = int(atom.identity.residue_number) - 1
        except ValueError as exc:
            raise RFdiffusionAdapterError("RFdiffusion output residue number is not an integer") from exc
        residue_atoms = atoms.setdefault(target_index, {})
        if atom.identity.atom_name in residue_atoms:
            raise RFdiffusionAdapterError("RFdiffusion output contains duplicate atom identities")
        residue_atoms[atom.identity.atom_name] = atom
    return atoms


def _partner_backbone_anchors(
    synthetic_text: str,
    raw_text: str,
    partner_chains: tuple[PartnerChainMapping, ...],
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    synthetic_atoms = {atom.identity: atom for atom in _parse_atoms(synthetic_text)}
    raw_atoms = {atom.identity: atom for atom in _parse_atoms(raw_text)}
    mobile: list[np.ndarray] = []
    target: list[np.ndarray] = []
    for partner in partner_chains:
        for ordinal in range(1, partner.residue_count + 1):
            for atom_name in _BACKBONE_ATOMS:
                synthetic_identity = AtomIdentity(
                    partner.synthetic_chain,
                    str(ordinal),
                    "",
                    atom_name,
                )
                authoritative = synthetic_atoms.get(synthetic_identity)
                if authoritative is None:
                    continue
                raw_identity = AtomIdentity(
                    "B",
                    str(partner.output_start + ordinal - 1),
                    "",
                    atom_name,
                )
                sampled = raw_atoms.get(raw_identity)
                if sampled is None:
                    raise RFdiffusionAdapterError(
                        "RFdiffusion output is missing partner-context backbone atom "
                        f"{raw_identity.chain}/{raw_identity.residue_number}/{atom_name}"
                    )
                mobile.append(sampled.coordinate)
                target.append(authoritative.coordinate)
    return mobile, target


def _parse_atoms(text: str) -> tuple[_PDBAtom, ...]:
    atoms: list[_PDBAtom] = []
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            coordinate = np.asarray(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=np.float64,
            )
        except ValueError as exc:
            raise RFdiffusionAdapterError("PDB contains malformed coordinates") from exc
        atom_name = line[12:16].strip()
        atoms.append(
            _PDBAtom(
                AtomIdentity(
                    line[21],
                    line[22:26].strip(),
                    line[26].strip(),
                    atom_name,
                ),
                line[17:20].strip(),
                coordinate,
                line[76:78].strip() or atom_name[:1],
            )
        )
    return tuple(atoms)


def _insert_generated_lines(
    source_text: str,
    right_anchor: ResidueIdentity,
    generated_lines: list[str],
) -> str:
    output: list[str] = []
    inserted = False
    for line in source_text.splitlines(keepends=True):
        if (
            not inserted
            and line.startswith(("ATOM  ", "HETATM"))
            and _residue_identity(line) == right_anchor
        ):
            output.extend(generated_lines)
            inserted = True
        output.append(line)
    if not inserted:
        raise RFdiffusionAdapterError("right anchor is absent from normalized PDB")
    return "".join(output)


def _format_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    residue: ResidueIdentity,
    coordinate: np.ndarray,
    element: str,
) -> str:
    atom_field = f" {atom_name:<3}" if len(atom_name) < 4 and not atom_name[0].isdigit() else f"{atom_name:>4}"
    return (
        f"ATOM  {serial:5d} {atom_field} {residue_name:>3} {residue.chain}"
        f"{int(residue.residue_number):4d}{residue.insertion_code or ' '}   "
        f"{coordinate[0]:8.3f}{coordinate[1]:8.3f}{coordinate[2]:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )


def _replace_pdb_identity(
    line: str,
    *,
    serial: int,
    chain: str,
    residue_number: str,
    insertion_code: str,
) -> str:
    padded = line.ljust(80)
    return (
        padded[:6]
        + f"{serial:5d}"
        + padded[11:21]
        + chain
        + f"{int(residue_number):4d}"
        + (insertion_code or " ")
        + padded[27:80]
    )


def _residue_identity(line: str) -> ResidueIdentity:
    return ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())


def _maximum_serial(text: str) -> int:
    serials = [
        int(line[6:11])
        for line in text.splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and line[6:11].strip().isdigit()
    ]
    return max(serials, default=0)


def _backend_option(request: DiffusionRequest, name: str, default: str) -> str:
    return next((option.value for option in request.backend_options if option.name == name), default)


def _boolean_backend_option(
    request: DiffusionRequest,
    name: str,
    *,
    default: bool,
) -> bool:
    value = _backend_option(request, name, "true" if default else "false").lower()
    if value not in {"true", "false"}:
        raise RFdiffusionAdapterError(f"{name} must be true or false")
    return value == "true"


def _verify_config(config: RFdiffusionV1Config) -> None:
    if not config.python.is_file() or not os.access(config.python, os.X_OK):
        raise RFdiffusionAdapterError("RFdiffusion Python executable is unavailable")
    inference = config.repository / "run_inference.py"
    if not inference.is_file():
        raise RFdiffusionAdapterError("RFdiffusion run_inference.py is unavailable")
    revision = subprocess.run(
        ("git", "-C", str(config.repository), "rev-parse", "HEAD"),
        capture_output=True,
        text=True,
        check=False,
    )
    if revision.returncode != 0 or revision.stdout.strip() != RFDIFFUSION_REVISION:
        raise RFdiffusionAdapterError("RFdiffusion source revision does not match the pinned v1 commit")
    if not config.checkpoint.is_file():
        raise RFdiffusionAdapterError("RFdiffusion checkpoint is unavailable")
    actual = _sha256(config.checkpoint)
    if actual != config.checkpoint_sha256.lower():
        raise RFdiffusionAdapterError("RFdiffusion checkpoint SHA-256 mismatch")


def _run_inference(
    command: tuple[str, ...],
    *,
    cwd: Path,
) -> tuple[int, str, str, int]:
    process = subprocess.Popen(
        command,
        cwd=cwd,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stop = threading.Event()
    peak = [0]

    def sample_vram() -> None:
        while not stop.wait(0.2):
            peak[0] = max(peak[0], _reported_gpu_memory_bytes())

    monitor = threading.Thread(target=sample_vram, daemon=True)
    monitor.start()
    try:
        stdout, stderr = process.communicate()
    finally:
        stop.set()
        monitor.join()
    peak[0] = max(peak[0], _reported_gpu_memory_bytes())
    return process.returncode, stdout, stderr, peak[0]


def _reported_gpu_memory_bytes() -> int:
    try:
        completed = subprocess.run(
            (
                "nvidia-smi",
                "--query-compute-apps=used_gpu_memory",
                "--format=csv,noheader,nounits",
            ),
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return 0
    if completed.returncode != 0:
        return 0
    values = [
        int(line.strip())
        for line in completed.stdout.splitlines()
        if line.strip().isdigit()
    ]
    return sum(values) * 1024 * 1024


def _driver_version() -> str:
    try:
        completed = subprocess.run(
            ("nvidia-smi", "--query-gpu=driver_version", "--format=csv,noheader"),
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    return completed.stdout.splitlines()[0].strip() if completed.returncode == 0 else ""


def _backend_provenance(config: RFdiffusionV1Config) -> BackendProvenance:
    return BackendProvenance(
        backend="rfdiffusion-v1-backbone-benchmark",
        runner_protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
        engine_repository=RFDIFFUSION_REPOSITORY,
        engine_revision=RFDIFFUSION_REVISION,
        source_license="BSD-3-Clause",
        checkpoint_sha256=config.checkpoint_sha256,
        checkpoint_license="BSD-3-Clause; official upstream clarification postdates v1.0.0",
        environment_hash=RFDIFFUSION_ENVIRONMENT_SHA256,
        environment_identity=(
            f"{RFDIFFUSION_ADAPTER_REVISION}+{BOUNDARY_REFINEMENT_REVISION}"
        ),
        device="cuda:0",
        precision="float32",
        framework="pytorch",
        framework_version="1.9.0+cu111",
        cuda_version="11.1",
        driver_version=_driver_version(),
        deterministic_algorithms=False,
        deterministic_flags=(
            "inference.deterministic=True",
            "one-design-per-request-seed",
            "seeded hydrogen placement and OpenMM Reference boundary refinement",
        ),
        known_nondeterministic_operations=("legacy CUDA scatter/reduction kernels",),
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Internal pinned RFdiffusion v1 adapter")
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument(
        "--checkpoint-sha256",
        default=RFDIFFUSION_CHECKPOINT_SHA256,
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    request_name = os.environ.get("DVBFIXER_DIFFUSION_REQUEST", REQUEST_MANIFEST)
    result_name = os.environ.get("DVBFIXER_DIFFUSION_RESULT", RESULT_MANIFEST)
    request = DiffusionRequest.from_json(Path(request_name).read_text(encoding="utf-8"))
    config = RFdiffusionV1Config(
        repository=args.repository.resolve(),
        python=args.python.resolve(),
        checkpoint=args.checkpoint.resolve(),
        checkpoint_sha256=args.checkpoint_sha256,
    )
    try:
        result = run_adapter(request, config)
    except Exception as exc:
        result = RunnerResult(
            schema_version=DIFFUSION_SCHEMA_VERSION,
            status=DiffusionStatus.FAILED,
            candidates=(),
            runner_diagnostics=RunnerDiagnostics(
                exit_code=1,
                timed_out=False,
                stderr=f"{type(exc).__name__}: {exc}",
            ),
            backend_provenance=_backend_provenance(config),
            message=str(exc),
        )
    sys.stdout.write(result.runner_diagnostics.stdout)
    sys.stderr.write(result.runner_diagnostics.stderr)
    Path(result_name).write_text(result.to_json(), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
