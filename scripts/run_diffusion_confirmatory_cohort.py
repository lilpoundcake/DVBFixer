#!/usr/bin/env python3
"""Run and aggregate the locked diffusion confirmatory cohort.

This is research orchestration only. Scientific construction, inference,
refinement, MODELLER execution, validation, and metric calculation remain in
their existing entry points.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import shutil
import stat
import subprocess
import tempfile
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from statistics import median
from typing import Any

from dvbfixer.model.diffusion.benchmark import (
    BackendCaseResult,
    compare_paired_backends,
    summarize_backend_cases,
)

ROOT = Path(__file__).resolve().parents[1]
PROTENIX_REVISION = "85767b811c40ed46e73a9b39519cf6bfca8701ba"
PROTENIX_CHECKPOINT_SHA256 = (
    "2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04"
)
BOLTZ_REVISION = "b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc"
BOLTZ_CHECKPOINT_SHA256 = (
    "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1"
)
LOCKED_COHORT = ROOT / "docs/research/diffusion-confirmatory-cohort.json"
STAGE_SCHEMA_VERSION = 1
MATERIALIZED_FILES = (
    "request.json",
    "input/normalized.pdb",
    "reference.pdb",
    "target.fasta",
)


class CohortRunnerError(RuntimeError):
    """Raised when cohort configuration or a stage is invalid."""


@dataclass(frozen=True, slots=True)
class RunnerConfig:
    protenix_checkout: Path
    protenix_python: Path
    kalign: Path
    checkpoint: Path
    dvbfixer_python: Path
    modeller_python: Path
    boltz_checkout: Path
    boltz_python: Path
    boltz_cache: Path
    boltz_checkpoint: Path
    expected_protenix_revision: str = PROTENIX_REVISION
    expected_checkpoint_sha256: str = PROTENIX_CHECKPOINT_SHA256
    expected_boltz_revision: str = BOLTZ_REVISION
    expected_boltz_checkpoint_sha256: str = BOLTZ_CHECKPOINT_SHA256
    expected_kalign_sha256: str | None = None
    materialization_marker: str = "../.complete.json"
    cycles: int = 1
    steps: int = 200
    modeller_num_models: int = 1
    modeller_num_loops: int = 2
    modeller_num_output: int = 2


@dataclass(frozen=True, slots=True)
class StageSpec:
    name: str
    command: tuple[str, ...]
    inputs: tuple[Path, ...]
    outputs: tuple[Path, ...]
    cleanup: tuple[Path, ...]
    environment: tuple[tuple[str, str], ...] = ()
    peak_vram_reader: Callable[[], int | None] | None = None
    directories_to_create: tuple[Path, ...] = ()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


DirectoryManifest = tuple[tuple[str, str, int, int], ...]


def _directory_metadata_manifest(path: Path) -> DirectoryManifest:
    root_stat = path.lstat()
    if stat.S_ISLNK(root_stat.st_mode):
        raise CohortRunnerError(f"input directory must not be a symlink: {path}")
    if not stat.S_ISDIR(root_stat.st_mode):
        raise CohortRunnerError(f"input directory does not exist: {path}")
    entries = sorted(path.rglob("*"), key=lambda item: item.relative_to(path).as_posix())
    manifest: list[tuple[str, str, int, int]] = []
    for entry in entries:
        relative = entry.relative_to(path).as_posix()
        metadata = entry.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            raise CohortRunnerError(f"input directory contains a symlink: {entry}")
        if stat.S_ISDIR(metadata.st_mode):
            kind = "D"
        elif stat.S_ISREG(metadata.st_mode):
            kind = "F"
        else:
            raise CohortRunnerError(f"input directory contains a non-regular entry: {entry}")
        manifest.append((relative, kind, metadata.st_size, metadata.st_mtime_ns))
    return tuple(manifest)


@lru_cache(maxsize=16)
def _directory_content_sha256(root: str, manifest: DirectoryManifest) -> str:
    """Hash an immutable input tree after its metadata manifest changes.

    Cohort inputs, especially the Boltz cache, must remain immutable while a
    runner process is active. Every lookup still rebuilds the cheap metadata
    manifest; a path, kind, size, or mtime change invalidates this content cache.
    """
    path = Path(root)
    digest = hashlib.sha256()
    for relative, kind, _size, _mtime_ns in manifest:
        encoded = relative.encode("utf-8")
        digest.update(kind.encode("ascii") + b"\0" + encoded + b"\0")
        if kind == "F":
            digest.update(_sha256(path / relative).encode("ascii") + b"\0")
    return digest.hexdigest()


def _directory_tree_sha256(path: Path) -> str:
    absolute = path.absolute()
    manifest = _directory_metadata_manifest(absolute)
    return _directory_content_sha256(str(absolute), manifest)


def _atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    ) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def _load_config(path: Path) -> RunnerConfig:
    raw = json.loads(path.read_text(encoding="utf-8"))
    required = (
        "protenix_checkout",
        "protenix_python",
        "kalign",
        "checkpoint",
        "dvbfixer_python",
        "modeller_python",
        "boltz_checkout",
        "boltz_python",
        "boltz_cache",
        "boltz_checkpoint",
    )
    missing = [key for key in required if not raw.get(key)]
    if missing:
        raise CohortRunnerError(f"configuration is missing: {', '.join(missing)}")
    # Keep virtual-environment interpreter symlinks lexical: resolving them can
    # silently select the base interpreter and its unrelated site-packages.
    path_fields = {key: Path(raw.pop(key)).expanduser().absolute() for key in required}
    try:
        return RunnerConfig(**path_fields, **raw)
    except TypeError as exc:
        raise CohortRunnerError(f"invalid configuration: {exc}") from exc


def verify_config(config: RunnerConfig) -> dict[str, str]:
    for label, path in (
        ("Protenix checkout", config.protenix_checkout),
        ("Protenix Python", config.protenix_python),
        ("Kalign", config.kalign),
        ("checkpoint", config.checkpoint),
        ("DVBFixer Python", config.dvbfixer_python),
        ("MODELLER Python", config.modeller_python),
        ("Boltz checkout", config.boltz_checkout),
        ("Boltz Python", config.boltz_python),
        ("Boltz cache", config.boltz_cache),
        ("Boltz checkpoint", config.boltz_checkpoint),
    ):
        if not path.exists():
            raise CohortRunnerError(f"{label} does not exist: {path}")
    for label, path in (
        ("Protenix Python", config.protenix_python),
        ("Kalign", config.kalign),
        ("DVBFixer Python", config.dvbfixer_python),
        ("MODELLER Python", config.modeller_python),
        ("Boltz Python", config.boltz_python),
    ):
        if not path.is_file() or not os.access(path, os.X_OK):
            raise CohortRunnerError(f"{label} is not executable: {path}")
    completed = subprocess.run(
        ("git", "-C", str(config.protenix_checkout), "rev-parse", "HEAD"),
        check=False,
        capture_output=True,
        text=True,
        timeout=10,
    )
    revision = completed.stdout.strip().lower()
    if completed.returncode or revision != config.expected_protenix_revision:
        raise CohortRunnerError(
            f"Protenix revision mismatch: expected {config.expected_protenix_revision}, "
            f"found {revision or 'unavailable'}"
        )
    patch_check = subprocess.run(
        (
            "git",
            "-C",
            str(config.protenix_checkout),
            "apply",
            "--reverse",
            "--check",
            str(ROOT / "deploy/protenix-v1/per-step-callback.patch"),
        ),
        check=False,
        capture_output=True,
        text=True,
        timeout=10,
    )
    if patch_check.returncode:
        raise CohortRunnerError("the audited Protenix callback patch is not applied")
    boltz_completed = subprocess.run(
        ("git", "-C", str(config.boltz_checkout), "rev-parse", "HEAD"),
        check=False,
        capture_output=True,
        text=True,
        timeout=10,
    )
    boltz_revision = boltz_completed.stdout.strip().lower()
    if boltz_completed.returncode or boltz_revision != config.expected_boltz_revision:
        raise CohortRunnerError(
            f"Boltz revision mismatch: expected {config.expected_boltz_revision}, "
            f"found {boltz_revision or 'unavailable'}"
        )
    boltz_patch_check = subprocess.run(
        (
            "git",
            "-C",
            str(config.boltz_checkout),
            "apply",
            "--reverse",
            "--check",
            str(ROOT / "deploy/boltz-2/per-step-callback.patch"),
        ),
        check=False,
        capture_output=True,
        text=True,
        timeout=10,
    )
    if boltz_patch_check.returncode:
        raise CohortRunnerError("the audited Boltz callback patch is not applied")
    checkpoint_digest = _sha256(config.checkpoint)
    if checkpoint_digest != config.expected_checkpoint_sha256:
        raise CohortRunnerError("Protenix checkpoint digest mismatch")
    boltz_checkpoint_digest = _sha256(config.boltz_checkpoint)
    if boltz_checkpoint_digest != config.expected_boltz_checkpoint_sha256:
        raise CohortRunnerError("Boltz checkpoint digest mismatch")
    if not (config.boltz_cache / "ccd.pkl").is_file():
        raise CohortRunnerError("Boltz cache is missing ccd.pkl")
    if not (config.boltz_cache / "mols").is_dir() or not any(
        (config.boltz_cache / "mols").iterdir()
    ):
        raise CohortRunnerError("Boltz cache is missing extracted molecule data")
    kalign_digest = _sha256(config.kalign)
    if config.expected_kalign_sha256 and kalign_digest != config.expected_kalign_sha256:
        raise CohortRunnerError("Kalign digest mismatch")
    if config.cycles <= 0 or config.steps <= 0:
        raise CohortRunnerError("Protenix cycles and steps must be positive")
    if min(config.modeller_num_models, config.modeller_num_loops, config.modeller_num_output) <= 0:
        raise CohortRunnerError("MODELLER counts must be positive")
    if config.modeller_num_output > config.modeller_num_models * config.modeller_num_loops:
        raise CohortRunnerError("MODELLER output count exceeds its candidate pool")
    return {
        "protenix_revision": revision,
        "checkpoint_sha256": checkpoint_digest,
        "kalign_sha256": kalign_digest,
        "boltz_revision": boltz_revision,
        "boltz_checkpoint_sha256": boltz_checkpoint_digest,
    }


def _artifact_digest(value: object) -> str | None:
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        digest = value.get("sha256")
        return digest if isinstance(digest, str) else None
    return None


def load_materialization(workspace: Path, marker_name: str) -> tuple[Path, dict[str, Any]]:
    marker = workspace / marker_name
    if not marker.is_file():
        raise CohortRunnerError(f"materialization marker is missing: {marker}")
    completion = json.loads(marker.read_text(encoding="utf-8"))
    metadata = completion
    if isinstance(completion.get("record"), str):
        record_path = marker.parent / completion["record"]
        if (
            not record_path.is_file()
            or _sha256(record_path) != completion.get("record_sha256")
        ):
            raise CohortRunnerError("materialization completion record is missing or corrupt")
        metadata = json.loads(record_path.read_text(encoding="utf-8"))
    if metadata.get("status") not in {"accepted", "complete", "completed", "materialized"}:
        raise CohortRunnerError("materialization marker does not record an accepted completion")
    if metadata.get("selected_for_inference") is not True:
        raise CohortRunnerError("materialized case is not selected for inference")
    if not isinstance(metadata.get("case_id"), str) or not metadata["case_id"]:
        raise CohortRunnerError("materialization marker requires case_id")
    if not metadata.get("independence_group") or not isinstance(
        metadata.get("leakage_resolved"), bool
    ):
        cohort = json.loads(LOCKED_COHORT.read_text(encoding="utf-8"))
        matching = [
            case
            for case in cohort["cases"]
            if str(case["pdb_id"]).lower() == metadata["case_id"].lower()
        ]
        if len(matching) != 1:
            raise CohortRunnerError("case is not uniquely present in the locked cohort")
        metadata = {
            **metadata,
            "independence_group": matching[0]["independence_group"],
            "leakage_resolved": "protenix-v1" in cohort.get("eligible_backends", []),
        }
    artifacts = metadata.get("artifacts")
    if not isinstance(artifacts, dict):
        raise CohortRunnerError("materialization marker requires artifact digests")
    for relative in MATERIALIZED_FILES:
        path = workspace / relative
        expected = _artifact_digest(artifacts.get(relative))
        if not path.is_file() or expected is None or _sha256(path) != expected:
            raise CohortRunnerError(f"materialized artifact is missing or corrupt: {relative}")
    return marker, metadata


def _path_digests(paths: tuple[Path, ...], workspace: Path) -> dict[str, str]:
    digests: dict[str, str] = {}
    for path in paths:
        if not path.is_file() and not path.is_dir():
            raise CohortRunnerError(f"required artifact is missing: {path}")
        try:
            label = path.resolve().relative_to(workspace.resolve()).as_posix()
        except ValueError:
            label = str(path.resolve())
        if path.is_file():
            digests[label] = _sha256(path)
        else:
            digests[label] = _directory_tree_sha256(path)
    return digests


def _read_rss_tree(root_pid: int) -> int:
    parents: dict[int, int] = {}
    rss: dict[int, int] = {}
    proc = Path("/proc")
    if not proc.is_dir():
        return 0
    for entry in proc.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            fields = (entry / "stat").read_text().split()
            parents[int(entry.name)] = int(fields[3])
            for line in (entry / "status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    rss[int(entry.name)] = int(line.split()[1]) * 1024
                    break
        except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError):
            continue
    selected = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, parent in parents.items():
            if parent in selected and pid not in selected:
                selected.add(pid)
                changed = True
    return sum(rss.get(pid, 0) for pid in selected)


def _remove_owned(paths: tuple[Path, ...]) -> None:
    for path in paths:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        elif path.exists() or path.is_symlink():
            path.unlink()


def _marker_matches(marker: Path, spec: StageSpec, workspace: Path) -> bool:
    try:
        data = json.loads(marker.read_text(encoding="utf-8"))
        return (
            data["schema_version"] == STAGE_SCHEMA_VERSION
            and data["stage"] == spec.name
            and data["command"] == list(spec.command)
            and data.get("environment", {}) == dict(spec.environment)
            and data["exit_status"] == 0
            and data["inputs"] == _path_digests(spec.inputs, workspace)
            and data["outputs"] == _path_digests(spec.outputs, workspace)
        )
    except (CohortRunnerError, KeyError, OSError, ValueError, json.JSONDecodeError):
        return False


def run_stage(workspace: Path, spec: StageSpec) -> tuple[bool, dict[str, Any]]:
    state = workspace / ".confirmatory"
    marker = state / "stages" / f"{spec.name}.json"
    failure = state / "failures" / f"{spec.name}.json"
    if marker.is_file() and _marker_matches(marker, spec, workspace):
        return True, json.loads(marker.read_text(encoding="utf-8"))
    marker.unlink(missing_ok=True)
    _remove_owned(spec.cleanup)
    for directory in spec.directories_to_create:
        directory.mkdir(parents=True, exist_ok=True)
    inputs = _path_digests(spec.inputs, workspace)
    logs = state / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    stdout_path = logs / f"{spec.name}.stdout.log"
    stderr_path = logs / f"{spec.name}.stderr.log"
    environment = os.environ.copy()
    environment.update(spec.environment)
    start = time.perf_counter()
    peak_ram = 0
    exit_status = 127
    error: str | None = None
    try:
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            process = subprocess.Popen(
                spec.command,
                cwd=ROOT,
                env=environment,
                stdout=stdout,
                stderr=stderr,
            )
            while process.poll() is None:
                peak_ram = max(peak_ram, _read_rss_tree(process.pid))
                time.sleep(0.05)
            peak_ram = max(peak_ram, _read_rss_tree(process.pid))
            exit_status = process.returncode
    except OSError as exc:
        error = str(exc)
    wall_time = time.perf_counter() - start
    peak_vram = spec.peak_vram_reader() if exit_status == 0 and spec.peak_vram_reader else None
    resources = {
        "wall_time_seconds": wall_time,
        "peak_ram_bytes": peak_ram,
        "peak_vram_bytes": peak_vram,
    }
    base = {
        "schema_version": STAGE_SCHEMA_VERSION,
        "stage": spec.name,
        "command": list(spec.command),
        "environment": dict(spec.environment),
        "inputs": inputs,
        "exit_status": exit_status,
        "resources": resources,
        "stdout": stdout_path.relative_to(workspace).as_posix(),
        "stderr": stderr_path.relative_to(workspace).as_posix(),
    }
    try:
        outputs = _path_digests(spec.outputs, workspace) if exit_status == 0 else {}
    except CohortRunnerError as exc:
        outputs = {}
        error = str(exc)
    if exit_status != 0 or error is not None:
        record = {**base, "outputs": outputs, "error": error or "stage command failed"}
        _atomic_json(failure, record)
        return False, record
    record = {**base, "outputs": outputs}
    failure.unlink(missing_ok=True)
    _atomic_json(marker, record)
    return False, record


def _summary_vram(summary: Path) -> int | None:
    try:
        value = json.loads(summary.read_text(encoding="utf-8")).get("peak_vram_bytes")
    except (OSError, ValueError, json.JSONDecodeError):
        return None
    return value if isinstance(value, int) and value >= 0 else None


def _stage_specs(workspace: Path, config: RunnerConfig, marker: Path) -> list[StageSpec]:
    output = workspace / "confirmatory"
    protenix_input = output / "protenix-input"
    protenix_raw = output / "protenix-raw"
    protenix_refined = output / "protenix-refined"
    boltz_input = output / "boltz-input"
    boltz_raw = output / "boltz-raw"
    boltz_refined = output / "boltz-refined"
    modeller = output / "modeller"
    model_base = modeller / "model.pdb"
    model_paths = tuple(
        modeller / f"model_{index}.pdb" for index in range(1, config.modeller_num_output + 1)
    )
    model_sidecars = tuple(path.with_suffix(".dat") for path in model_paths)
    analysis = output / "modeller-analysis.json"
    protenix_python_path = os.pathsep.join(
        (str(config.protenix_checkout), str(ROOT / "src"), str(ROOT / "deploy/protenix-v1"))
    )
    boltz_python_path = os.pathsep.join(
        (str(ROOT / "src"), str(config.boltz_checkout / "src"), str(ROOT / "deploy/boltz-2"))
    )
    shared = (workspace / "request.json", workspace / "input/normalized.pdb", marker)
    return [
        StageSpec(
            "protenix-input",
            (
                str(config.dvbfixer_python),
                str(ROOT / "deploy/protenix-v1/build-template-input.py"),
                str(workspace / "request.json"),
                str(protenix_input),
            ),
            (*shared, ROOT / "deploy/protenix-v1/build-template-input.py"),
            (protenix_input / "input.json", protenix_input / "template.json"),
            (protenix_input,),
        ),
        StageSpec(
            "protenix-inference",
            (
                str(config.protenix_python),
                str(ROOT / "deploy/protenix-v1/checkpoint_gap_smoke.py"),
                str(workspace / "request.json"),
                str(protenix_input / "input.json"),
                str(protenix_raw),
                "--kalign",
                str(config.kalign),
                "--checkpoint",
                str(config.checkpoint),
                "--cycles",
                str(config.cycles),
                "--steps",
                str(config.steps),
            ),
            (
                *shared,
                protenix_input / "input.json",
                protenix_input / "template.json",
                config.checkpoint,
                config.kalign,
                ROOT / "deploy/protenix-v1/checkpoint_gap_smoke.py",
                ROOT / "deploy/protenix-v1/reinjection.py",
                ROOT / "deploy/protenix-v1/per-step-callback.patch",
            ),
            (protenix_raw / "candidate.pdb", protenix_raw / "summary.json"),
            (protenix_raw,),
            (("PYTHONPATH", protenix_python_path), ("PROTENIX_ROOT_DIR", str(config.protenix_checkout)), ("LAYERNORM_TYPE", "torch")),
            lambda: _summary_vram(protenix_raw / "summary.json"),
        ),
        StageSpec(
            "protenix-reference-refinement",
            (
                str(config.dvbfixer_python),
                str(ROOT / "deploy/protenix-v1/refine_candidate.py"),
                str(workspace / "request.json"),
                str(protenix_raw / "candidate.pdb"),
                str(protenix_refined),
                "--platform",
                "Reference",
                "--expected-candidate-sha256",
                "{raw_candidate_sha256}",
            ),
            (
                *shared,
                protenix_raw / "candidate.pdb",
                protenix_raw / "summary.json",
                ROOT / "deploy/protenix-v1/refine_candidate.py",
                ROOT / "src/dvbfixer/model/diffusion/boundary_refinement.py",
            ),
            (protenix_refined / "candidate.pdb", protenix_refined / "summary.json"),
            (protenix_refined,),
        ),
        StageSpec(
            "boltz-input",
            (
                str(config.dvbfixer_python),
                str(ROOT / "deploy/boltz-2/build-template-input.py"),
                str(workspace / "request.json"),
                str(boltz_input),
            ),
            (
                *shared,
                ROOT / "deploy/boltz-2/build-template-input.py",
                ROOT / "deploy/boltz-2/atom_mapping.py",
            ),
            (boltz_input / "input.yaml", boltz_input / "template.pdb"),
            (boltz_input,),
            (("PYTHONPATH", boltz_python_path),),
        ),
        StageSpec(
            "boltz-inference",
            (
                str(config.boltz_python),
                str(ROOT / "deploy/boltz-2/checkpoint_gap_smoke.py"),
                str(workspace / "request.json"),
                str(boltz_input / "input.yaml"),
                str(boltz_raw),
                "--cache",
                str(config.boltz_cache),
                "--checkpoint",
                str(config.boltz_checkpoint),
                "--recycling-steps",
                "1",
                "--sampling-steps",
                "200",
            ),
            (
                *shared,
                boltz_input / "input.yaml",
                boltz_input / "template.pdb",
                config.boltz_cache,
                config.boltz_checkpoint,
                ROOT / "deploy/boltz-2/checkpoint_gap_smoke.py",
                ROOT / "deploy/boltz-2/atom_mapping.py",
                ROOT / "deploy/boltz-2/per-step-callback.patch",
            ),
            (boltz_raw / "candidate.pdb", boltz_raw / "summary.json"),
            (boltz_raw,),
            (("PYTHONPATH", boltz_python_path),),
            lambda: _summary_vram(boltz_raw / "summary.json"),
        ),
        StageSpec(
            "boltz-reference-refinement",
            (
                str(config.dvbfixer_python),
                str(ROOT / "deploy/boltz-2/refine_candidate.py"),
                str(workspace / "request.json"),
                str(boltz_raw / "candidate.pdb"),
                str(boltz_refined),
                "--platform",
                "Reference",
                "--expected-candidate-sha256",
                "{raw_candidate_sha256}",
            ),
            (
                *shared,
                boltz_raw / "candidate.pdb",
                boltz_raw / "summary.json",
                ROOT / "deploy/boltz-2/refine_candidate.py",
                ROOT / "src/dvbfixer/model/diffusion/boundary_refinement.py",
            ),
            (boltz_refined / "candidate.pdb", boltz_refined / "summary.json"),
            (boltz_refined,),
        ),
        StageSpec(
            "modeller-comparator",
            (
                str(config.modeller_python),
                str(ROOT / "scripts/run_locked_modeller_comparator.py"),
                str(workspace / "request.json"),
                str(workspace / "input/normalized.pdb"),
                "--fasta",
                str(workspace / "target.fasta"),
                "--num-models",
                str(config.modeller_num_models),
                "--num-loops",
                str(config.modeller_num_loops),
                "--num-output",
                str(config.modeller_num_output),
                "--md-level",
                "fast",
                "--output",
                str(model_base),
            ),
            (
                *shared,
                workspace / "target.fasta",
                ROOT / "scripts/run_locked_modeller_comparator.py",
                ROOT / "src/dvbfixer/model/cli.py",
                ROOT / "src/dvbfixer/model/diffusion/contract.py",
                ROOT / "src/dvbfixer/model/pipeline.py",
                ROOT / "src/dvbfixer/model/modeller_run.py",
            ),
            (*model_paths, *model_sidecars),
            (modeller,),
            directories_to_create=(modeller,),
        ),
        StageSpec(
            "analysis",
            (
                str(config.dvbfixer_python),
                str(ROOT / "scripts/analyze_modeller_benchmark.py"),
                str(workspace),
                *(str(path) for path in model_paths),
                "--output",
                str(analysis),
            ),
            (
                *shared,
                workspace / "reference.pdb",
                protenix_refined / "summary.json",
                *model_paths,
                ROOT / "scripts/analyze_modeller_benchmark.py",
                ROOT / "src/dvbfixer/model/diffusion/benchmark.py",
                ROOT / "src/dvbfixer/model/diffusion/validate.py",
            ),
            (analysis,),
            (analysis,),
        ),
    ]


def _resolve_refinement_digest(spec: StageSpec) -> StageSpec:
    digest = _sha256(Path(spec.command[3]))
    return StageSpec(
        spec.name,
        tuple(digest if item == "{raw_candidate_sha256}" else item for item in spec.command),
        spec.inputs,
        spec.outputs,
        spec.cleanup,
        spec.environment,
        spec.peak_vram_reader,
        spec.directories_to_create,
    )


def run_case(workspace: Path, config: RunnerConfig) -> bool:
    workspace = workspace.resolve()
    try:
        marker, metadata = load_materialization(workspace, config.materialization_marker)
    except Exception as exc:
        _atomic_json(
            workspace / ".confirmatory/failures/materialization.json",
            {"stage": "materialization", "error": str(exc), "exit_status": None},
        )
        return False
    specs = _stage_specs(workspace, config, marker)
    outcomes: dict[str, bool] = {}
    for spec in specs:
        if spec.name == "protenix-inference" and not outcomes.get("protenix-input", False):
            continue
        if spec.name == "protenix-reference-refinement":
            if not outcomes.get("protenix-inference", False):
                continue
            spec = _resolve_refinement_digest(spec)
        if spec.name == "boltz-inference" and not outcomes.get("boltz-input", False):
            continue
        if spec.name == "boltz-reference-refinement":
            if not outcomes.get("boltz-inference", False):
                continue
            spec = _resolve_refinement_digest(spec)
        if spec.name == "analysis" and not (
            outcomes.get("protenix-reference-refinement", False)
            and outcomes.get("modeller-comparator", False)
        ):
            continue
        if spec.name == "analysis" and outcomes.get("boltz-reference-refinement", False):
            spec = StageSpec(
                spec.name,
                spec.command,
                (*spec.inputs, workspace / "confirmatory/boltz-refined/summary.json"),
                spec.outputs,
                spec.cleanup,
                spec.environment,
                spec.peak_vram_reader,
                spec.directories_to_create,
            )
        try:
            resumed, record = run_stage(workspace, spec)
        except Exception as exc:
            _atomic_json(
                workspace / ".confirmatory/failures" / f"{spec.name}.json",
                {
                    "schema_version": STAGE_SCHEMA_VERSION,
                    "stage": spec.name,
                    "command": list(spec.command),
                    "error": str(exc),
                    "exit_status": None,
                    "outputs": {},
                },
            )
            outcomes[spec.name] = False
            continue
        outcomes[spec.name] = resumed or (
            record["exit_status"] == 0 and "error" not in record
        )
    confirmatory_complete = all(
        outcomes.get(stage, False)
        for stage in (
            "protenix-input",
            "protenix-inference",
            "protenix-reference-refinement",
            "modeller-comparator",
            "analysis",
        )
    )
    proxy_complete = all(
        outcomes.get(stage, False)
        for stage in ("boltz-input", "boltz-inference", "boltz-reference-refinement")
    )
    summary = {
        "case_id": metadata["case_id"],
        "independence_group": metadata["independence_group"],
        "stages": outcomes,
        "confirmatory_complete": confirmatory_complete,
        "proxy_complete": proxy_complete,
        "complete": confirmatory_complete and proxy_complete,
    }
    _atomic_json(workspace / ".confirmatory/case.json", summary)
    return bool(summary["complete"])


def _resource_record(workspace: Path, stage: str) -> dict[str, Any] | None:
    state = workspace / ".confirmatory"
    for directory in ("stages", "failures"):
        path = state / directory / f"{stage}.json"
        if path.is_file():
            return json.loads(path.read_text(encoding="utf-8"))
    return None


def _backend_resources(workspace: Path, stages: tuple[str, ...]) -> tuple[float | None, int | None, int | None]:
    records = [_resource_record(workspace, stage) for stage in stages]
    observed = [record["resources"] for record in records if record and "resources" in record]
    if not observed:
        return None, None, None
    wall = sum(float(item["wall_time_seconds"]) for item in observed)
    ram = max(int(item["peak_ram_bytes"]) for item in observed)
    vrams = [int(item["peak_vram_bytes"]) for item in observed if item["peak_vram_bytes"] is not None]
    return wall, ram, max(vrams) if vrams else None


def _case_results(
    workspace: Path,
    marker_name: str,
) -> tuple[BackendCaseResult, BackendCaseResult, BackendCaseResult, dict[str, Any]]:
    _, metadata = load_materialization(workspace, marker_name)
    case_id = metadata["case_id"]
    group = metadata["independence_group"]
    leakage = metadata["leakage_resolved"]
    refined_path = workspace / "confirmatory/protenix-refined/summary.json"
    boltz_refined_path = workspace / "confirmatory/boltz-refined/summary.json"
    modeller_path = workspace / "confirmatory/modeller-analysis.json"
    protenix_wall, protenix_ram, protenix_vram = _backend_resources(
        workspace,
        ("protenix-input", "protenix-inference", "protenix-reference-refinement"),
    )
    boltz_wall, boltz_ram, boltz_vram = _backend_resources(
        workspace,
        ("boltz-input", "boltz-inference", "boltz-reference-refinement"),
    )
    modeller_wall, modeller_ram, _ = _backend_resources(workspace, ("modeller-comparator",))
    if refined_path.is_file():
        refined = json.loads(refined_path.read_text(encoding="utf-8"))
        validation = refined["validation"]
        quality = refined["quality"]
        protenix = BackendCaseResult(
            "protenix-v1",
            case_id,
            group,
            bool(validation["passed"]),
            float(quality["fixed_heavy_rmsd_angstrom"]) == 0.0,
            True,
            leakage,
            float(quality["gap_backbone_rmsd_angstrom"]),
            protenix_wall,
            protenix_vram,
            protenix_ram,
            True,
        )
    else:
        protenix = BackendCaseResult(
            "protenix-v1", case_id, group, False, False, True, leakage,
            wall_time_seconds=protenix_wall, peak_vram_bytes=protenix_vram,
            peak_ram_bytes=protenix_ram,
        )
    if boltz_refined_path.is_file():
        boltz_refined = json.loads(boltz_refined_path.read_text(encoding="utf-8"))
        boltz_validation = boltz_refined["validation"]
        boltz_quality = boltz_refined["quality"]
        boltz = BackendCaseResult(
            "boltz-2-proxy-only",
            case_id,
            group,
            bool(boltz_validation["passed"]),
            float(boltz_quality["fixed_heavy_rmsd_angstrom"]) == 0.0,
            True,
            False,
            float(boltz_quality["gap_backbone_rmsd_angstrom"]),
            boltz_wall,
            boltz_vram,
            boltz_ram,
            True,
        )
    else:
        boltz = BackendCaseResult(
            "boltz-2-proxy-only",
            case_id,
            group,
            False,
            False,
            True,
            False,
            wall_time_seconds=boltz_wall,
            peak_vram_bytes=boltz_vram,
            peak_ram_bytes=boltz_ram,
        )
    if modeller_path.is_file():
        modeller_summary = json.loads(modeller_path.read_text(encoding="utf-8"))
        top = modeller_summary["candidates"][0]
        modeller_result = BackendCaseResult(
            "modeller-10.8",
            case_id,
            group,
            bool(top["validation_passed"]),
            float(top["fixed_heavy_rmsd_angstrom"]) == 0.0,
            True,
            leakage,
            float(top["gap_backbone_rmsd_angstrom"]),
            modeller_wall,
            None,
            modeller_ram,
            False,
        )
    else:
        modeller_result = BackendCaseResult(
            "modeller-10.8", case_id, group, False, False, True, leakage,
            wall_time_seconds=modeller_wall, peak_ram_bytes=modeller_ram,
            gpu_memory_applicable=False,
        )
    resources = {
        stage: record["resources"]
        for stage in (
            "protenix-input",
            "protenix-inference",
            "protenix-reference-refinement",
            "boltz-input",
            "boltz-inference",
            "boltz-reference-refinement",
            "modeller-comparator",
            "analysis",
        )
        if (record := _resource_record(workspace, stage)) is not None and "resources" in record
    }
    return protenix, modeller_result, boltz, resources


def _proxy_descriptive_comparison(
    proxy_runs: tuple[BackendCaseResult, ...],
    comparison_runs: tuple[BackendCaseResult, ...],
) -> dict[str, Any]:
    proxy_by_case = {run.case_id: run for run in proxy_runs}
    comparison_by_case = {run.case_id: run for run in comparison_runs}
    if set(proxy_by_case) != set(comparison_by_case):
        raise CohortRunnerError("proxy comparison requires identical case IDs")
    pairs = [
        (proxy_by_case[case_id], comparison_by_case[case_id])
        for case_id in sorted(proxy_by_case)
    ]
    both_passed = sum(first.validation_passed and second.validation_passed for first, second in pairs)
    proxy_only = sum(first.validation_passed and not second.validation_passed for first, second in pairs)
    comparison_only = sum(
        not first.validation_passed and second.validation_passed for first, second in pairs
    )
    paired_differences = [
        first.gap_backbone_rmsd_angstrom - second.gap_backbone_rmsd_angstrom
        for first, second in pairs
        if first.validation_passed
        and second.validation_passed
        and first.gap_backbone_rmsd_angstrom is not None
        and second.gap_backbone_rmsd_angstrom is not None
    ]
    proxy_rmsds = [
        first.gap_backbone_rmsd_angstrom
        for first, _ in pairs
        if first.validation_passed and first.gap_backbone_rmsd_angstrom is not None
    ]
    comparison_rmsds = [
        second.gap_backbone_rmsd_angstrom
        for _, second in pairs
        if second.validation_passed and second.gap_backbone_rmsd_angstrom is not None
    ]
    return {
        "comparison_backend_id": comparison_runs[0].backend_id,
        "paired_case_count": len(pairs),
        "both_passed_count": both_passed,
        "proxy_only_passed_count": proxy_only,
        "comparison_only_passed_count": comparison_only,
        "neither_passed_count": len(pairs) - both_passed - proxy_only - comparison_only,
        "proxy_failure_count": sum(not first.validation_passed for first, _ in pairs),
        "comparison_failure_count": sum(not second.validation_passed for _, second in pairs),
        "both_valid_rmsd_count": len(paired_differences),
        "proxy_median_valid_gap_backbone_rmsd_angstrom": (
            float(median(proxy_rmsds)) if proxy_rmsds else None
        ),
        "comparison_median_valid_gap_backbone_rmsd_angstrom": (
            float(median(comparison_rmsds)) if comparison_rmsds else None
        ),
        "median_paired_rmsd_difference_angstrom": (
            float(median(paired_differences)) if paired_differences else None
        ),
    }


def _aggregation_inputs(cohort_root: Path, marker_name: str, output: Path) -> tuple[Path, ...]:
    inputs: list[Path] = [ROOT / "src/dvbfixer/model/diffusion/benchmark.py"]
    for workspace in _workspaces(cohort_root, marker_name):
        inputs.append((workspace / marker_name).resolve())
        for relative in (
            "confirmatory/protenix-refined/summary.json",
            "confirmatory/boltz-refined/summary.json",
            "confirmatory/modeller-analysis.json",
        ):
            path = workspace / relative
            if path.is_file():
                inputs.append(path)
        state = workspace / ".confirmatory"
        for directory in ("stages", "failures"):
            if (state / directory).is_dir():
                inputs.extend(sorted((state / directory).glob("*.json")))
    return tuple(path for path in inputs if path.resolve() != output.resolve())


def aggregate(cohort_root: Path, marker_name: str, output: Path, *, bootstrap_samples: int) -> dict[str, Any]:
    completion = output.with_suffix(output.suffix + ".complete.json")
    command = ["internal:aggregate", str(cohort_root), str(output), str(bootstrap_samples)]
    inputs = _path_digests(_aggregation_inputs(cohort_root, marker_name, output), cohort_root)
    if completion.is_file() and output.is_file():
        try:
            previous = json.loads(completion.read_text(encoding="utf-8"))
            if (
                previous["command"] == command
                and previous["inputs"] == inputs
                and previous["outputs"] == _path_digests((output,), cohort_root)
                and previous["exit_status"] == 0
            ):
                return json.loads(output.read_text(encoding="utf-8"))
        except (KeyError, OSError, ValueError, json.JSONDecodeError):
            pass
    completion.unlink(missing_ok=True)
    start = time.perf_counter()
    candidate: list[BackendCaseResult] = []
    baseline: list[BackendCaseResult] = []
    proxy: list[BackendCaseResult] = []
    resources: dict[str, Any] = {}
    failures: list[dict[str, Any]] = []
    for workspace in _workspaces(cohort_root, marker_name):
        try:
            protenix, modeller, boltz, case_resources = _case_results(workspace, marker_name)
        except Exception as exc:
            failures.append({"workspace": str(workspace), "error": str(exc)})
            continue
        candidate.append(protenix)
        baseline.append(modeller)
        proxy.append(boltz)
        resources[protenix.case_id] = case_resources
    if failures:
        raise CohortRunnerError(f"cannot aggregate {len(failures)} corrupt case workspace(s)")
    comparison = compare_paired_backends(
        tuple(candidate),
        tuple(baseline),
        minimum_independence_groups=180,
        maximum_independence_groups=300,
        minimum_both_valid_independence_groups=30,
        bootstrap_samples=bootstrap_samples,
    )
    proxy_summary = summarize_backend_cases(tuple(proxy))
    result = {
        "schema_version": 1,
        "preregistered_thresholds": {
            "minimum_independence_groups": 180,
            "maximum_independence_groups": 300,
            "minimum_both_valid_independence_groups": 30,
        },
        "comparison": asdict(comparison),
        "failure_rates": {
            comparison.candidate.backend_id: 1.0 - comparison.candidate.validation_pass_rate,
            comparison.baseline.backend_id: 1.0 - comparison.baseline.validation_pass_rate,
            proxy_summary.backend_id: 1.0 - proxy_summary.validation_pass_rate,
        },
        "case_results": {
            "protenix-v1": [asdict(item) for item in candidate],
            "modeller-10.8": [asdict(item) for item in baseline],
            "boltz-2-proxy-only": [asdict(item) for item in proxy],
        },
        "proxy_only": {
            "backend_id": "boltz-2-proxy-only",
            "purpose": "sensitivity evidence only; excluded from confirmatory selection",
            "exclusion_reason": "unresolved-training-leakage",
            "backend_summary": asdict(proxy_summary),
            "descriptive_comparisons": {
                "protenix-v1": _proxy_descriptive_comparison(
                    tuple(proxy), tuple(candidate)
                ),
                "modeller-10.8": _proxy_descriptive_comparison(
                    tuple(proxy), tuple(baseline)
                ),
            },
        },
        "resources_by_case_and_stage": resources,
    }
    _atomic_json(output, result)
    peak_ram = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    _atomic_json(
        completion,
        {
            "schema_version": STAGE_SCHEMA_VERSION,
            "stage": "aggregation",
            "command": command,
            "inputs": inputs,
            "outputs": _path_digests((output,), cohort_root),
            "exit_status": 0,
            "resources": {
                "wall_time_seconds": time.perf_counter() - start,
                "peak_ram_bytes": peak_ram,
                "peak_vram_bytes": None,
            },
        },
    )
    return json.loads(output.read_text(encoding="utf-8"))


def _workspaces(root: Path, marker_name: str) -> list[Path]:
    cases = root / "cases"
    if cases.is_dir():
        candidates = [
            case / "workspace"
            for case in cases.iterdir()
            if (case / "workspace").is_dir() and (case / ".complete.json").is_file()
        ]
    else:
        candidates = [
            path for path in root.iterdir()
            if path.is_dir() and (path / marker_name).is_file()
        ]
    selected: list[tuple[int, Path]] = []
    for workspace in candidates:
        try:
            _marker, metadata = load_materialization(workspace, marker_name)
            screening_index = metadata["screening_index"]
            if not isinstance(screening_index, int) or screening_index < 0:
                continue
            selected.append((screening_index, workspace))
        except (CohortRunnerError, KeyError, OSError, ValueError, json.JSONDecodeError):
            continue
    return [workspace for _index, workspace in sorted(selected)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("cohort_root", type=Path)
    run_parser.add_argument("--config", required=True, type=Path)
    run_parser.add_argument("--case", action="append", default=[])
    aggregate_parser = subparsers.add_parser("aggregate")
    aggregate_parser.add_argument("cohort_root", type=Path)
    aggregate_parser.add_argument("--output", required=True, type=Path)
    aggregate_parser.add_argument("--materialization-marker", default="../.complete.json")
    aggregate_parser.add_argument("--bootstrap-samples", default=10_000, type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.action == "aggregate":
        aggregate(
            args.cohort_root.resolve(),
            args.materialization_marker,
            args.output.resolve(),
            bootstrap_samples=args.bootstrap_samples,
        )
        return 0
    config = _load_config(args.config)
    verified = verify_config(config)
    print(json.dumps({"verified": verified}, sort_keys=True))
    root = args.cohort_root.resolve()
    workspaces = _workspaces(root, config.materialization_marker)
    if args.case:
        selected = set(args.case)
        workspaces = [workspace for workspace in workspaces if workspace.parent.name in selected]
    failed = 0
    for workspace in workspaces:  # Deliberately sequential: one GPU case at a time.
        print(f"CASE {workspace.name}", flush=True)
        if not run_case(workspace, config):
            failed += 1
    print(json.dumps({"cases": len(workspaces), "failed": failed}, sort_keys=True))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
