"""Provision isolated diffusion/QM runtimes through the shared manager adapters.

Weights and Protenix data/Kalign are operator-provided. Registration is staged,
checked by the existing runner handshake, and only then published.
"""

from __future__ import annotations

import ast
import hashlib
import json
import os
import shutil
import stat
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Any

from dvbfixer.install import Manager, child_environment, platform_subdir, reserve_prefix
from dvbfixer.installation_profiles import load_catalog

SNAPSHOT_CODE = (
    "import importlib.metadata as m,json,platform,sys; from pathlib import Path; "
    "print(json.dumps({'python':platform.python_version(),'prefix':sys.prefix,"
    "'packages':{d.metadata['Name']:d.version for d in m.distributions()},"
    "'conda':[ {k:r.get(k) for k in ('name','version','build','channel')} "
    "for p in Path(sys.prefix,'conda-meta').glob('*.json') "
    "for r in [json.loads(p.read_text())]]},sort_keys=True))"
)


def runtime_names() -> tuple[str, ...]:
    return tuple(name for name, entry in load_catalog()["external"].items()
                 if entry["kind"] in {"separate-runtime", "separate-environment"})


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def checked_file(path: Path, digest: str, label: str) -> Path:
    path = path.expanduser().absolute()
    if path.is_symlink() or not path.is_file() or sha256(path) != digest:
        raise ValueError(f"{label} must be a regular file matching the frozen SHA256")
    return path


def asset_requirements(entry: dict[str, Any], repository: Path) -> dict[str, Any]:
    if "profile_lock" in entry:
        lock = json.loads((repository / entry["profile_lock"]).read_text())
        return {"checkpoint": lock["checkpoint"], "kalign": lock["kalign"],
                "auxiliary": lock["auxiliary_artifact_sha256"],
                "sources": lock["patched_source_sha256"], "patch": lock["patch_sha256"]}
    def literal(relative: str, name: str) -> Any:
        # Reuse the runner's frozen literals without importing its scientific stack.
        tree = ast.parse((repository / relative).read_text())
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
                return ast.literal_eval(node.value)
        raise ValueError(f"Runner identity is unavailable: {relative}:{name}")

    sources = literal("deploy/protpardelle-1c/cuda_runner.py", "SOURCE_HASHES") if entry["device"] == "cuda" else {
        "src/protpardelle/core/models.py": literal("deploy/protpardelle-1c/production_runner.py", "PATCHED_MODELS_SHA256")}
    identity = "src/dvbfixer/model/diffusion/cc94.py"
    result = {"checkpoint": {"name": literal(identity, "CHECKPOINT_NAME"),
                             "sha256": literal(identity, "CHECKPOINT_SHA256")},
              "config": literal(identity, "CONFIG_SHA256"), "sources": sources}
    if entry["device"] == "mps":
        result["patch"] = literal("deploy/protpardelle-1c/production_runner.py", "APPLE_PATCH_SHA256")
    return result


def verify_assets(requirements: dict[str, Any], *, checkpoint: Path | None,
                  kalign: Path | None, auxiliary: Path | None) -> list[str]:
    """Validate supplied material before expensive installation; missing material is deferred."""
    missing = []
    if checkpoint is None:
        missing.append("Provide the frozen checkpoint with --checkpoint; weights are never downloaded.")
    else:
        checked_file(checkpoint, requirements["checkpoint"]["sha256"], "Checkpoint")
        if checkpoint.name != requirements["checkpoint"]["name"]:
            raise ValueError(f"Checkpoint must retain its frozen filename: {requirements['checkpoint']['name']}")
        if "config" in requirements:
            checked_file(checkpoint.parent.parent / "configs/cc94.yaml", requirements["config"], "cc94 config")
    if "kalign" in requirements:
        if kalign is None:
            missing.append("Provide the digest-pinned executable with --kalign.")
        else:
            checked_file(kalign, requirements["kalign"]["sha256"], "Kalign")
            if not os.access(kalign, os.X_OK):
                raise ValueError("Kalign must be executable")
        if auxiliary is None:
            missing.append("Provide the frozen common/ data with --auxiliary-root.")
        else:
            for name, digest in requirements["auxiliary"].items():
                checked_file(auxiliary / name, digest, f"Protenix data {name}")
    return missing


def runtime_environment(entry: dict[str, Any], repository: Path, engine: Path) -> dict[str, str]:
    env = child_environment()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["GIT_LFS_SKIP_SMUDGE"] = "1"
    env["PYTHONPATH"] = os.pathsep.join((str(engine / entry["source_subdir"]), str(repository / "src")))
    env["LAYERNORM_TYPE"] = "torch"
    env["PROTENIX_ROOT_DIR"] = str(engine)
    if entry["device"] == "mps":
        env["PYTORCH_ENABLE_MPS_FALLBACK"] = "0"
        for name in ("PYTORCH_MPS_FAST_MATH", "PYTORCH_MPS_PREFER_METAL",
                     "PYTORCH_MPS_HIGH_WATERMARK_RATIO", "PYTORCH_MPS_LOW_WATERMARK_RATIO"):
            env.pop(name, None)
    return env


def owned_child(root: Path, name: str) -> Path:
    path = root / name
    if path.is_symlink() or not path.resolve().is_relative_to(root):
        raise ValueError(f"Runtime-owned path must not escape or be a symlink: {name}")
    return path


def prepare_source(entry: dict[str, Any], repository: Path, engine: Path,
                   source: Path | None, log: Path, env: dict[str, str],
                   requirements: dict[str, Any]) -> None:
    """Clone into our new root; never checkout/reset/patch an operator's source tree."""
    def run(argv: list[str], *, capture: bool = False) -> str:
        with log.open("a") as stream:
            result = subprocess.run(argv, env=env, check=True, text=True, timeout=300,
                                    stdout=subprocess.PIPE if capture else stream, stderr=stream)
        return result.stdout if capture else ""

    if not engine.exists():
        run(["git", "clone", "--no-hardlinks", "--",
             str(source) if source else entry["repository"], str(engine)])
    # No reset/clean is used on resume: an unexpected source state is an error.
    if not (engine / ".git").is_dir():
        raise ValueError("Incomplete engine clone: retain it for inspection and choose a new runtime prefix")
    observed = run(["git", "-C", str(engine), "rev-parse", "HEAD"], capture=True).strip()
    changed = run(["git", "-C", str(engine), "diff", "--name-only", "HEAD", "--"], capture=True).strip()
    if observed != entry["revision"]:
        if changed:
            raise ValueError("Engine source changed before pinned checkout; refusing to reset it")
        run(["git", "-C", str(engine), "checkout", "--detach", entry["revision"]])
    patch = repository / entry["patch"]
    if not changed:
        run(["git", "-C", str(engine), "apply", "--check", "--recount", str(patch)])
        run(["git", "-C", str(engine), "apply", "--recount", str(patch)])
    # A reverse check distinguishes the maintained patch from arbitrary edits.
    run(["git", "-C", str(engine), "apply", "--reverse", "--check", "--recount", str(patch)])
    expected_changes = {line.removeprefix("+++ b/") for line in patch.read_text().splitlines()
                        if line.startswith("+++ b/")}
    actual_changes = set(run(["git", "-C", str(engine), "diff", "--name-only", "HEAD", "--"], capture=True).splitlines())
    untracked = run(["git", "-C", str(engine), "ls-files", "--others", "--", "*.py", "*.pyc", "*.so", "*.dylib", "*.pyd"], capture=True)
    if actual_changes != expected_changes or untracked:
        raise ValueError("Engine source contains changes outside the maintained patch")
    for relative, digest in requirements["sources"].items():
        checked_file(engine / relative, digest, f"Frozen source {relative}")


def publish_launcher(candidate: Path, directory: Path, *, replace: bool = False) -> Path:
    """Atomic publication on the destination filesystem, with no implicit replacement."""
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / candidate.name
    if destination.is_symlink():
        raise ValueError("Refusing to replace a symlink launcher")
    if destination.is_file() and destination.read_bytes() == candidate.read_bytes():
        metadata = destination.stat()
        if metadata.st_uid == os.getuid() and stat.S_IMODE(metadata.st_mode) == 0o700:
            return destination
    fd, temporary_name = tempfile.mkstemp(prefix=".dvbfixer-runner-", dir=directory)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(candidate.read_bytes())
            stream.flush()
            os.fchmod(stream.fileno(), 0o700)
            os.fsync(stream.fileno())
        if replace:
            temporary.replace(destination)
        else:
            os.link(temporary, destination)
        return destination
    finally:
        temporary.unlink(missing_ok=True)


def copy_checked_asset(source: Path, destination: Path, digest: str) -> None:
    """A partial copy never becomes the resume-visible destination."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=".dvbfixer-asset-", dir=destination.parent)
    temporary = Path(temporary_name)
    os.close(fd)
    try:
        shutil.copyfile(source, temporary)
        checked_file(temporary, digest, "Copied Protenix asset")
        os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def provision_psi4(*, entry: dict[str, Any], prefix: Path, manager: Manager,
                   repository: Path, resume: bool, dry_run: bool) -> dict[str, Any]:
    if platform_subdir() not in entry["platforms"]:
        raise ValueError(f"PSI4 runtime does not support {platform_subdir()}")
    if prefix.is_symlink():
        raise ValueError("Runtime prefix must not be a symlink")
    prefix = prefix.expanduser().resolve()
    if any(character.isspace() for character in str(prefix)):
        raise ValueError("Choose a runtime prefix without whitespace for native scripts")
    spec = repository / entry["spec"]
    report: dict[str, Any] = {"schema_version": 1, "runtime": "psi4", "state": "planned", "manager": manager.validate(),
        "prefix": str(prefix), "environment_prefix": str(prefix / "environment"),
        "spec": entry["spec"], "dependencies": entry["dependencies"], "platform": platform_subdir(),
        "spec_sha256": sha256(spec), "scientific_acceptance": "not-tested", "steps": []}
    if dry_run:
        return report
    reserve_prefix(prefix, "runtime:psi4", (), resume)
    prefix.mkdir(exist_ok=True)
    plan = {"runtime": "psi4", "manager": manager.name, "spec_sha256": sha256(spec)}
    plan_file = owned_child(prefix, "runtime-plan.json")
    if plan_file.exists():
        if json.loads(plan_file.read_text()) != plan:
            raise ValueError("PSI4 recipe/manager changed; choose a new prefix")
    else:
        with plan_file.open("x") as stream:
            json.dump(plan, stream, indent=2)
    evidence = owned_child(prefix, "evidence")
    evidence.mkdir(exist_ok=True)
    environment = owned_child(prefix, "environment")
    command = manager.update(environment, spec) if (environment / "conda-meta/history").is_file() else manager.create(environment, spec)
    code = ("import psi4,psiresp,math,json,platform; import importlib.metadata as m; "
            "psi4.core.be_quiet(); psi4.set_num_threads(1); psi4.set_memory('1 GB'); "
            "mol=psi4.geometry('0 1\\nH 0 0 0\\nH 0 0 0.74\\nsymmetry c1'); "
            "energy=psi4.energy('hf/sto-3g',molecule=mol); assert math.isfinite(energy) and energy<0; "
            "assert hasattr(psiresp.configs,'TwoStageRESP'); "
            "print(json.dumps({'python':platform.python_version(),'energy':energy,"
            "'packages':{d.metadata['Name']:d.version for d in m.distributions()}},sort_keys=True))")

    def save() -> None:
        temporary = evidence / f".installation-{uuid.uuid4().hex}.json"
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(evidence / "installation.json")

    try:
        for label, argv in (("install-environment", command),
                            ("native-qm-kernel", manager.run(environment, ["python", "-c", code])),
                            ("versions-environment", manager.run(environment, ["python", "-c", SNAPSHOT_CODE]))):
            log = evidence / f"{label}-{uuid.uuid4().hex}.log"
            step: dict[str, Any] = {"name": label, "command": argv, "log": log.name, "state": "running"}
            report["steps"].append(step)
            save()
            with log.open("w") as stream:
                result = subprocess.run(argv, env=child_environment(), stdout=stream,
                                        stderr=subprocess.STDOUT, check=False, timeout=3600)
            step.update(state="passed" if result.returncode == 0 else "failed", returncode=result.returncode)
            if result.returncode:
                raise ValueError(f"{label} failed; see {log}")
        report.update(state="ready", action=f"Select --qm-engine psi4 --psi4-env {environment}")
        save()
        return report
    except BaseException as exc:
        report.update(state="error", error=f"{type(exc).__name__}: {exc}")
        save()
        raise


def provision_runtime(*, name: str, prefix: Path, manager: Manager, repository: Path,
                      resume: bool = False, dry_run: bool = False,
                      source: Path | None = None, checkpoint: Path | None = None,
                      kalign: Path | None = None, auxiliary: Path | None = None,
                      runner_directory: Path | None = None, replace_runner: bool = False,
                      timeout: float = 120.0) -> dict[str, Any]:
    entry = load_catalog()["external"][name]
    if entry["kind"] == "separate-environment":
        if any((source, checkpoint, kalign, auxiliary, runner_directory, replace_runner)):
            raise ValueError("PSI4 provisioning has no engine/checkpoint/launcher inputs")
        return provision_psi4(entry=entry, prefix=prefix, manager=manager, repository=repository,
                              resume=resume, dry_run=dry_run)
    if platform_subdir() not in entry["platforms"]:
        raise ValueError(f"{name} requires {entry['platforms']}; actual host: {platform_subdir()}")
    if prefix.is_symlink():
        raise ValueError("Runtime prefix must not be a symlink")
    prefix = prefix.expanduser().resolve()
    if any(character.isspace() for character in str(prefix)):
        raise ValueError("Choose a runtime prefix without whitespace for native environment scripts")
    requirements = asset_requirements(entry, repository)
    if "patch" in requirements:
        checked_file(repository / entry["patch"], requirements["patch"], "Maintained patch")
    checkpoint = checkpoint.expanduser().absolute() if checkpoint else None
    kalign = kalign.expanduser().absolute() if kalign else None
    auxiliary = auxiliary.expanduser().absolute() if auxiliary else None
    missing = verify_assets(requirements, checkpoint=checkpoint, kalign=kalign, auxiliary=auxiliary)
    version = manager.validate()
    report: dict[str, Any] = {"schema_version": 1, "runtime": name, "profile": entry["profile"],
        "prefix": str(prefix), "manager": version, "state": "planned", "actions": missing,
        "scientific_acceptance": "not-tested", "platform": platform_subdir(), "steps": [],
        "sampling_spec": entry["sampling_spec"], "refinement_spec": entry.get("refinement_spec", "in-runtime"),
        "engine_revision": entry["revision"], "expected_assets": requirements}
    if dry_run:
        return report
    reserve_prefix(prefix, f"runtime:{name}", (), resume)
    prefix.mkdir(exist_ok=True)
    engine = owned_child(prefix, "engine")
    evidence = owned_child(prefix, "evidence")
    evidence.mkdir(exist_ok=True)
    env = runtime_environment(entry, repository, engine)
    if checkpoint and "config" in requirements:
        env["PROTPARDELLE_MODEL_PARAMS"] = str(checkpoint.parent.parent)
    roles = {"sampling": entry["sampling_spec"]}
    if "refinement_spec" in entry:
        roles["refinement"] = entry["refinement_spec"]
    plan = {"runtime": name, "manager": manager.name, "source_revision": entry["revision"],
            "patch_sha256": sha256(repository / entry["patch"]),
            "spec_sha256": {role: sha256(repository / spec) for role, spec in roles.items()}}
    plan_file = owned_child(prefix, "runtime-plan.json")
    if plan_file.exists():
        if json.loads(plan_file.read_text()) != plan:
            raise ValueError("Runtime composition/manager changed; choose a new prefix rather than mutating the frozen install")
    else:
        with plan_file.open("x") as stream:
            json.dump(plan, stream, indent=2)

    def save() -> None:
        target = evidence / "installation.json"
        temporary = evidence / f".installation-{uuid.uuid4().hex}.json"
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(target)

    def run(label: str, argv: list[str], *, extra_env: dict[str, str] | None = None) -> None:
        log = evidence / f"{label}-{uuid.uuid4().hex}.log"
        step = {"name": label, "command": argv, "log": log.name, "state": "running"}
        report["steps"].append(step)
        save()
        with log.open("w") as stream:
            result = subprocess.run(argv, env=extra_env or env, stdout=stream,
                                    stderr=subprocess.STDOUT, check=False, timeout=3600)
        step["returncode"] = result.returncode
        step["state"] = "passed" if result.returncode == 0 else "failed"
        save()
        if result.returncode:
            raise ValueError(f"{label} failed; see {log}")

    try:
        prepare_source(entry, repository, engine, source, evidence / "source.log", env, requirements)
        report["source_revision"] = entry["revision"]
        report["patch_sha256"] = sha256(repository / entry["patch"])
        for role, relative_spec in roles.items():
            target = owned_child(prefix, role)
            spec = repository / relative_spec
            report.setdefault("spec_sha256", {})[role] = sha256(spec)
            command = manager.update(target, spec) if (target / "conda-meta/history").is_file() else manager.create(target, spec)
            run(f"install-{role}", command)
        sampling = prefix / "sampling"
        device = entry["device"]
        if auxiliary:
            for relative, digest in requirements["auxiliary"].items():
                destination = owned_child(engine, relative)
                if destination.exists():
                    checked_file(destination, digest, f"Existing data {relative}")
                else:
                    copy_checked_asset(auxiliary / relative, destination, digest)
        # Upstream Protpardelle imports validate the model_params directory.
        # Without operator assets we can verify kernels, but must defer that import.
        engine_import = f"import {entry['engine_import']}; " if not missing else ""
        report["engine_import_checked"] = False
        code = ("from dvbfixer.model.diffusion.validate import validate_runner_result; "
                f"import torch; assert torch.__version__ == {entry['torch_version']!r}; "
                + engine_import +
                f"x = torch.ones((2, 2), device={device!r}); assert float((x @ x).sum().item()) == 8; "
                "print(torch.__version__, x.device)")
        run("sampling-kernels", manager.run(sampling, ["python", "-c", code]))
        report["engine_import_checked"] = not bool(missing)
        if "refinement_spec" in entry:
            run("refinement-kernels", manager.run(prefix / "refinement", ["python", str(repository / entry["runner"]),
                "--profile", entry["profile"], "--refinement-preflight"]))
        else:
            run("refinement-kernels", manager.run(sampling, ["python", "-c",
                "import openmm; s=openmm.System(); s.addParticle(12); "
                "f=openmm.CustomExternalForce('x*x+y*y+z*z'); f.addParticle(0,[]); s.addForce(f); "
                "i=openmm.VerletIntegrator(0.001); "
                "c=openmm.Context(s,i,openmm.Platform.getPlatformByName('Reference')); "
                "c.setPositions([[1,1,1]]); "
                "assert c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)==3"]))
        for role in roles:
            run(f"versions-{role}", manager.run(prefix / role, ["python", "-c", SNAPSHOT_CODE]))
        if missing:
            report["state"] = "action-required"
            save()
            return report
        staging = owned_child(prefix, "candidates") / uuid.uuid4().hex
        command = ["python", str(repository / "scripts/install_diffusion_backend.py"), entry["backend"],
                   str(checkpoint), "--engine-root", str(engine), "--runner-directory", str(staging)]
        if "refinement_spec" in entry:
            command.extend(["--refinement-python", str(prefix / "refinement/bin/python")])
        if kalign:
            command.extend(["--kalign", str(kalign)])
        run("stage-launcher", manager.run(sampling, command))
        candidates = list(staging.iterdir())
        if (len(candidates) != 1 or candidates[0].name != f"dvbfixer-diffusion-{entry['profile']}"
                or candidates[0].is_symlink() or not candidates[0].is_file()):
            raise ValueError("Registration did not create exactly one staged launcher")
        # The package's diffusion namespace imports NumPy. Run the canonical
        # handshake in the new backend, keeping this bootstrap stdlib-only.
        code = ("import json; from dvbfixer.model.diffusion.preflight import invoke_runner_preflight; "
                "from dvbfixer.model.diffusion.runtime import diffusion_runner_directory; "
                f"r=invoke_runner_preflight(profile={entry['profile']!r}, runner={str(candidates[0])!r}, "
                f"timeout_seconds={timeout!r}); "
                "print(json.dumps({'preflight':r.to_dict(),'runner_directory':str(diffusion_runner_directory())}))")
        probe = subprocess.run(manager.run(sampling, ["python", "-c", code]),
                               capture_output=True, text=True, timeout=timeout + 60, env=env, check=True)
        payload = json.loads(probe.stdout)
        report["preflight"] = payload["preflight"]
        if not report["preflight"]["passed"]:
            raise ValueError("Frozen runner preflight failed; no launcher was published")
        destination = publish_launcher(candidates[0], runner_directory or Path(payload["runner_directory"]), replace=replace_runner)
        report.update(state="ready", launcher=str(destination), launcher_sha256=sha256(destination), actions=[])
        save()
        return report
    except BaseException as exc:
        report.update(state="error", error=f"{type(exc).__name__}: {exc}")
        save()
        raise
