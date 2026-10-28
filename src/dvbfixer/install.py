"""Small, stdlib-only managed installation helper (also usable before bootstrap)."""

from __future__ import annotations

import argparse
import getpass
import json
import math
import os
import platform
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import uuid
from dataclasses import dataclass
from pathlib import Path

from dvbfixer.installation_profiles import composition, environment_spec, load_catalog


@dataclass(frozen=True)
class Manager:
    name: str
    executable: str

    def validate(self) -> str:
        result = subprocess.run([self.executable, "--version"], capture_output=True,
                                text=True, check=True, timeout=30)
        version = result.stdout.strip()
        match = re.search(r"(\d+)\.(\d+)", version)
        minimum = (23, 10) if self.name == "conda" else (2, 0)
        if match is None or tuple(map(int, match.groups())) < minimum:
            raise ValueError(f"{self.name} requires version {minimum[0]}.{minimum[1]} or newer")
        return version

    def create(self, prefix: Path, spec: Path) -> list[str]:
        subcommand = ["create"] if self.name == "micromamba" else ["env", "create"]
        return [self.executable, *subcommand, "--yes", "--prefix", str(prefix), "--file", str(spec)]

    def run(self, prefix: Path, command: list[str]) -> list[str]:
        # conda's default capture layer disconnects stdin, including the secret
        # pipe used for license setup. Capture, when needed, belongs to our caller.
        flags = ["--no-capture-output"] if self.name == "conda" else []
        return [self.executable, "run", *flags, "--prefix", str(prefix), *command]

    def update(self, prefix: Path, spec: Path) -> list[str]:
        subcommand = ["install", "--yes"] if self.name == "micromamba" else ["env", "update"]
        return [self.executable, *subcommand, "--prefix", str(prefix), "--file", str(spec)]


def child_environment() -> dict[str, str]:
    """Never leak another editable source, user packages or MODELLER key into a child."""
    env = os.environ.copy()
    for key in ("PYTHONPATH", "PYTHONHOME", "KEY_MODELLER"):
        env.pop(key, None)
    env["PYTHONNOUSERSITE"] = "1"
    env["CONDA_CHANNEL_PRIORITY"] = "strict"
    env["CONDA_ALWAYS_YES"] = "true"
    return env


def platform_subdir() -> str:
    value = (platform.system(), platform.machine())
    return {("Linux", "x86_64"): "linux-64", ("Darwin", "arm64"): "osx-arm64"}.get(
        value, "/".join(value)
    )


def reserve_prefix(prefix: Path, profile: str, addons: tuple[str, ...], resume: bool) -> Path:
    """Refuse existing environments; resume only installations recorded by this helper."""
    if prefix.is_symlink():
        raise ValueError("The target prefix must not be a symlink")
    if not prefix.parent.is_dir():
        raise ValueError(f"Create the target parent directory first: {prefix.parent}")
    marker = prefix.with_name(prefix.name + ".dvbfixer-install.json")
    expected = {"schema_version": 1, "prefix": str(prefix), "profile": profile,
                "addons": sorted(set(addons))}
    if resume:
        if marker.is_symlink() or not marker.is_file():
            raise ValueError("Cannot resume: no installation record owned by this helper")
        if json.loads(marker.read_text()) != expected:
            raise ValueError("Cannot resume: prefix/profile/add-ons differ from installation record")
    else:
        if prefix.exists():
            raise ValueError(f"Target already exists; choose a NEW prefix: {prefix}")
        with marker.open("x", encoding="utf-8") as stream:
            json.dump(expected, stream, indent=2)
    return marker


def validate_prefix_path(prefix: Path) -> None:
    if any(character.isspace() for character in str(prefix)):
        raise ValueError("Environment prefixes containing whitespace are unsupported by MAFFT's bundled "
                         "shell scripts; choose a new prefix without spaces. Input/source paths may contain spaces.")


MODELLER_SETUP = r'''
import ast, os, pathlib, sys, tempfile
prefix = pathlib.Path(sys.prefix).resolve()
configs = list(prefix.glob("lib/modeller-*/modlib/modeller/config.py"))
if len(configs) != 1:
    raise SystemExit(2)
path = configs[0]
if not path.resolve().is_relative_to(prefix) or path.is_symlink():
    raise SystemExit(2)
key = sys.stdin.read().strip()
if not key or "\n" in key or "\r" in key:
    raise SystemExit(2)
text = path.read_text()
tree = ast.parse(text)
assignments = [n for n in tree.body if isinstance(n, ast.Assign)
               and any(isinstance(t, ast.Name) and t.id == "license" for t in n.targets)]
if len(assignments) != 1:
    raise SystemExit(2)
node = assignments[0]
lines = text.splitlines(keepends=True)
lines[node.lineno-1:node.end_lineno] = ["license = " + repr(key) + "\n"]
fd, temp = tempfile.mkstemp(dir=path.parent)
try:
    os.fchmod(fd, 0o600)
    with os.fdopen(fd, "w") as stream:
        stream.write("".join(lines))
    os.replace(temp, path)
finally:
    if os.path.exists(temp):
        os.unlink(temp)
'''


def configure_modeller(manager: Manager, prefix: Path, key: str) -> None:
    # Secret travels only over stdin, never argv or inherited environment.
    result = subprocess.run(manager.run(prefix, ["python", "-c", MODELLER_SETUP]),
                            input=key, text=True, capture_output=True, timeout=30,
                            env=child_environment(), check=False)
    if result.returncode:
        raise ValueError("MODELLER configuration failed in target environment (details suppressed)")
    result = subprocess.run(manager.run(prefix, ["python", "-c", "import modeller"]),
                            capture_output=True, timeout=30, env=child_environment(), check=False)
    if result.returncode:
        raise ValueError("MODELLER license/import check failed (details suppressed)")


def run_checked(command: list[str], *, cwd: Path | None = None) -> None:
    subprocess.run(command, cwd=cwd, env=child_environment(), check=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("standard", "desktop"))
    parser.add_argument("--manager", choices=("micromamba", "mamba", "conda"))
    parser.add_argument("--prefix", type=Path)
    parser.add_argument("--addon", action="append", default=[], choices=tuple(load_catalog()["addons"]))
    parser.add_argument("--non-interactive", action="store_true")
    parser.add_argument("--resume", action="store_true", help="Resume only a helper-owned installation")
    parser.add_argument("--configure-only", action="store_true",
                        help="With --resume, finish license setup without reinstalling packages or rebuilding GUI")
    parser.add_argument("--modeller-key-fd", type=int,
                        help="Read license from an already-open protected file descriptor")
    parser.add_argument("--skip-modeller", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Print composition without creating files")
    from dvbfixer.runtime_install import runtime_names
    parser.add_argument("--runtime", choices=runtime_names(),
                        help="Provision a separate diffusion runtime in --prefix (without --profile)")
    parser.add_argument("--engine-source", type=Path, help="Optional local Git clone source (read-only)")
    parser.add_argument("--checkpoint", type=Path, help="Operator-provided frozen weights; never downloaded")
    parser.add_argument("--kalign", type=Path, help="Operator-provided digest-pinned Protenix Kalign executable")
    parser.add_argument("--auxiliary-root", type=Path, help="Protenix data root containing frozen common/ assets")
    parser.add_argument("--runner-directory", type=Path, help="Publish the checked launcher to this explicit directory")
    parser.add_argument("--replace-runner", action="store_true", help="Explicitly replace an existing launcher after preflight")
    parser.add_argument("--preflight-timeout", type=float, default=120.0)
    args = parser.parse_args(argv)
    if args.modeller_key_fd is not None and (args.skip_modeller or args.modeller_key_fd < 3):
        parser.error("Use --modeller-key-fd >=3, without --skip-modeller")
    if args.configure_only and not args.resume:
        parser.error("--configure-only requires --resume")
    interactive = not args.non_interactive and sys.stdin.isatty()
    if interactive:
        if not args.profile and not args.runtime:
            choice = input(f"Profile [standard/desktop] or runtime [{'/'.join(runtime_names())}]: ").strip()
            if choice in runtime_names():
                args.runtime = choice
            else:
                args.profile = choice
        args.manager = args.manager or input("Manager [micromamba/mamba/conda]: ").strip()
        args.prefix = args.prefix or Path(input("NEW environment prefix: ").strip())
        if not args.addon and not args.runtime:
            selected = input("Add-ons (comma-separated msa-alternatives,probe; empty for none): ").strip()
            args.addon = [value.strip() for value in selected.split(",") if value.strip()]
        if args.runtime:
            if not args.checkpoint:
                selected = input("Operator checkpoint path (empty to configure later): ").strip()
                args.checkpoint = Path(selected) if selected else None
            if args.runtime == "protenix-cuda":
                if not args.kalign:
                    selected = input("Pinned Kalign executable (empty to configure later): ").strip()
                    args.kalign = Path(selected) if selected else None
                if not args.auxiliary_root:
                    selected = input("Protenix data root containing common/ (empty to configure later): ").strip()
                    args.auxiliary_root = Path(selected) if selected else None
    if args.manager is not None and args.manager not in ("micromamba", "mamba", "conda"):
        parser.error("Unknown manager")
    if not args.runtime and any((args.engine_source, args.checkpoint, args.kalign, args.auxiliary_root,
                                 args.runner_directory, args.replace_runner)):
        parser.error("Runtime source/checkpoint/registration options require --runtime")
    if args.runtime:
        if args.profile or args.addon or args.configure_only or args.modeller_key_fd is not None:
            parser.error("Runtime provisioning uses its own --prefix, without --profile/--addon/license setup")
        if not args.manager or not args.prefix:
            parser.error("Runtime provisioning requires explicit --manager and --prefix")
        if not math.isfinite(args.preflight_timeout) or args.preflight_timeout <= 0:
            parser.error("--preflight-timeout must be positive and finite")
        try:
            from dvbfixer.runtime_install import provision_runtime
            executable = shutil.which(args.manager)
            if executable is None:
                raise ValueError(f"Install the selected manager first: {args.manager}")
            repository = Path(__file__).resolve().parents[2]
            if not (repository / "scripts/install_diffusion_backend.py").is_file():
                raise ValueError("Run runtime provisioning from a source checkout containing deploy/ and scripts/")
            report = provision_runtime(name=args.runtime, prefix=args.prefix,
                manager=Manager(args.manager, executable), repository=repository,
                resume=args.resume, dry_run=args.dry_run, source=args.engine_source,
                checkpoint=args.checkpoint, kalign=args.kalign, auxiliary=args.auxiliary_root,
                runner_directory=args.runner_directory, replace_runner=args.replace_runner,
                timeout=args.preflight_timeout)
            print(json.dumps(report, indent=2))
            return 2 if report["state"] == "action-required" else 0
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
    if not all((args.profile, args.manager, args.prefix)):
        parser.error("Explicit --profile, --manager and --prefix are required without an interactive terminal")
    if args.manager not in ("micromamba", "mamba", "conda"):
        parser.error("Unknown manager")
    try:
        packages = composition(args.profile, tuple(args.addon))
        if platform_subdir() not in load_catalog()["platforms"]:
            raise ValueError(f"Unsupported platform: {platform_subdir()}; use a supported Linux/macOS host")
        prefix = args.prefix.expanduser().absolute()
        if prefix.is_symlink():
            raise ValueError("The target prefix must not be a symlink")
        prefix = prefix.resolve()
        validate_prefix_path(prefix)
        executable = shutil.which(args.manager)
        if executable is None:
            raise ValueError(f"Install the selected manager first: {args.manager}")
        manager = Manager(args.manager, executable)
        version = manager.validate()
        print(json.dumps({"profile": args.profile, "addons": args.addon, "prefix": str(prefix),
                          "manager": version, "packages": [p["spec"] for p in packages]}, indent=2), flush=True)
        if args.dry_run:
            return 0
        reserve_prefix(prefix, args.profile, tuple(args.addon), args.resume)
        root = Path(__file__).resolve().parents[2]
        if not (root / "pyproject.toml").is_file():
            raise ValueError("Run the helper from a source checkout containing pyproject.toml and gui/")
        # Manager installs native dependencies; pip uses normal resolution afterwards.
        if not args.configure_only:
            with tempfile.TemporaryDirectory(prefix="dvbfixer-install-") as temp:
                spec = Path(temp) / "environment.yml"
                spec.write_text(environment_spec(args.profile, tuple(args.addon)).split("  - pip:\n")[0])
                if (prefix / "conda-meta/history").is_file():
                    # History may exist after an interrupted native transaction.
                    run_checked(manager.update(prefix, spec))
                else:
                    run_checked(manager.create(prefix, spec))
            pip_specs = [p["spec"] for p in packages if p.get("source") == "pip"]
            run_checked(manager.run(prefix, ["python", "-m", "pip", "install", *pip_specs, str(root)]))
        if args.modeller_key_fd is not None:
            metadata = os.fstat(args.modeller_key_fd)
            if stat.S_ISREG(metadata.st_mode) and (
                stat.S_IMODE(metadata.st_mode) & 0o077 or metadata.st_uid != os.getuid()
            ):
                raise ValueError("MODELLER key file must be owned by this user with permissions 0600")
            with os.fdopen(os.dup(args.modeller_key_fd), encoding="utf-8") as stream:
                key = stream.read(4097)
            if len(key) > 4096:
                raise ValueError("MODELLER secret exceeds supported length")
            configure_modeller(manager, prefix, key)
        elif interactive and not args.skip_modeller:
            key = getpass.getpass("MODELLER key (empty to configure later): ")
            if key:
                configure_modeller(manager, prefix, key)
        if args.profile == "desktop" and not args.configure_only:
            # Isolate GUI artifacts from the checkout and other installations.
            gui = prefix / "share/dvbfixer/gui"
            gui.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(root / "gui", gui, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns("node_modules", "dist", "dist-server", "structures", ".env*"))
            run_checked(manager.run(prefix, ["npm", "ci"]), cwd=gui)
            run_checked(manager.run(prefix, ["npm", "run", "build"]), cwd=gui)
        run_checked(manager.run(prefix, ["python", str(root / "scripts/check_scientific_linalg.py")]))
        license_probe = subprocess.run(manager.run(prefix, ["python", "-c", "import modeller"]),
                                       env=child_environment(), capture_output=True, timeout=30, check=False)
        if not args.configure_only or license_probe.returncode == 0:
            evidence = prefix / "share/dvbfixer/installation-smoke" / uuid.uuid4().hex
            smoke = ["python", str(root / "scripts/smoke_installation_profile.py"),
                     "--output", str(evidence), "--root", str(root)]
            if license_probe.returncode == 0:
                smoke.append("--with-modeller")
            if args.configure_only:
                smoke.extend(["--only", "model", "--only", "homology"])
            run_checked(manager.run(prefix, smoke))
            print(f"Functional smoke evidence: {evidence}", flush=True)
        result = subprocess.run(manager.run(prefix, ["dvbfixer", "doctor", "--profile", args.profile,
                                *[arg for addon in args.addon for arg in ("--addon", addon)]]),
                                env=child_environment(), check=False)
        return result.returncode
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
