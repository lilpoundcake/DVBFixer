"""Small, stdlib-only managed installation helper (also usable before bootstrap)."""

from __future__ import annotations

import argparse
import getpass
import json
import math
import os
import platform
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import uuid
from dataclasses import dataclass
from pathlib import Path

from dvbfixer.install_progress import logged_run, stage, waiting_for_input
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
        # libmamba otherwise merges channels from user rc files even when
        # the environment spec declares its own channels and nodefaults.
        flags = ["--override-channels"] if self.name in ("micromamba", "mamba") else []
        return [self.executable, *subcommand, *flags, "--yes", "--prefix", str(prefix), "--file", str(spec)]

    def run(self, prefix: Path, command: list[str]) -> list[str]:
        # conda's default capture layer disconnects stdin, including the secret
        # pipe used for license setup. Capture, when needed, belongs to our caller.
        flags = ["--no-capture-output"] if self.name == "conda" else []
        return [self.executable, "run", *flags, "--prefix", str(prefix), *command]

    def update(self, prefix: Path, spec: Path) -> list[str]:
        subcommand = ["install", "--yes"] if self.name == "micromamba" else ["env", "update"]
        flags = ["--override-channels"] if self.name in ("micromamba", "mamba") else []
        return [self.executable, *subcommand, *flags, "--prefix", str(prefix), "--file", str(spec)]


def child_environment() -> dict[str, str]:
    """Never leak another editable source, user packages or MODELLER key into a child."""
    env = os.environ.copy()
    for key in ("PYTHONPATH", "PYTHONHOME", "KEY_MODELLER"):
        env.pop(key, None)
    # A caller may prepend a tools prefix to PATH while another conda env is
    # active. Keeping that activation state makes conda run replace the old
    # prefix in-place, leaving tools/python ahead of the requested interpreter.
    # Start child activation from an inactive shell; never mutate the caller.
    for key in tuple(env):
        if key in ("CONDA_PREFIX", "CONDA_SHLVL", "CONDA_DEFAULT_ENV", "CONDA_PROMPT_MODIFIER") or key.startswith("CONDA_PREFIX_"):
            env.pop(key, None)
    env["PYTHONNOUSERSITE"] = "1"
    env["CONDA_CHANNEL_PRIORITY"] = "strict"
    env["CONDA_ALWAYS_YES"] = "true"
    # Conda must not report ready after a corrupted cache or conflicting package
    # files were accepted by a permissive caller policy.
    env["CONDA_SAFETY_CHECKS"] = "enabled"
    env["CONDA_PATH_CONFLICT"] = "prevent"
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

MODELLER_READ = r'''
import ast, pathlib, sys
prefix = pathlib.Path(sys.prefix).resolve()
configs = list(prefix.glob('lib/modeller-*/modlib/modeller/config.py'))
if len(configs) != 1 or configs[0].is_symlink() or not configs[0].resolve().is_relative_to(prefix):
    raise SystemExit(1)
tree = ast.parse(configs[0].read_text())
values = [ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
          and any(isinstance(t, ast.Name) and t.id == 'license' for t in n.targets)]
if len(values) != 1 or not isinstance(values[0], str) or not values[0] or '\n' in values[0] or '\r' in values[0]:
    raise SystemExit(1)
sys.stdout.write(values[0])
'''


def configure_modeller(manager: Manager, prefix: Path, key: str, *, verify_import: bool = True) -> None:
    # Secret travels only over stdin, never argv or inherited environment.
    result = subprocess.run(manager.run(prefix, ["python", "-c", MODELLER_SETUP]),
                            input=key, text=True, capture_output=True, timeout=30,
                            env=child_environment(), check=False)
    if result.returncode:
        raise ValueError("MODELLER configuration failed in target environment (details suppressed)")
    if not verify_import:
        return
    result = subprocess.run(manager.run(prefix, ["python", "-c", "import modeller"]),
                            capture_output=True, timeout=30, env=child_environment(), check=False)
    if result.returncode:
        raise ValueError("MODELLER license/import check failed (details suppressed)")


def run_checked(command: list[str], *, cwd: Path | None = None, log: Path | None = None) -> None:
    if log is None:
        subprocess.run(command, cwd=cwd, env=child_environment(), check=True)
    else:
        logged_run(command, cwd=cwd, env=child_environment(), log=log)


def recorded_installation(prefix: Path) -> dict:
    """Read ownership without adopting an arbitrary conda environment."""
    prefix = prefix.expanduser().absolute()
    if prefix.is_symlink():
        raise ValueError("The target prefix must not be a symlink")
    marker = prefix.with_name(prefix.name + ".dvbfixer-install.json")
    if marker.is_symlink() or not marker.is_file():
        raise ValueError("No helper-owned installation record; for a manual installation choose package-only update "
                         "(--update --package-only) or explicitly attach/update a profile with --extend --profile")
    record = json.loads(marker.read_text())
    if (not isinstance(record, dict) or record.get("schema_version") != 1
            or record.get("prefix") != str(prefix.resolve())
            or not isinstance(record.get("profile"), str)
            or not isinstance(record.get("addons"), list)
            or not all(isinstance(value, str) for value in record["addons"])):
        raise ValueError("Invalid installation record; preserve it for inspection")
    return record


def runtime_summary(report: dict) -> None:
    print(f"\nBackend: {report['runtime']}\nInstallation directory: {report['prefix']}", flush=True)
    if report["state"] == "planned":
        print("Dry run only; no installation was changed.", flush=True)
        return
    print(f"Detailed report: {Path(report['prefix']) / 'evidence/installation.json'}", flush=True)
    if report["state"] == "action-required":
        print("ACTION REQUIRED: dependencies installed, but the backend is not ready.", flush=True)
        for action in report.get("actions", []):
            print(f"  - {action}", flush=True)
        print("Resume this same backend and directory to resolve the listed actions. Model downloads require explicit consent.", flush=True)
    else:
        print("Backend provisioning checks passed. This does not establish scientific qualification.", flush=True)
        if report.get("action"):
            print(report["action"], flush=True)
    print(f"Scientific acceptance: {report['scientific_acceptance']}", flush=True)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("standard", "desktop"))
    parser.add_argument("--manager", choices=("micromamba", "mamba", "conda"))
    target = parser.add_mutually_exclusive_group()
    target.add_argument("--name", help="Environment name; select or create it through the chosen manager")
    target.add_argument("--prefix", type=Path,
                        help="Advanced: explicit installation directory instead of name/number selection")
    parser.add_argument("--addon", action="append", default=[], choices=tuple(load_catalog()["addons"]))
    parser.add_argument("--non-interactive", action="store_true")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true", help="Continue an interrupted helper-owned installation")
    mode.add_argument("--update", action="store_true", help="Select and update an existing DVBFixer installation")
    mode.add_argument("--extend", action="store_true", help="Add Standard/Desktop or GUI to an existing environment after preview and confirmation")
    parser.add_argument("--gui-only", action="store_true", help="With --extend, add/update GUI without changing the scientific stack")
    parser.add_argument("--package-only", action="store_true",
                        help="With --update, reinstall only DVBFixer without changing scientific dependencies or adopting a manual environment")
    parser.add_argument("--json", action="store_true", help="Also print the machine-readable plan/result")
    parser.add_argument("--configure-only", action="store_true",
                        help="With --resume, finish license setup without reinstalling packages or rebuilding GUI")
    parser.add_argument("--modeller-key-fd", type=int,
                        help="Read license from an already-open protected file descriptor")
    parser.add_argument("--skip-modeller", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Preview without changing the target environment (solver caches/preview logs may be written)")
    parser.add_argument("--desktop-shortcut", action="store_true",
                        help="Create a per-user macOS application/Linux desktop entry (Desktop or GUI-only)")
    from dvbfixer.runtime_install import runtime_names
    parser.add_argument("--runtime", choices=runtime_names(),
                        help="Provision a separate diffusion or PSI4 runtime using --name or advanced --prefix (without --profile)")
    parser.add_argument("--engine-source", type=Path, help="Optional local Git clone source (read-only)")
    parser.add_argument("--checkpoint", type=Path, help="Existing digest-matching model weights")
    parser.add_argument("--download-weights", action="store_true", help="Download official pinned weights into this checkout's shared models/ cache")
    parser.add_argument("--accept-model-terms", action="store_true", help="Confirm reviewed model terms for non-interactive --download-weights")
    parser.add_argument("--kalign", type=Path, help="Operator-provided digest-pinned Protenix Kalign executable")
    parser.add_argument("--auxiliary-root", type=Path, help="Protenix data root containing frozen common/ assets")
    parser.add_argument("--runner-directory", type=Path, help="Publish the checked launcher to this explicit directory")
    parser.add_argument("--replace-runner", action="store_true", help="Explicitly replace an existing launcher after preflight")
    parser.add_argument("--preflight-timeout", type=float, default=120.0)
    args = parser.parse_args(argv)
    from dvbfixer import install_result

    install_result.configure(args)
    if args.modeller_key_fd is not None and (args.skip_modeller or args.modeller_key_fd < 3):
        parser.error("Use --modeller-key-fd >=3, without --skip-modeller")
    if args.configure_only and not args.resume:
        parser.error("--configure-only requires --resume")
    if args.gui_only and not args.extend:
        parser.error("--gui-only requires --extend")
    if args.gui_only and args.profile is not None:
        parser.error("Choose either --gui-only or --profile, not both")
    if args.extend and (args.runtime or args.package_only or args.configure_only):
        parser.error("--extend cannot be combined with a runtime, --package-only or --configure-only")
    interactive = not args.non_interactive and sys.stdin.isatty()
    args.named_parent = None
    from dvbfixer.install_discovery import resolve_selection

    try:
        resolve_selection(args, interactive=interactive, runtimes=runtime_names(),
                          repository=Path(__file__).resolve().parents[2])
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        install_result.failure(exc)
        print(f"error: {exc}", file=sys.stderr)
        return 1
    install_result.configure(args)
    if (args.download_weights or args.accept_model_terms) and not args.runtime:
        parser.error("Model download options require a diffusion --runtime")
    if args.accept_model_terms and not args.download_weights:
        parser.error("--accept-model-terms requires --download-weights")
    if interactive:
        from dvbfixer.install_ui import addons, ask, confirm

        if args.profile == "desktop" and not args.desktop_shortcut and not args.package_only:
            args.desktop_shortcut = confirm("Create a user application shortcut?")
        if not args.addon and not args.runtime and not (args.update or args.resume or args.gui_only or args.full_update_selected):
            args.addon = addons([(name, name) for name in load_catalog()["addons"]])
        if args.runtime:
            from dvbfixer.runtime_install import asset_instructions

            instructions = asset_instructions(args.runtime)
            if instructions:
                print("\nRequired model assets (local file or confirmed managed download):", flush=True)
                for instruction in instructions:
                    print(f"  - {instruction}", flush=True)
                print("Empty asset paths defer registration: dependencies can be installed, but the backend will not be ready.", flush=True)
            try:
                from dvbfixer.install_weights import select

                select(args, Path(__file__).resolve().parents[2], interactive=True)
            except (OSError, ValueError) as exc:
                install_result.failure(exc)
                print(f"error: {exc}", file=sys.stderr)
                return 1
            if args.runtime == "protenix-cuda":
                if not args.kalign:
                    selected = ask("Pinned Kalign executable (empty to configure later): ").strip()
                    args.kalign = Path(selected) if selected else None
                if not args.auxiliary_root:
                    selected = ask("Protenix data root containing common/ (empty to configure later): ").strip()
                    args.auxiliary_root = Path(selected) if selected else None
    elif args.runtime:
        try:
            from dvbfixer.install_weights import select

            select(args, Path(__file__).resolve().parents[2], interactive=False)
        except (OSError, ValueError) as exc:
            install_result.failure(exc)
            print(f"error: {exc}", file=sys.stderr)
            return 1
    if args.extend:
        if args.gui_only:
            args.profile = "desktop"
        if args.prefix is None or args.profile not in {"standard", "desktop"}:
            parser.error("--extend requires an existing environment and --profile or --gui-only outside an interactive terminal")
        try:
            from dvbfixer.install_extend import extension_composition

            args.profile, args.addon = extension_composition(args)
        except (OSError, ValueError) as exc:
            install_result.failure(exc)
            print(f"error: {exc}", file=sys.stderr)
            return 1
    elif args.update:
        if args.runtime or args.configure_only:
            parser.error("--update applies only to an existing Standard/Desktop installation")
        if args.package_only:
            from dvbfixer.install_update import update_package_only

            try:
                install_result.configure(args)
                return update_package_only(args, interactive=interactive)
            except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
                install_result.failure(exc)
                print(f"error: {exc}", file=sys.stderr)
                return 1
        if args.prefix is None:
            parser.error("--update requires --name or advanced --prefix outside an interactive terminal")
        try:
            record = recorded_installation(args.prefix)
            if record["profile"] not in ("standard", "desktop"):
                raise ValueError("Selected installation is a backend, not DVBFixer Standard/Desktop")
            if args.profile and args.profile != record["profile"]:
                raise ValueError("Update profile differs from the installation record")
            if args.addon and sorted(set(args.addon)) != record["addons"]:
                raise ValueError("Update must retain the recorded add-ons")
            args.profile, args.addon = record["profile"], record["addons"]
            if not (args.prefix.expanduser() / "conda-meta/history").is_file():
                raise ValueError("Existing environment is incomplete; use --resume, not --update")
        except (OSError, ValueError) as exc:
            install_result.failure(exc)
            print(f"error: {exc}", file=sys.stderr)
            return 1
    elif args.package_only:
        parser.error("--package-only requires --update")
    elif args.resume and args.prefix is not None:
        try:
            record = recorded_installation(args.prefix)
            if record["profile"].startswith("runtime:"):
                runtime = record["profile"].removeprefix("runtime:")
                if args.profile or (args.runtime and args.runtime != runtime):
                    raise ValueError("Selected backend differs from the existing installation")
                args.runtime = runtime
            else:
                if args.runtime or (args.profile and args.profile != record["profile"]):
                    raise ValueError("Selected profile differs from the existing installation")
                args.profile = record["profile"]
                if not args.addon:
                    args.addon = record["addons"]
        except (OSError, ValueError) as exc:
            install_result.failure(exc)
            print(f"error: {exc}", file=sys.stderr)
            return 1
    if args.manager is not None and args.manager not in ("micromamba", "mamba", "conda"):
        parser.error("Unknown manager")
    if args.runtime is not None and args.runtime not in runtime_names():
        parser.error("Unknown additional backend")
    if args.desktop_shortcut and (args.profile != "desktop" or args.runtime):
        parser.error("--desktop-shortcut requires --profile desktop, without --runtime")
    if not args.runtime and any((args.engine_source, args.checkpoint, args.kalign, args.auxiliary_root,
                                 args.runner_directory, args.replace_runner)):
        parser.error("Runtime source/checkpoint/registration options require --runtime")
    if args.runtime:
        install_result.configure(args)
        if args.profile or args.addon or args.configure_only or args.modeller_key_fd is not None:
            parser.error("Runtime provisioning uses its own --prefix, without --profile/--addon/license setup")
        if not args.manager or not args.prefix:
            parser.error("Runtime provisioning requires --manager and --name (or advanced --prefix)")
        if not math.isfinite(args.preflight_timeout) or args.preflight_timeout <= 0:
            parser.error("--preflight-timeout must be positive and finite")
        try:
            from dvbfixer.runtime_install import confirm_runtime_replacement, provision_runtime
            executable = shutil.which(args.manager)
            if executable is None:
                raise ValueError(f"Install the selected manager first: {args.manager}")
            repository = Path(__file__).resolve().parents[2]
            if not (repository / "scripts/install_diffusion_backend.py").is_file():
                raise ValueError("Run runtime provisioning from a source checkout containing deploy/ and scripts/")
            print(f"Additional backend: {args.runtime} (not a Standard/Desktop installation)", flush=True)
            if args.named_parent is not None and not args.dry_run:
                if platform_subdir() not in load_catalog()["external"][args.runtime]["platforms"]:
                    raise ValueError("Additional backend is unsupported on this platform")
                args.named_parent.mkdir(parents=True, exist_ok=True)
            with stage("Provision additional backend"):
                if not args.dry_run:
                    install_result.started()
                report = provision_runtime(name=args.runtime, prefix=args.prefix,
                    manager=Manager(args.manager, executable), repository=repository,
                    resume=args.resume, dry_run=args.dry_run, source=args.engine_source,
                    checkpoint=args.checkpoint, kalign=args.kalign, auxiliary=args.auxiliary_root,
                     runner_directory=args.runner_directory, replace_runner=args.replace_runner,
                    confirm_replace=confirm_runtime_replacement if interactive else None,
                     timeout=args.preflight_timeout)
            runtime_summary(report)
            install_result.required(report.get("actions", []))
            if args.json:
                print(json.dumps(report, indent=2))
            return 2 if report["state"] == "action-required" else 0
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            install_result.failure(exc)
            print(f"error: {exc}", file=sys.stderr)
            return 1
    if not all((args.profile, args.manager, args.prefix)):
        parser.error("Explicit --profile, --manager and --name (or advanced --prefix) are required without an interactive terminal")
    if args.manager not in ("micromamba", "mamba", "conda"):
        parser.error("Unknown manager")
    extension = None
    install_result.configure(args)
    try:
        packages = composition(args.profile, tuple(args.addon))
        if args.gui_only:
            packages = [package for package in packages if package["group"] == "gui"]
        catalog = load_catalog()
        subdir = platform_subdir()
        if subdir not in catalog["platforms"]:
            raise ValueError(f"Unsupported platform: {subdir}; use a supported Linux/macOS host")
        for addon in args.addon:
            entry = catalog["addons"][addon]
            if subdir not in entry.get("managed_platforms", catalog["platforms"]):
                route = entry.get("manual_route", "docs/installation.md")
                raise ValueError(f"Managed add-on {addon!r} is unavailable on {subdir}; "
                                 f"see {route}. Omit this helper add-on when using the manual route.")
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
        scope = "GUI only" if args.gui_only else args.profile
        operation = "Update profile" if args.full_update_selected else "Extend" if args.extend else "Update" if args.update else "Resume" if args.resume else "Install"
        print(f"\n{operation} DVBFixer {scope}\n"
              f"Directory: {prefix}\nManager: {args.manager} {version}\n"
              f"Add-ons: {', '.join(args.addon) or 'none'}", flush=True)
        if args.json:
            print(json.dumps({"profile": None if args.gui_only else args.profile, "gui_only": args.gui_only,
                              "addons": args.addon, "prefix": str(prefix),
                              "desktop_shortcut": args.desktop_shortcut,
                              "manager": version, "packages": [p["spec"] for p in packages]}, indent=2), flush=True)
        if args.update:
            with stage("Read installed version"):
                installed = subprocess.run(manager.run(prefix, ["python", "-c",
                    "from importlib.metadata import version; print(version('dvbfixer'))"]),
                    env=child_environment(), capture_output=True, text=True, check=True, timeout=30)
            from dvbfixer import __version__
            print(f"Version: {installed.stdout.strip()} -> {__version__} (current checkout)", flush=True)
            print("Close this installation's GUI and jobs before updating. Updates are in-place; no automatic rollback.", flush=True)
            from dvbfixer.install_ui import confirm

            if interactive and not confirm("Update this installation?"):
                install_result.cancelled()
                print("Cancelled; nothing changed.", flush=True)
                return 0
        if args.extend:
            from dvbfixer.install_extend import prepare_extension_interactive
            from dvbfixer.install_ui import confirm

            extension = prepare_extension_interactive(args, manager, interactive=interactive)
            if extension is None:
                install_result.cancelled()
                print("Cancelled; nothing changed.", flush=True)
                return 0
            if args.gui_only:
                install_result.configure(args)
                packages = [package for package in composition("desktop") if package["group"] == "gui"]
            if not args.dry_run and interactive and not confirm("Apply this profile/GUI plan to the selected environment?"):
                install_result.cancelled()
                print("Cancelled; nothing changed.", flush=True)
                return 0
        if args.dry_run:
            print("Dry run only; no installation was changed.", flush=True)
            return 0
        if args.named_parent is not None:
            args.named_parent.mkdir(parents=True, exist_ok=True)
        if not args.extend:
            reserve_prefix(prefix, args.profile, tuple(args.addon), args.resume or args.update)
        logs = prefix.parent / f"{prefix.name}.installation-logs" / uuid.uuid4().hex
        logs.mkdir(mode=0o700, parents=True, exist_ok=False)
        print(f"Installation logs: {logs}", flush=True)
        install_result.log_directory(logs)
        if extension is not None:
            extension.retain(logs)
            extension.verify_before_apply()

        def run(label: str, command: list[str], *, cwd: Path | None = None) -> None:
            run_checked(command, cwd=cwd, log=logs / f"{label}.log")

        saved_license = None
        saved_license_ready = False
        preserve_existing = args.update or (args.extend and not args.gui_only and bool(list(prefix.glob("lib/modeller-*/modlib/modeller/config.py"))))
        preserved_license = False
        if preserve_existing and args.modeller_key_fd is None:
            with stage("Preserve existing MODELLER configuration (protected; output not logged)"):
                saved = subprocess.run(manager.run(prefix, ["python", "-c", MODELLER_READ]),
                                       env=child_environment(), capture_output=True, text=True, timeout=30, check=False)
                if saved.returncode or not saved.stdout or len(saved.stdout) > 4096:
                    raise ValueError("Cannot safely preserve MODELLER configuration; no packages were updated")
                saved_license = saved.stdout
                preserved_license = True
                saved_license_ready = subprocess.run(manager.run(prefix, ["python", "-c", "import modeller"]),
                    env=child_environment(), capture_output=True, timeout=30, check=False).returncode == 0
        root = Path(__file__).resolve().parents[2]
        if not (root / "pyproject.toml").is_file():
            raise ValueError("Run the helper from a source checkout containing pyproject.toml and gui/")
        # Manager installs native dependencies; pip uses normal resolution afterwards.
        install_result.started()
        try:
            if not args.configure_only:
                if extension is not None:
                    run("Install confirmed native profile components", extension.native_command)
                    extension.verify_pip_after_native()
                else:
                    with tempfile.TemporaryDirectory(prefix="dvbfixer-install-") as temp:
                        spec = Path(temp) / "environment.yml"
                        spec.write_text(environment_spec(args.profile, tuple(args.addon)).split("  - pip:\n")[0])
                        if (prefix / "conda-meta/history").is_file():
                            # History may exist after an interrupted native transaction.
                            run("Update scientific dependencies", manager.update(prefix, spec))
                        else:
                            run("Install scientific dependencies", manager.create(prefix, spec))
                pip_specs = [p["spec"] for p in packages if p.get("source") == "pip"]
                pip_flags = ["--upgrade"] if args.update else []
                command = extension.pip_command if extension is not None else manager.run(prefix, ["python", "-m", "pip", "install", *pip_flags, *pip_specs, str(root)])
                if args.gui_only:
                    from dvbfixer.install_native_tools import preserve_native_openbabel_tool

                    with preserve_native_openbabel_tool(prefix):
                        run("Install DVBFixer package", command)
                else:
                    run("Install DVBFixer package", command)
        except BaseException:
            if saved_license is not None:
                try:
                    with stage("Recover MODELLER configuration after interrupted update (protected)"):
                        configure_modeller(manager, prefix, saved_license, verify_import=False)
                except (OSError, ValueError, subprocess.SubprocessError):
                    print("WARNING: Could not restore MODELLER configuration; keep your original key for recovery.",
                          file=sys.stderr, flush=True)
            raise
        if saved_license is not None:
            with stage("Restore existing MODELLER configuration (protected; output not logged)"):
                configure_modeller(manager, prefix, saved_license, verify_import=saved_license_ready)
            saved_license = None
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
            with stage("Configure MODELLER license (protected; output not logged)"):
                configure_modeller(manager, prefix, key)
        elif interactive and not args.skip_modeller and not args.update and not preserved_license and not args.gui_only:
            with waiting_for_input():
                key = getpass.getpass("MODELLER key (empty to configure later): ")
            if key:
                with stage("Configure MODELLER license (protected; output not logged)"):
                    configure_modeller(manager, prefix, key)
        if args.profile == "desktop" and not args.configure_only:
            # Isolate GUI artifacts from the checkout and other installations.
            gui = prefix / "share/dvbfixer/gui"
            gui.parent.mkdir(parents=True, exist_ok=True)
            with stage("Prepare Desktop assets (preserving workspaces and local settings)"):
                shutil.copytree(root / "gui", gui, dirs_exist_ok=True,
                                ignore=shutil.ignore_patterns("node_modules", "dist", "dist-server", "structures", ".env*"))
            run("Install GUI dependencies", manager.run(prefix, ["npm", "ci"]), cwd=gui)
            run("Build Desktop GUI", manager.run(prefix, ["npm", "run", "build"]), cwd=gui)
        if not args.gui_only:
            run("Check scientific linear algebra", manager.run(prefix, ["python", str(root / "scripts/check_scientific_linalg.py")]))
        if not args.gui_only:
            with stage("Check MODELLER license (protected; output not logged)"):
                license_probe = subprocess.run(manager.run(prefix, ["python", "-c", "import modeller"]),
                                               env=child_environment(), capture_output=True, timeout=30, check=False)
        if not args.gui_only and (not args.configure_only or license_probe.returncode == 0):
            evidence = prefix / "share/dvbfixer/installation-smoke" / uuid.uuid4().hex
            smoke = ["python", str(root / "scripts/smoke_installation_profile.py"),
                     "--output", str(evidence), "--root", str(root)]
            if license_probe.returncode == 0:
                smoke.append("--with-modeller")
            if args.configure_only:
                smoke.extend(["--only", "model", "--only", "homology"])
            run("Run functional science checks", manager.run(prefix, smoke))
            print(f"Functional smoke evidence: {evidence}", flush=True)
        doctor_options = [] if args.gui_only else ["--profile", args.profile, *[arg for addon in args.addon for arg in ("--addon", addon)]]
        result = logged_run(manager.run(prefix, ["dvbfixer", "doctor", *doctor_options]),
                            env=child_environment(), log=logs / "Verify installation readiness.log",
                            allowed_codes=(0, 1, 2))
        install_result.readiness(result.returncode, logs / "Verify installation readiness.log")
        if args.profile == "desktop" and result.returncode in (0, 2):
            run("Check Desktop launcher", manager.run(prefix, ["python", "-m", "dvbfixer.gui_launcher", "--check"]))
            if args.desktop_shortcut:
                from dvbfixer.desktop_shortcuts import install_desktop_shortcut
                shortcut = install_desktop_shortcut(prefix)
                print(f"Desktop application shortcut: {shortcut}", flush=True)
            print("Start the GUI: dvbfixer gui\nStop the GUI: dvbfixer gui --stop", flush=True)
        if extension is not None and result.returncode in (0, 2):
            extension.commit()
        elif extension is not None:
            print("Extension did not complete; profile ownership was not changed. Inspect the logs and repeat --extend.", flush=True)
        if result.returncode == 0:
            label = "GUI only (scientific stack unchanged)" if args.gui_only else args.profile
            print(f"\nREADY: DVBFixer {label}. Installation checks passed.", flush=True)
        elif result.returncode == 2:
            print("\nACTION REQUIRED: installation is not fully ready; follow the Doctor findings above.", flush=True)
        else:
            print("\nFAILED: installation readiness checks failed; inspect the findings above.", flush=True)
        activation = shlex.quote(args.activation_name or str(prefix))
        if args.gui_only:
            print("GUI-only operation: no Standard/Desktop scientific-profile qualification is claimed.", flush=True)
        print(f"Logs: {logs}\nActivate: {args.manager} activate {activation}", flush=True)
        return result.returncode
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        install_result.failure(exc)
        print(f"error: {exc}", file=sys.stderr)
        if extension is not None:
            print("Extension did not complete. Native changes may already be present; inspect the retained logs "
                  "and repeat --extend to review a fresh plan. Profile ownership was not published.", file=sys.stderr)
        return 1
    finally:
        if extension is not None:
            extension.close()


def main(argv: list[str] | None = None) -> int:
    from dvbfixer import install_result

    result = install_result.InstallResult()
    token = install_result.current.set(result)
    code = 1
    report = True
    try:
        code = _main(argv)
        return code
    except (KeyboardInterrupt, EOFError):
        code = 130
        install_result.cancelled()
        print("\nCancelled. Any existing diagnostics are retained. Rerun the same operation; "
              "--resume is for interrupted helper-managed installations.",
              file=sys.stderr, flush=True)
        return 130
    except SystemExit as exc:
        # argparse help is not an installation; usage errors keep argparse's status.
        code = 0 if exc.code == 0 else 1
        report = code != 0
        if report:
            install_result.failure("Invalid installer arguments; correct the usage error above and consult --help.")
        raise
    except Exception as exc:
        install_result.failure(exc)
        raise
    finally:
        try:
            if report:
                install_result.finish(result, code)
        finally:
            install_result.current.reset(token)


if __name__ == "__main__":
    raise SystemExit(main())
