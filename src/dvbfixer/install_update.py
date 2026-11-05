"""Explicit package-only updates for pre-profile and manually maintained installs."""

from __future__ import annotations

import argparse
import shlex
import shutil
import uuid
from pathlib import Path

from dvbfixer import install_result
from dvbfixer.install import Manager, child_environment, recorded_installation, run_checked
from dvbfixer.install_discovery import MARKER_SUFFIX, package_details, runtime_component
from dvbfixer.install_native_tools import preserve_native_openbabel_tool
from dvbfixer.install_progress import logged_run
from dvbfixer.install_ui import confirm


def update_package_only(args: argparse.Namespace, *, interactive: bool) -> int:
    """Update just the selected package, without adopting or solving its environment."""
    if args.prefix is None or args.manager not in {"micromamba", "mamba", "conda"}:
        raise ValueError("Select an installation and environment manager first")
    if any((args.runtime, args.configure_only, args.addon, args.desktop_shortcut, args.modeller_key_fd is not None)):
        raise ValueError("Package-only updates cannot change runtimes, profiles/add-ons, shortcuts or licenses")
    prefix = args.prefix.expanduser().absolute()
    if prefix.is_symlink():
        raise ValueError("The target prefix must not be a symlink")
    prefix = prefix.resolve()
    if not (prefix / "conda-meta/history").is_file():
        raise ValueError("Selected installation is not an existing conda-family environment")
    if runtime_component(prefix):
        raise ValueError("A frozen backend interpreter is not a package-update target")
    marker = prefix.with_name(prefix.name + MARKER_SUFFIX)
    if marker.exists() or marker.is_symlink():
        record = recorded_installation(prefix)
        if record["profile"] not in {"standard", "desktop"}:
            raise ValueError("A frozen backend is not a package-update target")
    version, editable, found = package_details(prefix)
    if not found:
        raise ValueError("No existing DVBFixer package or CLI entry point found; refusing to adopt this environment")
    executable = shutil.which(args.manager)
    if executable is None:
        raise ValueError(f"Install the selected manager first: {args.manager}")
    manager = Manager(args.manager, executable)
    manager.validate()
    repository = Path(__file__).resolve().parents[2]
    if not (repository / "pyproject.toml").is_file():
        raise ValueError("Run the updater from the intended DVBFixer source checkout")
    from dvbfixer import __version__

    print(f"\nPackage-only update: {prefix.name}\nVersion: {version} -> {__version__}\n"
          f"Source checkout: {repository}\n"
          "Only DVBFixer will be reinstalled. Scientific dependencies, MODELLER configuration,\n"
          "workspace data and unrelated packages are not update targets.\n"
          "This does not adopt the environment as a managed Standard/Desktop profile.\n"
          "Close this installation's GUI and jobs first; package updates have no automatic rollback.", flush=True)
    if editable:
        print("Editable mode will be retained, pointing to the source checkout above.", flush=True)
    print("GUI assets will not be rebuilt by this package-only operation.", flush=True)
    if args.dry_run:
        print("Dry run only; no installation was changed.", flush=True)
        return 0
    if interactive and not confirm("Apply the package-only update plan? Dependencies and GUI assets will not change."):
        install_result.cancelled()
        print("Cancelled; nothing changed.", flush=True)
        return 0
    logs = prefix.parent / f"{prefix.name}.installation-logs" / uuid.uuid4().hex
    logs.mkdir(mode=0o700, parents=True, exist_ok=False)
    print(f"Update logs: {logs}", flush=True)
    install_result.log_directory(logs)
    install_result.started()
    mode = ["-e"] if editable else []
    command = manager.run(prefix, ["python", "-m", "pip", "install", "--upgrade", "--no-deps", *mode, str(repository)])
    with preserve_native_openbabel_tool(prefix):
        run_checked(command, log=logs / "Update DVBFixer package only.log")
    result = logged_run(manager.run(prefix, ["python", "-m", "dvbfixer.cli", "doctor"]),
                        env=child_environment(), log=logs / "Check existing environment.log",
                        allowed_codes=(0, 1, 2))
    install_result.readiness(result.returncode, logs / "Check existing environment.log")
    if result.returncode == 0:
        print("\nDVBFixer package updated. Existing-environment diagnostic checks completed.", flush=True)
    else:
        print("\nDVBFixer package updated, but environment checks need attention; inspect the findings above.", flush=True)
    print("Scientific dependencies were not upgraded; no Standard/Desktop profile readiness is claimed.", flush=True)
    activation = shlex.quote(args.activation_name or str(prefix))
    print(f"Logs: {logs}\nActivate: {args.manager} activate {activation}", flush=True)
    return result.returncode
