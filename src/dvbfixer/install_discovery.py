"""Read-only manager discovery and name/number-based installer selection."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from importlib import metadata as package_metadata
from pathlib import Path

from dvbfixer.install import Manager, child_environment, recorded_installation
from dvbfixer.install_progress import stage
from dvbfixer.install_ui import choose
from dvbfixer.installation_profiles import load_catalog

MARKER_SUFFIX = ".dvbfixer-install.json"


@dataclass(frozen=True)
class Installation:
    prefix: Path
    profile: str
    version: str
    managed: bool = True
    editable: bool = False
    package_present: bool = True

    @property
    def backend(self) -> bool:
        return self.profile.startswith("runtime:")

    @property
    def label(self) -> str:
        if not self.managed:
            if not self.package_present:
                return f"Existing environment — {self.prefix.name} (DVBFixer not installed)"
            return f"DVBFixer {self.version} — {self.prefix.name} (existing/manual)"
        if self.backend:
            return f"Additional backend {self.profile.removeprefix('runtime:')} — {self.prefix.name}"
        return f"DVBFixer {self.profile.title()} {self.version} — {self.prefix.name}"


def package_details(prefix: Path) -> tuple[str, bool, bool]:
    """Find wheel, legacy egg-info/egg-link and metadata-incomplete CLI installs.

    Reading distribution metadata does not execute .pth files, linked sources or
    the discovered environment's Python.
    """
    versions = set()
    editable = False
    found = False
    for site in prefix.glob("lib/python*/site-packages"):
        if any(site.glob("__editable__.dvbfixer*.pth")):
            found, editable = True, True
        locations = [site]
        for link in site.glob("*.egg-link"):
            if link.stem.lower().replace("_", "-") != "dvbfixer":
                continue
            lines = link.read_text(encoding="utf-8").splitlines()
            if lines and lines[0].strip():
                linked = Path(lines[0].strip())
                locations.append(linked if linked.is_absolute() else site / linked)
                found, editable = True, True
        for distribution in package_metadata.distributions(path=[str(path) for path in locations]):
            name = distribution.metadata.get("Name") or ""
            if re.sub(r"[-_.]+", "-", name).lower() != "dvbfixer":
                continue
            found = True
            if distribution.version:
                versions.add(distribution.version)
            direct = distribution.read_text("direct_url.json")
            if direct:
                try:
                    payload = json.loads(direct)
                except ValueError:
                    print(f"WARNING: Editable-install metadata is unreadable in {prefix.name}; "
                          "the existing package is still listed.", flush=True)
                    continue
                if isinstance(payload, dict) and isinstance(payload.get("dir_info"), dict):
                    editable |= payload["dir_info"].get("editable") is True
    found |= (prefix / "bin/dvbfixer").is_file()
    version = next(iter(versions)) if len(versions) == 1 else "(version unknown)"
    return version, editable, found


def installed_version(prefix: Path) -> str:
    return package_details(prefix)[0]


def runtime_component(prefix: Path) -> bool:
    """A frozen backend or its child interpreter is not a manual CLI target."""
    if (prefix / "runtime-plan.json").exists() or (prefix.parent / "runtime-plan.json").exists():
        return True
    parent_marker = prefix.parent.with_name(prefix.parent.name + MARKER_SUFFIX)
    if parent_marker.exists() or parent_marker.is_symlink():
        try:
            return recorded_installation(prefix.parent)["profile"].startswith("runtime:")
        except (OSError, ValueError):
            return True
    return False


def manager_inventory(manager: Manager) -> tuple[list[Path], list[Path]]:
    def query(*args: str) -> dict:
        result = subprocess.run([manager.executable, *args, "--json"], env=child_environment(),
                                capture_output=True, text=True, check=True, timeout=30)
        payload = json.loads(result.stdout)
        if not isinstance(payload, dict):
            raise ValueError("Environment manager returned an invalid inventory")
        return payload

    with stage(f"Find installations using {manager.name}"):
        info = query("info")
        listing = query("env", "list")
    directories = info.get("envs_dirs", info.get("envs directories"))
    if directories is None:
        root = info.get("root_prefix", info.get("base environment"))
        directories = [str(Path(root) / "envs")] if isinstance(root, str) else []
    envs = listing.get("envs")
    if (not isinstance(directories, list) or not directories
            or not isinstance(envs, list)
            or not all(isinstance(value, str) and Path(value).expanduser().is_absolute()
                       for value in [*directories, *envs])):
        raise ValueError("Environment manager did not provide valid environment directories; use advanced --prefix")
    return ([Path(value).expanduser() for value in directories], [Path(value).expanduser() for value in envs])


def manager_base(manager: Manager) -> Path:
    result = subprocess.run([manager.executable, "info", "--json"], env=child_environment(),
                            capture_output=True, text=True, check=True, timeout=30)
    info = json.loads(result.stdout)
    root = info.get("root_prefix", info.get("base environment")) if isinstance(info, dict) else None
    if not isinstance(root, str) or not Path(root).is_absolute():
        raise ValueError("Cannot identify the manager base environment safely; no changes made")
    return Path(root).resolve()


def protected_environment(prefix: Path, base: Path | None = None, envs: list[Path] | None = None) -> bool:
    """Keep application bundles and all discoverable manager bases out of extension targets."""
    prefix = prefix.expanduser().resolve()
    bases = {env.parent.parent.resolve() for env in envs or [] if env.parent.name == "envs"}
    if base is not None:
        bases.add(base.resolve())
    if prefix in bases:
        return True
    if (prefix / "condabin").is_dir() and (prefix / "bin/conda").is_file():
        return True
    return any(parent.suffix == ".app" and (parent / "Contents" == prefix or parent / "Contents" in prefix.parents)
               for parent in prefix.parents)


def discover_installations(directories: list[Path], envs: list[Path], *, repository: Path,
                           report_manual: bool = False, include_manual: bool = False,
                           include_empty: bool = False) -> list[Installation]:
    # No recursive filesystem scan. Composite runtimes own their sampling/
    # refinement/ environments, so inspect their immediate parent as well.
    roots = {*directories, Path.cwd(), repository}
    candidates = set(envs) | {path.parent for path in envs}
    for root in roots:
        for marker in root.glob(f"*{MARKER_SUFFIX}"):
            candidates.add(marker.with_name(marker.name.removesuffix(MARKER_SUFFIX)))
    installations = {}
    for candidate in sorted(candidates):
        marker = candidate.with_name(candidate.name + MARKER_SUFFIX)
        if not marker.exists() and not marker.is_symlink():
            if (report_manual or include_manual or include_empty) and candidate in envs and not candidate.is_symlink():
                try:
                    if runtime_component(candidate):
                        continue
                    version, editable, found = package_details(candidate)
                    if (found and include_manual) or (include_empty and (candidate / "conda-meta/history").is_file()):
                        prefix = candidate.resolve()
                        installations[prefix] = Installation(prefix, "manual", version, managed=False, editable=editable, package_present=found)
                    elif found:
                        print(f"Existing/manual DVBFixer: {candidate.name} — {version}; select update for a package-only update.", flush=True)
                except (OSError, ValueError) as exc:
                    print(f"WARNING: Could not read environment metadata for {candidate.name}: {exc}", flush=True)
            continue
        try:
            record = recorded_installation(candidate)
            profile = record["profile"]
            if profile not in {"standard", "desktop"} and not profile.startswith("runtime:"):
                raise ValueError("Unrecognized installation kind")
            prefix = candidate.resolve()
            version, editable, found = ("", False, True) if profile.startswith("runtime:") else package_details(prefix)
            installations[prefix] = Installation(prefix, profile, version, editable=editable, package_present=found)
        except (OSError, ValueError) as exc:
            print(f"WARNING: Ignoring unusable installation record {marker}: {exc}", flush=True)
    return sorted(installations.values(), key=lambda entry: (entry.profile, entry.prefix.name, str(entry.prefix)))


def runtime_status(prefix: Path) -> str:
    """Describe retained provisioning evidence, not scientific qualification."""
    evidence = prefix / "evidence/installation.json"
    try:
        if evidence.stat().st_size > 1024 * 1024:
            return "status unavailable (oversized evidence)"
        record = json.loads(evidence.read_text())
        state = record.get("state") if isinstance(record, dict) else None
        states = {
            "ready": "provisioned; revalidation required on resume",
            "action-required": "action required (assets/setup incomplete)",
            "error": "previous provisioning failed",
            "planned": "provisioning started; not ready",
        }
        return states.get(state, "status not recorded; not verified ready") if isinstance(state, str) else "status not recorded; not verified ready"
    except FileNotFoundError:
        return "provisioning record only; readiness not verified"
    except (OSError, ValueError):
        return "status unavailable (unreadable evidence)"


def select_existing(entries: list[Installation], *, other_entries: list[Installation] | None = None,
                    backend: bool = False) -> Installation:
    if not entries and other_entries:
        print("No existing DVBFixer installation found; choose an environment to add it to.", flush=True)
        return select_existing(other_entries)
    if not entries:
        raise ValueError("No matching DVBFixer installations found. No environment was changed. "
                         "Check the selected manager; custom locations can use advanced --prefix.")
    title = "Choose a separate backend runtime to continue" if backend else "Choose an existing installation"
    print(f"\n{title} (0 cancels):", flush=True)
    duplicates = {entry.prefix.name for entry in entries if sum(other.prefix.name == entry.prefix.name for other in entries) > 1}
    for number, entry in enumerate(entries, 1):
        print(f"  {number}. {entry.label}", flush=True)
        if backend or entry.prefix.name in duplicates or not entry.package_present:
            print(f"     Location: {entry.prefix}", flush=True)
        if backend:
            print(f"     State: {runtime_status(entry.prefix)}", flush=True)
    if other_entries:
        print(f"  {len(entries) + 1}. Other existing environments…", flush=True)
    while True:
        choice = input("Runtime number: " if backend else "Installation number: ").strip()
        if choice == "0":
            raise EOFError()
        if choice.isdecimal() and len(choice) < 8 and 1 <= int(choice) <= len(entries):
            return entries[int(choice) - 1]
        if other_entries and choice == str(len(entries) + 1):
            return select_existing(other_entries)
        print("Choose a number from the list, or 0 to cancel.", flush=True)


def validate_name(name: str) -> None:
    if (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name)
            or name in {"base", "root"}):
        raise ValueError("Use an environment name (letters, numbers, '.', '_' or '-'), not a path; base/root are reserved")


def creation_directory(directories: list[Path]) -> Path:
    """Use manager directory preference order without creating anything on discovery."""
    for directory in directories:
        parent = directory
        while not parent.exists() and parent != parent.parent:
            parent = parent.parent
        if parent.is_dir() and os.access(parent, os.W_OK | os.X_OK):
            return directory
    raise ValueError("The manager has no writable environment directory; configure its envs_dirs or use advanced --prefix")


def new_named_prefix(directory: Path, envs: list[Path], *, default: str,
                     name: str | None, interactive: bool) -> Path:
    def occupied(value: str) -> bool:
        path = directory / value
        return (path.exists() or path.is_symlink()
                or path.with_name(path.name + MARKER_SUFFIX).exists()
                or path.with_name(path.name + MARKER_SUFFIX).is_symlink()
                or any(env.name == value for env in envs))

    suggested = default
    suffix = 2
    while occupied(suggested):
        suggested = f"{default}-{suffix}"
        suffix += 1
    while True:
        value = name if name is not None else input(f"Environment name [{suggested}]: ").strip() or suggested
        try:
            validate_name(value)
            if occupied(value):
                raise ValueError(f"Environment name '{value}' is already in use; choose a different name")
            return directory / value
        except ValueError as exc:
            if name is not None or not interactive:
                raise
            print(str(exc), flush=True)


def select_diffusion_for_application(args: argparse.Namespace, *, runtimes: tuple[str, ...], repository: Path) -> None:
    """Expose diffusion here without co-installing Torch into the chosen app prefix."""
    from dvbfixer.install import platform_subdir

    if args.addon or args.modeller_key_fd is not None or args.desktop_shortcut:
        raise ValueError("Diffusion provisioning is separate from scientific add-ons, license changes and GUI shortcuts; "
                         "configure those in a separate operation")
    catalog = load_catalog()["external"]
    labels = {
        "protenix-cuda": "Protenix v1 — NVIDIA CUDA",
        "protpardelle-cuda": "Protpardelle cc94 — NVIDIA CUDA",
        "protpardelle-mps": "Protpardelle cc94 — Apple MPS (experimental/case-tested)",
    }
    available = [name for name in runtimes if catalog[name]["kind"] == "separate-runtime"
                 and platform_subdir() in catalog[name]["platforms"]]
    if not available:
        raise ValueError("No maintained diffusion runtime is available for this platform")
    application = args.prefix.expanduser().resolve()
    print(f"\nScientific installation remains unchanged: {application.name}\n"
          "Diffusion uses a separate runtime environment, not PyTorch inside Standard/Desktop.\n"
          "Validated launchers use existing per-user discovery and may be visible to other compatible DVBFixer installations.", flush=True)
    runtime = choose("Diffusion backend (separate environment)", [(name, labels.get(name, name)) for name in available])
    print(f"Using the selected environment manager: {args.manager}", flush=True)
    while True:
        operation = choose("Diffusion runtime installation", [
            ("install", "Create a new isolated runtime"),
            ("resume", "Continue an existing helper-recorded runtime"),
            ("manager", f"Change runtime manager (currently {args.manager})"),
        ])
        if operation != "manager":
            break
        args.manager = choose("Runtime environment manager", [(name, name) for name in ("micromamba", "mamba", "conda")], default=args.manager)
    args.application_prefix = application
    args.extend, args.update, args.gui_only, args.package_only = False, False, False, False
    args.resume = operation == "resume"
    args.runtime, args.profile, args.prefix, args.name, args.named_parent = runtime, None, None, None, None
    resolve_selection(args, interactive=True, runtimes=runtimes, repository=repository)
    if args.prefix.resolve() == application:
        raise ValueError("The diffusion runtime must not reuse the scientific installation directory")


def resolve_selection(args: argparse.Namespace, *, interactive: bool,
                      runtimes: tuple[str, ...], repository: Path) -> None:
    """Fill identity/prefix from existing records or a manager-owned environment name."""
    args.activation_name = None
    args.full_update_selected = False
    package_only_requested = args.package_only
    selected = None
    kind = "application"
    if interactive and not any((args.profile, args.runtime, args.resume, args.update, args.extend)):
        choice = choose("DVBFixer setup", [
            ("install", "New Standard/Desktop installation"),
            ("update", "Update an existing DVBFixer installation"),
            ("extend", "Add a profile or GUI to an existing environment"),
            ("backend", "Configure an additional model/QM backend"),
        ])
        if choice == "update":
            args.update = True
        elif choice == "extend":
            args.extend = True
        elif choice == "backend":
            args.runtime = choose("Additional backend", [(name, name) for name in runtimes])
    if (interactive and args.resume and not any((args.prefix, args.name, args.profile, args.runtime, args.configure_only))):
        kind = choose("Resume an existing installation", [("application", "DVBFixer application"), ("backend", "Additional backend")])
    if args.runtime:
        kind = "backend"
        if args.runtime not in runtimes:
            raise ValueError("Unknown additional backend")
    if interactive and not any((args.profile, args.runtime, args.update, args.resume, args.extend)):
        args.profile = choose("Application profile", [("standard", "Standard — scientific CLI"), ("desktop", "Desktop — Standard plus GUI")])
    if args.profile is not None and args.profile not in {"standard", "desktop"}:
        raise ValueError("Choose the Standard or Desktop application profile")
    if interactive:
        args.manager = args.manager or choose("Environment manager", [(name, name) for name in ("micromamba", "mamba", "conda")])
    if args.prefix is None and (interactive or args.name is not None):
        if args.manager not in {"micromamba", "mamba", "conda"}:
            raise ValueError("Choose micromamba, mamba or conda")
        executable = shutil.which(args.manager)
        if executable is None:
            raise ValueError(f"Install the selected manager first: {args.manager}")
        manager = Manager(args.manager, executable)
        manager.validate()
        directories, envs = manager_inventory(manager)
        if args.update or args.resume or args.extend:
            with stage("Inspect existing environments (package metadata only)"):
                entries = discover_installations(directories, envs, repository=repository,
                                                 include_manual=args.update or args.extend, include_empty=args.extend)
            if args.extend:
                base = manager_base(manager)
                entries = [entry for entry in entries if not protected_environment(entry.prefix, base, envs)
                           and (entry.prefix / "conda-meta/history").is_file()]
            entries = [entry for entry in entries if (not (args.update or args.extend) or not entry.backend)
                       and (args.extend or not args.profile or entry.profile == args.profile)
                       and (not args.runtime or entry.profile == f"runtime:{args.runtime}")]
            if args.extend and args.gui_only:
                entries = [entry for entry in entries if entry.package_present]
            if args.name is not None:
                entries = [entry for entry in entries if entry.prefix.name == args.name]
                if len(entries) != 1:
                    raise ValueError("Environment name is missing or ambiguous; use interactive selection or advanced --prefix")
                selected = entries[0]
            else:
                entries = [entry for entry in entries if entry.backend == (kind == "backend")]
                if args.extend:
                    primary = [entry for entry in entries if entry.package_present]
                    others = [entry for entry in entries if not entry.package_present]
                    selected = select_existing(primary, other_entries=others)
                else:
                    selected = select_existing(entries, backend=kind == "backend")
            args.prefix = selected.prefix
            if not selected.managed and not args.extend:
                if not interactive and not args.package_only:
                    raise ValueError("This is an existing/manual installation. Add --package-only to explicitly update "
                                     "DVBFixer without changing its scientific dependencies or adopting the environment.")
                args.package_only = True
        else:
            from dvbfixer import __version__

            identity = args.runtime or args.profile
            if identity is None:
                raise ValueError("Choose a profile or additional backend")
            if not interactive and args.name is None:
                raise ValueError("Supply --name or advanced --prefix")
            args.prefix = new_named_prefix(creation_directory(directories), envs,
                default=f"dvbfixer-{identity}-{__version__}", name=args.name, interactive=interactive)
            args.named_parent = args.prefix.parent
        if args.prefix.parent.resolve() in {directory.resolve() for directory in directories}:
            args.activation_name = args.prefix.name
    if args.prefix is not None:
        label = "Selected separate runtime" if kind == "backend" else "Selected installation"
        print(f"{label}: {args.prefix.expanduser().absolute()}", flush=True)
    if interactive and args.update and not package_only_requested:
        candidate = args.prefix.expanduser().absolute()
        marker = candidate.with_name(candidate.name + MARKER_SUFFIX)
        if not marker.exists() and not marker.is_symlink() and package_details(candidate)[2] and not runtime_component(candidate):
            update_scope = choose("What should be updated?", [
                ("package", "Only DVBFixer — keep dependencies and GUI assets unchanged"),
                ("all", "Full installation — DVBFixer, profile dependencies and GUI when Desktop is selected"),
            ])
            if update_scope == "all":
                args.update, args.extend, args.package_only = False, True, False
                args.full_update_selected = True
            else:
                args.package_only = True
    if interactive and args.extend and not args.profile and not args.gui_only:
        title = "Choose the full update profile" if args.full_update_selected else "Add to the selected environment"
        component = choose(title, [
            ("standard", "Standard — complete the scientific CLI profile"),
            ("desktop", "Desktop — Standard plus built GUI"),
            ("gui", "GUI only — keep the current scientific stack"),
            ("diffusion", "Diffusion backend — provision or resume a separate model runtime"),
        ])
        if component == "diffusion":
            select_diffusion_for_application(args, runtimes=runtimes, repository=repository)
            return
        args.gui_only = component == "gui"
        args.profile = "desktop" if args.gui_only else component
    if args.extend and args.gui_only:
        args.profile = "desktop"
    if interactive and (args.resume or args.update) and args.prefix is not None:
        candidate = args.prefix.expanduser().absolute()
        marker = candidate.with_name(candidate.name + MARKER_SUFFIX)
        if args.update and not marker.exists() and not marker.is_symlink():
            _, _, found = package_details(candidate)
            if found and not runtime_component(candidate):
                args.package_only = True
        if args.update and args.package_only:
            return
        record = recorded_installation(args.prefix)
        if record["profile"].startswith("runtime:"):
            if args.update:
                raise ValueError("--update is for DVBFixer; use backend --resume for frozen runtimes")
            runtime = record["profile"].removeprefix("runtime:")
            if args.profile or (args.runtime and args.runtime != runtime):
                raise ValueError("Selected backend differs from the existing installation")
            args.runtime = runtime
        else:
            if args.runtime or (args.profile and args.profile != record["profile"]):
                raise ValueError("Selected profile differs from the existing installation")
            args.profile = record["profile"]
