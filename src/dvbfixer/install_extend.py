"""Preview and explicitly attach profile/GUI components to existing environments.

Native and pip solvers remain separate. Both are previewed before confirmation;
native plan/environment drift is rejected, and pip is rechecked after the native
transaction. This is not a transactional rollback mechanism.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import uuid
from importlib import metadata as package_metadata
from pathlib import Path

from dvbfixer.install import Manager, child_environment, recorded_installation, validate_prefix_path
from dvbfixer.install_discovery import (
    MARKER_SUFFIX,
    manager_base,
    package_details,
    protected_environment,
    runtime_component,
)
from dvbfixer.install_progress import stage
from dvbfixer.installation_profiles import composition, load_catalog

NATIVE_PYPI_NAMES = {"openff-toolkit-base": "openff-toolkit"}


class ProfileStackConflict(ValueError):
    """A full scientific profile conflicts with a stack GUI-only can leave intact."""


def native_records(prefix: Path) -> dict[str, dict]:
    records = {}
    for path in (prefix / "conda-meta").glob("*.json"):
        record = json.loads(path.read_text())
        if not isinstance(record, dict) or not isinstance(record.get("name"), str) or not isinstance(record.get("version"), str):
            raise ValueError(f"Invalid native package metadata: {path.name}")
        records[record["name"]] = record
    return records


def fingerprint(prefix: Path) -> str:
    digest = hashlib.sha256()
    paths = [*(prefix / "conda-meta").glob("*"), *prefix.glob("lib/python*/site-packages/*.dist-info/METADATA"),
             *prefix.glob("lib/python*/site-packages/*.dist-info/direct_url.json"),
             *prefix.glob("lib/python*/site-packages/*.egg-link"),
             *prefix.glob("lib/python*/site-packages/*.egg-info/PKG-INFO"),
             *prefix.glob("lib/python*/site-packages/*.egg-info")]
    for path in sorted(paths):
        if path.is_file():
            digest.update(str(path.relative_to(prefix)).encode())
            digest.update(path.read_bytes())
    return digest.hexdigest()


def extension_composition(args: argparse.Namespace) -> tuple[str, list[str]]:
    prefix = args.prefix.expanduser().absolute()
    if prefix.is_symlink():
        raise ValueError("The existing environment must not be a symlink")
    prefix = prefix.resolve()
    if not (prefix / "conda-meta/history").is_file():
        raise ValueError("Select an existing conda-family environment")
    if runtime_component(prefix):
        raise ValueError("Frozen backend roots and interpreters cannot be converted to application profiles")
    profile, addons = args.profile, set(args.addon)
    marker = prefix.with_name(prefix.name + MARKER_SUFFIX)
    if marker.exists() or marker.is_symlink():
        record = recorded_installation(prefix)
        if record["profile"] not in {"standard", "desktop"}:
            raise ValueError("A frozen backend cannot be converted to an application profile")
        if not args.gui_only:
            # Adding Standard must not erase Desktop ownership or existing add-ons.
            profile = "desktop" if record["profile"] == "desktop" else profile
            addons.update(record["addons"])
    if args.gui_only and (args.addon or args.modeller_key_fd is not None):
        raise ValueError("GUI-only installation cannot change scientific add-ons or the MODELLER license")
    return profile, sorted(addons)


def native_request(manager: Manager, prefix: Path, specs: list[str]) -> list[str]:
    command = [manager.executable, "install", "--yes", "--prefix", str(prefix), "--override-channels"]
    for channel in load_catalog()["channels"]:
        command.extend(["--channel", channel])
    return [*command, *specs]


def solve_native(command: list[str]) -> dict:
    with stage("Preview native dependency changes (no packages installed)"):
        result = subprocess.run([*command, "--dry-run", "--json"], env=child_environment(),
                                text=True, capture_output=True, check=False)
    try:
        payload = json.loads(result.stdout)
    except ValueError:
        raise ValueError("The native solver did not return a valid plan; no packages were changed")
    if result.returncode or not isinstance(payload, dict) or payload.get("success") is not True:
        message = payload.get("message", payload.get("error", "dependency conflict")) if isinstance(payload, dict) else "invalid plan"
        raise ValueError(f"Native dependency plan failed: {message}. Use a new environment if this stack is incompatible.")
    actions = payload.get("actions", {})
    if not isinstance(actions, dict):
        raise ValueError("Malformed native transaction plan")
    normalized = {}
    for action in ("UNLINK", "LINK"):
        entries = actions.get(action, [])
        if not isinstance(entries, list) or any(not isinstance(entry, dict)
                or not isinstance(entry.get("name"), str) or not entry["name"]
                or not isinstance(entry.get("version"), str) or not entry["version"]
                or not isinstance(entry.get("build_string", entry.get("build", "")), str) for entry in entries):
            raise ValueError("Incomplete native transaction metadata; refusing to apply")
        normalized[action] = sorted(
            [{"name": entry["name"], "version": entry["version"],
              "build": entry.get("build_string", entry.get("build", "")),
              "python_package": any(dependency.split()[0] in {"python", "python_abi"}
                                    for dependency in entry.get("depends", []) if isinstance(dependency, str) and dependency.split())}
             for entry in entries],
            key=lambda entry: (entry["name"], entry["version"], entry["build"]),
        )
    return normalized


def pip_installations(report: Path) -> set[tuple[str, str]]:
    payload = json.loads(report.read_text())
    if not isinstance(payload, dict) or not isinstance(payload.get("install"), list):
        raise ValueError("Pip did not return an installation report")
    installs = set()
    for item in payload["install"]:
        metadata = item.get("metadata", {}) if isinstance(item, dict) else {}
        if not isinstance(metadata.get("name"), str) or not isinstance(metadata.get("version"), str):
            raise ValueError("Invalid pip installation metadata")
        installs.add((metadata["name"].lower().replace("_", "-"), metadata["version"]))
    return installs


class ExtensionPlan:
    def __init__(self, args: argparse.Namespace, manager: Manager):
        self.args, self.manager = args, manager
        self.prefix = args.prefix.expanduser().resolve()
        self.repository = Path(__file__).resolve().parents[2]
        self.temporary = tempfile.TemporaryDirectory(prefix="dvbfixer-extension-")
        self.directory = Path(self.temporary.name)
        self.logs: Path | None = None
        self.marker = self.prefix.with_name(self.prefix.name + MARKER_SUFFIX)
        self.original_marker = self.marker.read_bytes() if self.marker.exists() else None
        self.before = fingerprint(self.prefix)

    def close(self) -> None:
        if self.logs is None and any(self.directory.glob("*.log")):
            retained = self.repository / "tmp" / f"installer-profile-preview-{uuid.uuid4().hex}"
            retained.mkdir(mode=0o700, parents=True)
            for log in self.directory.glob("*.log"):
                shutil.copy2(log, retained / log.name)
            from dvbfixer import install_result

            install_result.log_directory(retained)
            result = install_result.current.get()
            if result is not None and result.stage_log is not None and result.stage_log.parent == self.directory:
                result.stage_log = retained / result.stage_log.name
            print(f"Preview diagnostics retained: {retained}", flush=True)
        self.temporary.cleanup()

    def preview(self) -> None:
        args = self.args
        validate_prefix_path(self.prefix)
        # Never attach a profile to a manager's base environment by accident.
        if protected_environment(self.prefix, manager_base(self.manager)):
            raise ValueError("Do not attach application profiles to the manager base or an application bundle; select a named environment")
        self.editable = False
        if args.gui_only:
            _, self.editable, found = package_details(self.prefix)
            if not found:
                raise ValueError("GUI-only needs an existing DVBFixer installation; select the full Desktop profile instead")
        records = native_records(self.prefix)
        sites = [str(path) for path in self.prefix.glob("lib/python*/site-packages")]
        installed_names = {(d.metadata.get("Name") or "").lower().replace("_", "-")
                           for d in package_metadata.distributions(path=sites)}
        if not args.gui_only and ({"pytorch", "pytorch-gpu", "torch"} & (set(records) | installed_names)):
            raise ProfileStackConflict("Keep PyTorch/diffusion runtimes separate from Standard/Desktop; "
                                       "select GUI-only or create a new scientific profile environment")
        python = records.get("python", {}).get("version", "")
        if not python or len(python.split(".")) < 2:
            raise ValueError("Existing environment has no identifiable Python; create a new profile environment")
        packages = composition(args.profile, tuple(args.addon))
        native = [p for p in packages if p.get("source") != "pip" and (not args.gui_only or p["group"] == "gui")]
        if not args.gui_only and "openbabel" in records:
            raise ProfileStackConflict("The complete Standard/Desktop profile uses openbabel-wheel, but this environment "
                                       "has conda Open Babel. They must not be overlaid. For GUI only replace --profile with --gui-only; "
                                       "for the complete profile choose a new environment.")
        specs = [p["channel"] + "::" + p["spec"] for p in native]
        # Keep the Python ABI stable so pip can be previewed in this interpreter.
        # Unsupported Python/profile combinations fail in the solver, not midway
        # through replacing the environment's interpreter.
        if not args.gui_only:
            specs.append(f"conda-forge::python={'.'.join(python.split('.')[:2])}.*")
        if args.gui_only:
            # Ask the solver for a Node version compatible with the existing
            # stack, rather than letting it eagerly update Python/shared libs
            # and then rejecting an otherwise avoidable change.
            for name, record in sorted(records.items()):
                if name == "nodejs":
                    continue
                version = record["version"]
                build = record.get("build_string", record.get("build", ""))
                if (not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name)
                        or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.+!_-]*", version)
                        or not isinstance(build, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", build)):
                    raise ValueError("Cannot safely pin the existing native stack for GUI-only; use a new environment")
                specs.append(f"{name}={version}={build}")
        self.native_command = native_request(self.manager, self.prefix, specs)
        self.native_plan = solve_native(self.native_command)
        projected = dict(records)
        for entry in self.native_plan["UNLINK"]:
            projected.pop(entry["name"], None)
        for entry in self.native_plan["LINK"]:
            projected[entry["name"]] = entry
        old_minor = python.split(".")[:2]
        if projected.get("python", {}).get("version", "").split(".")[:2] != old_minor:
            raise ValueError("This extension would change the Python minor version; use a new profile environment")
        self.projected_native = {name: (record["version"], record.get("build_string", record.get("build", "")))
                                 for name, record in projected.items()}
        if args.gui_only:
            for name, record in records.items():
                original = (record["version"], record.get("build_string", record.get("build", "")))
                if name != "nodejs" and self.projected_native.get(name) != original:
                    raise ValueError(f"GUI-only plan would change existing native package {name}; "
                                     "choose a full Desktop plan or a new environment instead")
        if not args.gui_only and "openbabel" in projected:
            raise ValueError("Native plan introduces conda Open Babel; refusing a wheel/native overlay")
        self.constraints = self.directory / "native-constraints.txt"
        mapping = NATIVE_PYPI_NAMES
        constraints = []
        if not args.gui_only:
            names = {mapping.get(p["spec"].split()[0], p["spec"].split()[0])
                     for p in packages if p.get("source") != "pip" and p.get("module")}
            names.update(installed_names)
            names.update(mapping.get(name, name) for name, record in projected.items()
                         if record.get("python_package") or any(
                             dependency.split()[0] in {"python", "python_abi"}
                             for dependency in record.get("depends", []) if isinstance(dependency, str) and dependency.split()))
            for name, entry in sorted(projected.items()):
                mapped = mapping.get(name, name)
                if mapped.lower().replace("_", "-") not in names:
                    continue
                if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", mapped):
                    continue
                if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.+!_-]*", entry["version"]):
                    raise ValueError(f"Native version for {mapped} cannot safely constrain pip; no changes made")
                constraints.append(f"{mapped}=={entry['version']}")
        self.constraints.write_text("\n".join(constraints) + "\n")
        pip_specs = [] if args.gui_only else [p["spec"] for p in packages if p.get("source") == "pip"]
        flags = ["--no-deps"] if args.gui_only else ["--constraint", str(self.constraints)]
        editable_flags = ["-e"] if self.editable else []
        self.pip_command = self.manager.run(self.prefix, ["python", "-m", "pip", "install", "--upgrade", *flags,
                                                       *pip_specs, *editable_flags, str(self.repository)])
        self.pip_plan = self.preview_pip("pip-preview")
        print("\nNative package changes (UNLINK old, LINK new):", flush=True)
        for action in ("UNLINK", "LINK"):
            for entry in self.native_plan[action]:
                print(f"  {action}: {entry['name']} {entry['version']} {entry['build']}", flush=True)
        if not any(self.native_plan.values()):
            print("  No native package changes.", flush=True)
        print("Pip packages proposed (rechecked after native installation):", flush=True)
        for name, version in sorted(self.pip_plan):
            print(f"  {name} {version}", flush=True)
        print("Python minor version is preserved. This is an in-place, two-stage operation, not an atomic rollback.\n"
              "If a solver reports incompatibility, use a new environment. Close GUI/jobs before applying.", flush=True)
        if args.gui_only:
            print("GUI-only: scientific dependencies and profile ownership will not be changed.", flush=True)
            if self.editable:
                print(f"Editable DVBFixer mode will be retained, using this checkout: {self.repository}", flush=True)

    def preview_pip(self, label: str) -> set[tuple[str, str]]:
        report = self.directory / f"{label}.json"
        log = self.directory / f"{label}.log"
        with stage("Preview pip dependencies (no packages installed)", log), log.open("w") as stream:
            result = subprocess.run([*self.pip_command, "--dry-run", "--report", str(report)], env=child_environment(),
                                    stdout=stream, stderr=subprocess.STDOUT, check=False)
            if result.returncode:
                raise ValueError("Pip dependency preview failed; inspect the output above. No pip packages were installed.")
        return pip_installations(report)

    def retain(self, logs: Path) -> None:
        self.logs = logs
        for path in self.directory.iterdir():
            # pip reports may contain authenticated download URLs; retain a
            # reduced identity-only plan, never raw report URLs.
            if path.suffix != ".json":
                shutil.copy2(path, logs / path.name)
        (logs / "extension-plan.json").write_text(json.dumps({
            "profile": None if self.args.gui_only else self.args.profile, "gui_only": self.args.gui_only, "addons": self.args.addon,
            "prefix": str(self.prefix), "native": self.native_plan, "pip": sorted(self.pip_plan),
            "before": self.before, "state": "confirmed-not-completed",
        }, indent=2) + "\n")

    def verify_before_apply(self) -> None:
        if fingerprint(self.prefix) != self.before or solve_native(self.native_command) != self.native_plan:
            raise ValueError("Environment or native plan changed after preview; rerun and confirm the new plan")

    def verify_pip_after_native(self) -> None:
        try:
            actual_native = {name: (record["version"], record.get("build_string", record.get("build", "")))
                             for name, record in native_records(self.prefix).items()}
            if actual_native != self.projected_native:
                raise ValueError("Applied native packages differ from the confirmed plan; stopping before pip and profile publication")
            if not self.args.gui_only:
                # Compact conda dry-run LINK rows may omit Python dependency
                # metadata. Identify newly installed Python distributions now,
                # so pip cannot overlay a native transitive package unnoticed.
                sites = [str(path) for path in self.prefix.glob("lib/python*/site-packages")]
                names = {(d.metadata.get("Name") or "").lower().replace("_", "-")
                         for d in package_metadata.distributions(path=sites)}
                mapping = NATIVE_PYPI_NAMES
                with self.constraints.open("a") as stream:
                    for name, (version, _) in sorted(self.projected_native.items()):
                        name = mapping.get(name, name)
                        if name.lower().replace("_", "-") not in names:
                            continue
                        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.+!_-]*", version):
                            raise ValueError("Native distribution metadata cannot safely constrain pip")
                        stream.write(f"{name}=={version}\n")
            actual = self.preview_pip("pip-after-native")
            if not actual <= self.pip_plan:
                raise ValueError("Pip plan changed after native installation; no unconfirmed pip packages were installed. "
                                 "Native changes may already be present; rerun --extend to review the new plan")
            with self.constraints.open("a") as stream:
                for name, version in sorted(actual):
                    stream.write(f"{name}=={version}\n")
        finally:
            if self.logs is not None:
                log = self.directory / "pip-after-native.log"
                if log.exists():
                    shutil.copy2(log, self.logs / log.name)
                shutil.copy2(self.constraints, self.logs / "applied-pip-constraints.txt")

    def commit(self) -> None:
        if self.marker.is_symlink():
            raise ValueError("Ownership record changed during extension; refusing to replace it")
        current = self.marker.read_bytes() if self.marker.exists() else None
        if current != self.original_marker:
            raise ValueError("Ownership record changed during extension; refusing to replace it")
        if not self.args.gui_only:
            record = {"schema_version": 1, "prefix": str(self.prefix), "profile": self.args.profile,
                      "addons": sorted(set(self.args.addon))}
            temporary = self.marker.with_name(f".{self.marker.name}.{uuid.uuid4().hex}")
            try:
                with temporary.open("x") as stream:
                    json.dump(record, stream, indent=2)
                if self.original_marker is None:
                    os.link(temporary, self.marker)
                else:
                    temporary.replace(self.marker)
            finally:
                temporary.unlink(missing_ok=True)
        if self.logs is not None:
            (self.logs / "extension-completed.json").write_text(json.dumps({
                "profile": None if self.args.gui_only else self.args.profile,
                "gui_only": self.args.gui_only, "checks_completed": True,
            }, indent=2) + "\n")


def prepare_extension(args: argparse.Namespace, manager: Manager) -> ExtensionPlan:
    plan = ExtensionPlan(args, manager)
    try:
        plan.preview()
        return plan
    except BaseException as exc:
        (plan.directory / "preview-error.log").write_text(f"{type(exc).__name__}: {exc}\n")
        plan.close()
        raise


def prepare_extension_interactive(args: argparse.Namespace, manager: Manager, *, interactive: bool) -> ExtensionPlan | None:
    """Recover an incompatible full-profile choice without silently changing scope."""
    try:
        return prepare_extension(args, manager)
    except ProfileStackConflict as exc:
        if not interactive or args.gui_only or not package_details(args.prefix)[2]:
            raise
        from dvbfixer.install_ui import choose

        print(f"\nThe complete profile cannot be attached safely: {exc}\n"
              "No packages or profile ownership were changed. Adding only the GUI does not require replacing this stack.",
              flush=True)
        action = choose("How would you like to continue?", [
            ("gui", "Add GUI only — retain Open Babel/scientific dependencies; skip requested scientific add-ons/license changes"),
            ("cancel", "Stop without changes; create a separate environment for the complete profile"),
        ], default="cancel")
        if action == "cancel":
            return None
        args.gui_only, args.profile, args.addon, args.modeller_key_fd = True, "desktop", [], None
        print("\nScope changed explicitly to GUI only. A new dependency preview and confirmation follow.", flush=True)
        return prepare_extension(args, manager)
