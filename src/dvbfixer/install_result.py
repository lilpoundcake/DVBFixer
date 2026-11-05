"""One final, actionable terminal result for every installer operation."""

from __future__ import annotations

import re
import shlex
from contextvars import ContextVar
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class InstallResult:
    scope: str = "Setup"
    retry: str = "Repeat the same operation and select the same environment."
    prefix: Path | None = None
    logs: Path | None = None
    preview: bool = False
    cancelled: bool = False
    started: bool = False
    stage: str | None = None
    failed_stage: str | None = None
    stage_log: Path | None = None
    errors: list[str] = field(default_factory=list)
    actions: list[str] = field(default_factory=list)


current: ContextVar[InstallResult | None] = ContextVar("installer_result", default=None)


def clean(text: str) -> str:
    text = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", text)
    return re.sub(r"(https?://)[^/\s@]+@", r"\1[redacted]@", text)


def configure(args) -> None:
    result = current.get()
    if result is None:
        return
    result.prefix = args.prefix.expanduser().absolute() if args.prefix is not None else None
    result.preview = args.dry_run
    if args.runtime:
        result.scope = f"Additional backend {args.runtime} (provisioning, not scientific qualification)"
        command = ["python3", "scripts/install_dvbfixer.py", "--runtime", args.runtime, "--resume"]
        if args.manager:
            command.extend(["--manager", args.manager])
        if result.prefix is not None:
            command.extend(["--prefix", str(result.prefix)])
        result.retry = f"After resolving the listed actions, continue this runtime: {shlex.join(command)} (interactive setup follows)."
    elif args.gui_only:
        result.scope = "GUI only (scientific stack unchanged)"
        result.retry = "Repeat Add profile/GUI > GUI only, or --extend --gui-only, for this environment."
    elif args.package_only:
        result.scope = "DVBFixer package only (dependencies and GUI assets not upgraded)"
        result.retry = "Repeat the package-only update (--update --package-only) in this environment."
    elif args.extend:
        result.scope = f"DVBFixer {args.profile or 'profile'} extension/update"
        result.retry = "Repeat --extend, choose the same environment and review a fresh profile plan."
    else:
        result.scope = f"DVBFixer {args.profile or 'installation'}"
        result.retry = "For an interrupted helper-recorded installation choose --resume and select the same installation."


def log_directory(path: Path) -> None:
    if result := current.get():
        result.logs = path


def started() -> None:
    if result := current.get():
        result.started = True


def cancelled() -> None:
    if result := current.get():
        result.cancelled = True


def failure(exc: BaseException | str) -> None:
    if result := current.get():
        message = clean(str(exc))
        if message and message not in result.errors:
            result.errors.append(message)


def stage_started(label: str, log: Path | None) -> None:
    if result := current.get():
        if result.failed_stage is None:
            result.stage, result.stage_log = label, log


def stage_failed(label: str, log: Path | None = None) -> None:
    if result := current.get():
        # Preserve the innermost failure, not an enclosing orchestration stage.
        if result.failed_stage is None:
            result.failed_stage, result.stage_log = label, log


def readiness(code: int, log: Path) -> None:
    result = current.get()
    if result is None or code == 0:
        return
    result.stage_log = log
    result.failed_stage = "Installation readiness checks"
    if log.is_file():
        for line in log.read_text(errors="replace").splitlines():
            text = clean(line).strip()
            if re.match(r"(?:error|action-required):", text, re.IGNORECASE):
                target = result.actions if code == 2 else result.errors
                if text not in target:
                    target.append(text)


def required(actions: list[str]) -> None:
    if result := current.get():
        result.actions.extend(clean(action) for action in actions if action not in result.actions)


def remedies(text: str, result: InstallResult) -> list[str]:
    value = text.lower()
    hints = []
    if "pdbfixer-resources" in value or "soft.xml" in value:
        hints.append("Preserve this damaged prefix and create a fresh complete profile; do not copy individual PDBFixer resources into it.")
    elif "openbabel" in value or "open babel" in value:
        hints.append("Keep conda Open Babel and openbabel-wheel separate. Choose GUI only to retain this stack, or a new environment for the full profile.")
    if "modeller" in value or "license" in value:
        hints.append("Configure an appropriate MODELLER license using hidden entry or the protected key descriptor; never put the key in command-line arguments.")
    if any(word in value for word in ("checkpoint", "--kalign", "auxiliary", "frozen sha256")):
        hints.append("Use the managed download after reviewing the model terms, or supply the exact local checkpoint/config and other listed frozen assets, then resume this runtime.")
    if "launcher" in value and any(word in value for word in ("already exists", "kept unchanged", "replace-runner")):
        hints.append("Resume and explicitly confirm shared launcher replacement (--replace-runner for unattended use), or choose another --runner-directory; existing launchers are never silently overwritten.")
    if any(word in value for word in ("clobber", "path conflict", "path collision", "safety checks")):
        hints.append("Do not use --clobber or disable safety checks. Preserve the diagnostics and choose a clean prefix/supported manager route.")
    if any(word in value for word in ("unsatisfiable", "resolutionimpossible", "dependency conflict", "different python", "requires-python", "no matching distribution", "incompatible python")):
        hints.append("Use the declared compatible Python/dependency profile or a new environment; do not relax the project's version constraints to force the solve.")
    if any(word in value for word in ("no space", "disk quota", "enospc")):
        hints.append("Free disk space/check the quota, then retry; keep the failed prefix and diagnostics until inspected.")
    if any(word in value for word in ("permission denied", "not writable", "access is denied")):
        hints.append("Choose a user-writable named environment and check its ownership; do not force a shared/system prefix to become writable.")
    if any(word in value for word in ("connection", "could not resolve", "certificate", "sslerror", "timed out")):
        hints.append("Check network/proxy/certificate configuration and retry; do not disable TLS or checksum verification.")
    if any(word in value for word in ("node", "npm", "gui build", "build desktop gui")):
        hints.append("Review the GUI build log and rerun GUI-only/Desktop setup using its selected manager's Node runtime (22+); do not use an unrelated global Node installation.")
    if "ownership record" in value or "installation record" in value:
        hints.append("Keep the existing record for inspection. Do not overwrite/adopt another installation or replace a changed record to bypass the check.")
    if not hints:
        hints.append("Inspect the failed step's log below. Retain it when reporting a reproducible installer failure; do not assume partial package changes were rolled back.")
    hints.append(result.retry)
    return list(dict.fromkeys(hints))


def finish(result: InstallResult, code: int) -> None:
    print("\n=== Installation result ===", flush=True)
    if result.cancelled or code == 130:
        print("CANCELLED — operation stopped.", flush=True)
    elif code == 0 and result.preview:
        print("PREVIEW COMPLETED — no target environment installation was applied.", flush=True)
    elif code == 0:
        print(f"SUCCESS — {result.scope}. Requested operation completed.", flush=True)
    elif code == 2 and result.actions:
        print("ACTION REQUIRED — setup is not fully ready.", flush=True)
    elif code == 2:
        print("ACTION REQUIRED — check the required setup findings above; this is not a completed licensed installation.", flush=True)
    else:
        print("ERROR — installation did not complete successfully.", flush=True)
    if result.prefix is not None:
        print(f"Installation: {result.prefix}", flush=True)
    if code not in (0, 130):
        # Action-required is not a subprocess failure. Successful version snapshots
        # contain arbitrary package metadata, not evidence of a network problem.
        step = result.failed_stage or (result.stage if code != 2 else None)
        if step:
            print(f"Uncompleted step: {step}", flush=True)
        elif code == 2:
            print("Pending: required setup actions below; no failed command is reported.", flush=True)
        findings = list(dict.fromkeys([*result.errors, *result.actions]))
        tail = ""
        if code != 2 and result.stage_log is not None and result.stage_log.is_file():
            # Native solver logs can be large; summary generation must be bounded.
            with result.stage_log.open("rb") as stream:
                stream.seek(0, 2)
                stream.seek(max(0, stream.tell() - 16384))
                lines = stream.read().decode("utf-8", errors="replace").splitlines()
            tail = "\n".join(clean(line) for line in lines[-15:])
            if not findings or all("returned non-zero exit status" in line for line in findings):
                meaningful = [line.strip() for line in tail.splitlines() if line.strip()]
                findings.extend(meaningful[-5:])
        for finding in findings:
            print(f"  - {finding}", flush=True)
        print("Next steps:", flush=True)
        for hint in remedies("\n".join([*findings, tail, step or ""]), result):
            print(f"  - {hint}", flush=True)
    if result.started and code != 0:
        print("Some steps may already have changed this installation; no automatic rollback is promised.", flush=True)
    if result.stage_log is not None and (code != 2 or result.failed_stage is not None):
        print(f"Step log: {result.stage_log}", flush=True)
    if result.logs is not None:
        print(f"Diagnostics: {result.logs}", flush=True)
