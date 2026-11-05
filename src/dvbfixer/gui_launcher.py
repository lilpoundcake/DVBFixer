"""Launch the installed Desktop GUI without manager/npm command composition."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import hmac
import json
import math
import os
import re
import secrets
import selectors
import signal
import socket
import stat
import subprocess
import sys
import time
import webbrowser
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

from dvbfixer.batch import add_runtime_help


class GuiLaunchError(RuntimeError):
    """An actionable startup/instance-management failure."""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="dvbfixer gui", description=__doc__)
    parser.add_argument("--no-browser", action="store_true", help="Print the local URL without opening a browser")
    parser.add_argument("--port", type=int, default=0, help="Loopback port (default: choose a free port)")
    parser.add_argument("--data-dir", type=Path, help="Workspace data directory (default: installed GUI structures/)")
    parser.add_argument("--startup-timeout", type=float, default=30.0, help="Startup/readiness timeout in seconds")
    operation = parser.add_mutually_exclusive_group()
    operation.add_argument("--stop", action="store_true", help="Stop this installation's local GUI and its managed jobs")
    operation.add_argument("--check", action="store_true", help="Check installed Desktop assets and Node without starting")
    add_runtime_help(parser)
    args = parser.parse_args(argv)
    if not 0 <= args.port <= 65535:
        parser.error("--port must be between 0 and 65535")
    if not math.isfinite(args.startup_timeout) or not 0 < args.startup_timeout <= 600:
        parser.error("--startup-timeout must be positive and at most 600 seconds")
    return args


def installed_gui(prefix: Path) -> Path:
    """Never fall back to a checkout or an unrelated environment."""
    root = prefix / "share/dvbfixer/gui"
    for relative in ("dist/index.html", "dist-server/server.js", "package.json"):
        if not (root / relative).is_file():
            raise GuiLaunchError(
                "Desktop GUI is not installed in this environment. "
                "Install the Desktop profile with scripts/install_dvbfixer.py "
                "(Standard does not include the GUI). Missing: " + str(root / relative)
            )
    for executable in ("node", "python"):
        if not (prefix / "bin" / executable).is_file() or not os.access(prefix / "bin" / executable, os.X_OK):
            raise GuiLaunchError(f"Desktop {executable} runtime is missing: {prefix / 'bin' / executable}; repair the Desktop installation")
    return root


def launcher_environment(prefix: Path, root: Path, data: Path, port: int, ready_fd: int) -> dict[str, str]:
    from dvbfixer.install import child_environment

    env = child_environment()
    if env.get("DVBFIXER_HOST", "127.0.0.1") not in {"127.0.0.1", "localhost", "::1"}:
        raise GuiLaunchError("dvbfixer gui is loopback-only; unset DVBFIXER_HOST for local Desktop launch")
    env.update({
        "PATH": str(prefix / "bin") + os.pathsep + env.get("PATH", ""),
        "CONDA_PREFIX": str(prefix), "CONDA_DEFAULT_ENV": str(prefix), "CONDA_SHLVL": "1",
        "DVBFIXER_GUI_ROOT": str(root), "DVBFIXER_STATIC_DIR": str(root / "dist"),
        "DVBFIXER_GUI_DATA_DIR": str(data), "DVBFIXER_HOST": "127.0.0.1", "DVBFIXER_PORT": str(port),
        "DVBFIXER_EXECUTABLE": str(prefix / "bin/python"),
        "DVBFIXER_ARGS": json.dumps(["-m", "dvbfixer.cli"]),
        "DVBFIXER_DESKTOP_READY_FD": str(ready_fd),
    })
    # Conda's AmberTools activation normally supplies this. Shortcuts do not
    # activate a shell and must not inherit a different installation's value.
    if (prefix / "bin/tleap").is_file():
        env["AMBERHOME"] = str(prefix)
    return env


def _private_directory(path: Path) -> None:
    path.mkdir(mode=0o700, parents=True, exist_ok=True)
    info = path.lstat()
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
        raise GuiLaunchError(f"Launcher state directory must be user-owned, private and not a symlink: {path}")


def state_directory(prefix: Path, data: Path) -> Path:
    if sys.platform == "darwin":
        base = Path.home() / "Library/Caches"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache"))).expanduser()
        if not base.is_absolute():
            raise GuiLaunchError("XDG_CACHE_HOME must be an absolute path")
    key = hashlib.sha256(f"{prefix}\0{data}".encode()).hexdigest()[:24]
    return base / "dvbfixer-desktop" / key


def _open_private(path: Path, *, append: bool = False) -> int:
    fd = os.open(path, os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW | (os.O_APPEND if append else 0), 0o600)
    info = os.fstat(fd)
    if not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077 or info.st_nlink != 1:
        os.close(fd)
        raise GuiLaunchError(f"Launcher file must be private, user-owned and not linked: {path}")
    return fd


def _write_state(fd: int, state: dict[str, Any]) -> None:
    payload = json.dumps(state).encode()
    os.lseek(fd, 0, os.SEEK_SET)
    os.ftruncate(fd, 0)
    os.write(fd, payload)
    os.fsync(fd)


def _read_state(fd: int) -> dict[str, Any]:
    os.lseek(fd, 0, os.SEEK_SET)
    try:
        state = json.loads(os.read(fd, 8192))
        if (not isinstance(state, dict) or state.get("version") != 1
                or type(state.get("control_port")) is not int
                or not 1 <= state["control_port"] <= 65535
                or not isinstance(state.get("token"), str) or len(state["token"]) != 64):
            raise ValueError("invalid state")
        return state
    except (ValueError, KeyError, TypeError) as exc:
        raise GuiLaunchError("GUI is starting or its private state is unavailable; retry shortly") from exc


def _socket_message(connection: socket.socket, limit: int) -> Any:
    payload = b""
    while b"\n" not in payload:
        chunk = connection.recv(min(512, limit + 1 - len(payload)))
        if not chunk:
            raise ValueError("incomplete control message")
        payload += chunk
        if len(payload) > limit:
            raise ValueError("oversized control message")
    return json.loads(payload.split(b"\n", 1)[0])


def _control(state: dict[str, Any], action: str) -> dict[str, Any]:
    try:
        with socket.create_connection(("127.0.0.1", state["control_port"]), timeout=2) as connection:
            connection.sendall(json.dumps({"token": state["token"], "action": action}).encode() + b"\n")
            response = _socket_message(connection, 4096)
            if not isinstance(response, dict) or response.get("ok") is not True:
                raise ValueError("control rejected")
            return response
    except (OSError, ValueError) as exc:
        raise GuiLaunchError("The existing GUI launcher is not responding; inspect its local launcher log") from exc


def _open_browser(url: str, enabled: bool) -> None:
    print(f"DVBFixer GUI: {url}", flush=True)
    if enabled:
        try:
            opened = webbrowser.open(url)
        except webbrowser.Error:
            opened = False
        if not opened:
            print("WARNING: Could not open a browser automatically; open the local URL above.", file=sys.stderr)


def _auth_fingerprint(environment: Mapping[str, str]) -> str:
    settings = {key: environment.get(key) for key in ("DVBFIXER_AUTH_PRINCIPALS", "DVBFIXER_LEGACY_WORKSPACE_OWNER")}
    return hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()


def _wait_ready(fd: int, process: subprocess.Popen[bytes], timeout: float, stopped: Callable[[], bool] = lambda: False) -> str:
    deadline = time.monotonic() + timeout
    payload = b""
    with selectors.DefaultSelector() as selector:
        selector.register(fd, selectors.EVENT_READ)
        while time.monotonic() < deadline:
            if stopped():
                raise GuiLaunchError("GUI startup cancelled")
            if process.poll() is not None:
                raise GuiLaunchError("GUI server exited before readiness")
            if not selector.select(min(0.1, max(0, deadline - time.monotonic()))):
                continue
            chunk = os.read(fd, 4096)
            if not chunk:
                raise GuiLaunchError("GUI server closed its readiness channel before startup")
            payload += chunk
            if len(payload) > 4096:
                raise GuiLaunchError("GUI server sent an invalid readiness message")
            if b"\n" in payload:
                try:
                    ready = json.loads(payload.split(b"\n", 1)[0])
                    port = ready["port"]
                    if (ready["version"] != 1 or ready["pid"] != process.pid or ready["host"] != "127.0.0.1"
                            or type(port) is not int or not 1 <= port <= 65535):
                        raise ValueError("invalid readiness")
                except (ValueError, KeyError, TypeError) as exc:
                    raise GuiLaunchError("GUI server sent an invalid readiness message") from exc
                return f"http://127.0.0.1:{port}"
    raise GuiLaunchError("Timed out waiting for GUI readiness")


def _shutdown(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    process.send_signal(signal.SIGTERM)
    try:
        process.wait(timeout=20)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait(timeout=5)


def launch(prefix: Path, args: argparse.Namespace) -> int:
    import fcntl

    root = prefix / "share/dvbfixer/gui"
    if not args.stop:
        installed_gui(prefix)
        if os.environ.get("DVBFIXER_HOST", "127.0.0.1") not in {"127.0.0.1", "localhost", "::1"}:
            raise GuiLaunchError("dvbfixer gui is loopback-only; unset DVBFIXER_HOST for local Desktop launch")
    raw_data = args.data_dir or Path(os.environ.get("DVBFIXER_GUI_DATA_DIR", str(root / "structures")))
    data = (root / raw_data.expanduser()).resolve()
    directory = state_directory(prefix, data)
    _private_directory(directory)
    lock_fd = _open_private(directory / "instance.json")
    try:
        try:
            fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            deadline = time.monotonic() + args.startup_timeout
            while True:
                try:
                    state = _read_state(lock_fd)
                    reply = _control(state, "stop" if args.stop else "status")
                    break
                except GuiLaunchError:
                    if time.monotonic() >= deadline:
                        raise
                    time.sleep(0.1)
            if args.stop:
                print("GUI shutdown requested; the server will stop managed jobs before exiting.")
            else:
                if reply.get("auth_sha256") != _auth_fingerprint(os.environ):
                    raise GuiLaunchError("GUI is already running with different authentication settings; stop it before changing authentication")
                url = reply.get("url", "")
                if not isinstance(url, str) or not re.fullmatch(r"http://127\.0\.0\.1:[1-9][0-9]{0,4}", url):
                    raise GuiLaunchError("Existing launcher returned an invalid local URL")
                if args.port and url != f"http://127.0.0.1:{args.port}":
                    raise GuiLaunchError("This GUI is already running on another port; stop it before changing --port")
                _open_browser(url, not args.no_browser)
            return 0
        if args.stop:
            print("No GUI is running for this installation/data directory.")
            return 0
        root = installed_gui(prefix)
        _write_state(lock_fd, {})
        return _serve(prefix, root, data, directory, lock_fd, args)
    finally:
        os.close(lock_fd)


def _serve(prefix: Path, root: Path, data: Path, directory: Path, lock_fd: int, args: argparse.Namespace) -> int:
    stop = False

    def request_stop(_signal: int, _frame: Any) -> None:
        nonlocal stop
        stop = True

    previous = {sig: signal.signal(sig, request_stop) for sig in (signal.SIGINT, signal.SIGTERM)}
    process: subprocess.Popen[bytes] | None = None
    read_fd, write_fd = os.pipe()
    try:
        env = launcher_environment(prefix, root, data, args.port, write_fd)
        log = directory / "server.log"
        with os.fdopen(_open_private(log, append=True), "ab") as stream, socket.socket() as control:
            control.bind(("127.0.0.1", 0))
            control.listen(8)
            control.settimeout(0.2)
            print(f"GUI launcher log: {log}", flush=True)
            process = subprocess.Popen([str(prefix / "bin/node"), str(root / "dist-server/server.js")],
                cwd=root, env=env, stdout=stream, stderr=stream, pass_fds=(write_fd,), start_new_session=True)
            os.close(write_fd)
            write_fd = -1
            try:
                url = _wait_ready(read_fd, process, args.startup_timeout, lambda: stop)
            except GuiLaunchError as exc:
                raise GuiLaunchError(f"{exc}; see {log}") from exc
            os.close(read_fd)
            read_fd = -1
            token = secrets.token_hex(32)
            auth_sha256 = _auth_fingerprint(env)
            _write_state(lock_fd, {"version": 1, "control_port": control.getsockname()[1], "token": token})
            _open_browser(url, not args.no_browser)
            print("Stop with Ctrl+C or dvbfixer gui --stop (same data directory).", flush=True)
            while not stop and process.poll() is None:
                try:
                    connection, _ = control.accept()
                except TimeoutError:
                    continue
                with connection:
                    connection.settimeout(1)
                    try:
                        message = _socket_message(connection, 512)
                        candidate = message.get("token") if isinstance(message, dict) else None
                        authorized = isinstance(candidate, str) and candidate.isascii() and len(candidate) == 64 and hmac.compare_digest(candidate, token)
                        if not authorized or message.get("action") not in ("status", "stop"):
                            connection.sendall(b'{"ok":false}\n')
                            continue
                        connection.sendall(json.dumps({"ok": True, "url": url, "auth_sha256": auth_sha256}).encode() + b"\n")
                        stop = message["action"] == "stop"
                    except (OSError, ValueError):
                        continue
            if not stop and process.returncode:
                raise GuiLaunchError(f"GUI server stopped unexpectedly; see {log}")
        return 0
    finally:
        try:
            if process is not None:
                _shutdown(process)
        finally:
            if read_fd >= 0:
                os.close(read_fd)
            if write_fd >= 0:
                os.close(write_fd)
            _write_state(lock_fd, {})
            for sig, handler in previous.items():
                signal.signal(sig, handler)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if sys.platform not in {"darwin", "linux"}:
            raise GuiLaunchError("Desktop launch currently supports macOS and Linux")
        prefix = Path(sys.prefix).resolve()
        if args.check:
            print(f"Desktop GUI ready: {installed_gui(prefix)}")
            return 0
        return launch(prefix, args)
    except (GuiLaunchError, OSError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        main()
