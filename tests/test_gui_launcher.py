"""Desktop launch contracts, including a real supervised fake-server process."""

from __future__ import annotations

import json
import os
import signal
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest

from dvbfixer import gui_launcher as gui
from dvbfixer.desktop_shortcuts import _desktop_quote, install_desktop_shortcut

ROOT = Path(__file__).resolve().parents[1]


def desktop(tmp_path: Path) -> Path:
    prefix = tmp_path / "desktop"
    root = prefix / "share/dvbfixer/gui"
    for name, content in (("dist/index.html", "<title>DVBFixer</title>"),
                          ("dist-server/server.js", "// fake test server"), ("package.json", '{}')):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    (prefix / "bin").mkdir()
    (prefix / "bin/python").symlink_to(sys.executable)
    node = prefix / "bin/node"
    node.write_text(f"#!{sys.executable}\n" + '''
import http.server, json, os, signal, sys, time
from pathlib import Path
root = Path(os.environ["DVBFIXER_GUI_ROOT"])
if os.environ.get("FAKE_GUI_FAILURE"):
    print("deliberate startup failure", file=sys.stderr, flush=True)
    sys.exit(7)
time.sleep(float(os.environ.get("FAKE_GUI_DELAY", "0")))
server = http.server.HTTPServer(("127.0.0.1", int(os.environ["DVBFIXER_PORT"])), http.server.BaseHTTPRequestHandler)
with (root / "starts").open("a") as stream: stream.write("started\\n")
def stop(*_):
    (root / "stopped").write_text("graceful")
    sys.exit(0)
signal.signal(signal.SIGTERM, stop)
fd = int(os.environ["DVBFIXER_DESKTOP_READY_FD"])
os.write(fd, (json.dumps({"version": 1, "pid": os.getpid(), "host": "127.0.0.1", "port": server.server_port}) + "\\n").encode())
os.close(fd)
server.serve_forever()
''')
    node.chmod(0o755)
    return prefix


def command(prefix: Path, state: Path, *args: str) -> list[str]:
    code = (
        "import sys; from pathlib import Path; from dvbfixer import gui_launcher as g; "
        f"sys.prefix={str(prefix)!r}; g.state_directory=lambda *_:Path({str(state)!r}); "
        f"g.webbrowser.open=lambda url: Path({str(prefix / 'browser-opened')!r}).write_text(url); "
        f"g.main({list(args)!r})"
    )
    return [sys.executable, "-c", code]


def environment(**extra: str) -> dict[str, str]:
    return {**{key: value for key, value in os.environ.items() if not key.startswith("DVBFIXER_")},
            "PYTHONPATH": str(ROOT / "src"), **extra}


def wait_state(path: Path, process: subprocess.Popen[str]) -> dict:
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise AssertionError(process.communicate())
        try:
            value = json.loads((path / "instance.json").read_text())
            if value.get("version") == 1:
                return value
        except (OSError, ValueError):
            pass
        time.sleep(0.02)
    raise AssertionError("launcher never became ready")


def test_installed_assets_never_fall_back_to_current_directory(tmp_path, monkeypatch):
    prefix = desktop(tmp_path)
    monkeypatch.chdir(prefix / "share/dvbfixer/gui")
    with pytest.raises(gui.GuiLaunchError, match="Standard does not include"):
        gui.installed_gui(tmp_path / "standard")
    assert gui.installed_gui(prefix) == prefix / "share/dvbfixer/gui"


def test_environment_is_bound_to_selected_prefix_and_preserves_auth(tmp_path, monkeypatch):
    prefix = desktop(tmp_path)
    monkeypatch.setenv("DVBFIXER_AUTH_PRINCIPALS", "private-auth-configuration")
    monkeypatch.setenv("DVBFIXER_EXECUTABLE", "/other/bin/python")
    monkeypatch.setenv("PYTHONPATH", "/other/source")
    root = gui.installed_gui(prefix)
    env = gui.launcher_environment(prefix, root, root / "structures", 0, 7)
    assert env["PATH"].split(os.pathsep)[0] == str(prefix / "bin")
    assert env["DVBFIXER_EXECUTABLE"] == str(prefix / "bin/python")
    assert env["DVBFIXER_HOST"] == "127.0.0.1"
    assert env["DVBFIXER_AUTH_PRINCIPALS"] == "private-auth-configuration"
    assert "PYTHONPATH" not in env
    monkeypatch.setenv("DVBFIXER_HOST", "0.0.0.0")
    with pytest.raises(gui.GuiLaunchError, match="loopback-only"):
        gui.launcher_environment(prefix, root, root / "structures", 0, 7)


def test_missing_runtime_and_nonprivate_state_fail_closed(tmp_path):
    prefix = desktop(tmp_path)
    (prefix / "bin/node").unlink()
    with pytest.raises(gui.GuiLaunchError, match="node runtime is missing"):
        gui.installed_gui(prefix)
    directory = tmp_path / "public-state"
    directory.mkdir(mode=0o755)
    directory.chmod(0o755)
    with pytest.raises(gui.GuiLaunchError, match="private"):
        gui._private_directory(directory)
    original = tmp_path / "original"
    original.write_text("user data")
    linked = tmp_path / "linked"
    linked.symlink_to(original)
    with pytest.raises(OSError):
        gui._open_private(linked)
    assert original.read_text() == "user data"


def test_readiness_reuse_and_explicit_stop(tmp_path):
    prefix = desktop(tmp_path)
    state = tmp_path / "runtime"
    process = subprocess.Popen(command(prefix, state), env=environment(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        first = wait_state(state, process)
        with pytest.raises(gui.GuiLaunchError, match="not responding"):
            gui._control({**first, "token": "0" * 64}, "stop")
        for malformed in ({"token": "я" * 64, "action": "stop"}, {"token": first["token"], "action": []}):
            with socket.create_connection(("127.0.0.1", first["control_port"]), timeout=2) as connection:
                connection.sendall(json.dumps(malformed).encode() + b"\n")
                assert json.loads(connection.recv(100))["ok"] is False
        assert process.poll() is None
        browser = prefix / "browser-opened"
        deadline = time.monotonic() + 2
        while not browser.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        assert browser.read_text().startswith("http://127.0.0.1:")
        second = subprocess.run(command(prefix, state, "--no-browser"), env=environment(), capture_output=True, text=True, timeout=10)
        assert second.returncode == 0, second.stderr
        changed = subprocess.run(command(prefix, state, "--no-browser"), env=environment(DVBFIXER_AUTH_PRINCIPALS="private-config"),
                                 capture_output=True, text=True, timeout=10)
        assert changed.returncode == 1
        assert "different authentication settings" in changed.stderr
        assert "private-config" not in changed.stderr
        port = browser.read_text().rsplit(":", 1)[1]
        occupied = subprocess.run(command(prefix, tmp_path / "busy-runtime", "--no-browser", "--port", port,
            "--data-dir", str(tmp_path / "other-data")), env=environment(), capture_output=True, text=True, timeout=10)
        assert occupied.returncode == 1
        assert "server.log" in occupied.stderr
        assert process.poll() is None
        assert json.loads((state / "instance.json").read_text()) == first
        assert (prefix / "share/dvbfixer/gui/starts").read_text().splitlines() == ["started"]
        stopped = subprocess.run(command(prefix, state, "--stop"), env=environment(), capture_output=True, text=True, timeout=10)
        assert stopped.returncode == 0, stopped.stderr
        process.wait(timeout=10)
        assert process.returncode == 0
        assert (prefix / "share/dvbfixer/gui/stopped").read_text() == "graceful"
        assert first["token"] not in process.communicate()[0]
        assert (state / "instance.json").stat().st_mode & 0o077 == 0
    finally:
        if process.poll() is None:
            process.terminate()
            process.wait(timeout=10)


@pytest.mark.parametrize("settings", [{"FAKE_GUI_FAILURE": "1"}, {"FAKE_GUI_DELAY": "5"}])
def test_failed_or_timed_out_startup_never_opens_browser(tmp_path, settings):
    prefix = desktop(tmp_path)
    result = subprocess.run(command(prefix, tmp_path / "runtime", "--startup-timeout", "0.3"),
                            env=environment(**settings), capture_output=True, text=True, timeout=10)
    assert result.returncode == 1
    assert "ERROR:" in result.stderr
    assert "server.log" in result.stderr
    assert not (prefix / "browser-opened").exists()


def test_sigterm_stops_supervised_server(tmp_path):
    prefix = desktop(tmp_path)
    process = subprocess.Popen(command(prefix, tmp_path / "runtime", "--no-browser"), env=environment(),
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        wait_state(tmp_path / "runtime", process)
        process.send_signal(signal.SIGTERM)
        process.wait(timeout=10)
        assert process.returncode == 0
        assert (prefix / "share/dvbfixer/gui/stopped").exists()
    finally:
        if process.poll() is None:
            process.terminate()
            process.wait(timeout=10)


@pytest.mark.parametrize("platform", ["darwin", "linux"])
def test_shortcuts_are_prefix_bound_idempotent_and_do_not_overwrite(tmp_path, monkeypatch, platform):
    prefix = desktop(tmp_path)
    home = tmp_path / "home"
    monkeypatch.delenv("XDG_DATA_HOME", raising=False)
    target = install_desktop_shortcut(prefix, home=home, platform=platform)
    assert install_desktop_shortcut(prefix, home=home, platform=platform) == target
    script = target / "Contents/MacOS/DVBFixer" if platform == "darwin" else prefix / "share/dvbfixer/desktop-launcher"
    assert str(prefix / "bin/python") in script.read_text()
    assert "-m dvbfixer.cli gui" in script.read_text()
    assert os.access(script, os.X_OK)
    if platform == "linux":
        assert "--stop" in target.read_text()
    conflict = target / "Contents/Info.plist" if platform == "darwin" else target
    conflict.write_text("user-owned contents")
    with pytest.raises(gui.GuiLaunchError, match="not overwritten"):
        install_desktop_shortcut(prefix, home=home, platform=platform)
    assert conflict.read_text() == "user-owned contents"


def test_shortcut_symlinks_are_rejected(tmp_path, monkeypatch):
    prefix = desktop(tmp_path)
    monkeypatch.delenv("XDG_DATA_HOME", raising=False)
    script = prefix / "share/dvbfixer/desktop-launcher"
    script.symlink_to(tmp_path / "user-file")
    with pytest.raises(gui.GuiLaunchError, match="symlink"):
        install_desktop_shortcut(prefix, home=tmp_path / "home", platform="linux")


def test_linux_exec_quoting_handles_reserved_characters_without_shell_expansion():
    quoted = _desktop_quote('/prefix/with $value/"quotes"/back\\slash/100%/launcher')
    assert quoted.startswith('"') and quoted.endswith('"')
    assert "\\\\$value" in quoted
    assert '\\\\"quotes\\\\"' in quoted
    assert "back\\\\\\\\slash" in quoted
    assert "100%%" in quoted


@pytest.mark.parametrize("argv", [["--port", "65536"], ["--startup-timeout", "nan"], ["--startup-timeout", "-1"]])
def test_invalid_startup_options_are_rejected(argv):
    with pytest.raises(SystemExit):
        gui.parse_args(argv)
