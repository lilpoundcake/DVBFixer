# `dvbfixer gui`

[User guide](../README.md) · [Installation](../installation.md) · [CLI reference](../reference/gui.md)

Launch the installed Desktop GUI with one command, without remembering Node/npm
paths. This command is not a scientific workflow and cannot be submitted as a
GUI-managed job.

## Quick start

From the active Desktop environment:

```bash
dvbfixer gui
```

The server binds IPv4 loopback on an available port. Its local URL is printed,
and the browser opens only after the server confirms readiness over an inherited
private pipe. Keep this terminal session running. Stop with Ctrl+C or, from
another terminal in the same environment:

```bash
dvbfixer gui --stop
```

Stopping also cancels managed child processes through the server's normal
shutdown path. Closing only the browser does not stop the server.

Repeated launch reopens the running instance for that installation and workspace
data directory. An unrelated service occupying a port is never assumed to be
DVBFixer. No API credential is passed in a browser URL.
Changing authentication settings requires stopping the existing instance first;
the launcher refuses to silently reuse a server with different authentication.

## Application shortcuts

Opt in when installing Desktop:

```bash
python3 scripts/install_dvbfixer.py --profile desktop --manager micromamba \
  --prefix "$HOME/micromamba/envs/dvbfixer-desktop" --desktop-shortcut
```

The macOS shortcut is a user-local `.app` under `~/Applications`; Linux uses
`$XDG_DATA_HOME/applications` or `~/.local/share/applications`. Names include a
prefix-derived identifier so separate installations do not replace each other.
Shortcuts invoke their own environment's Python without shell activation.
The Linux entry also offers a **Stop DVBFixer** action. Native Finder/desktop-menu
acceptance is recorded separately from file-generation and command-line tests.

## Options and data

- `--no-browser`: print the URL without opening a browser.
- `--port 5173`: request a particular local port; an occupied port fails clearly.
- `--data-dir PATH`: use a separate workspace root. Repeat the same option when
  stopping or reopening that instance.
- `--startup-timeout 60`: allow more time for startup (maximum 600 seconds).
- `--check`: check that the selected environment contains GUI assets and Node,
  without starting a server. This is not a scientific-backend readiness test.

By default, workspace files stay in the installed GUI's `structures/` directory
below `<prefix>/share/dvbfixer/gui/`, matching the previous installed npm route.
`DVBFIXER_GUI_DATA_DIR` is also respected. Relative data paths are relative to
the installed GUI, never the terminal's current directory.

Launcher state and `server.log` are under `~/Library/Caches/dvbfixer-desktop/`
on macOS or `$XDG_CACHE_HOME/dvbfixer-desktop/` (default `~/.cache`) on Linux.
The selected log path is printed on startup. Instance-control capabilities are
kept in private user-owned state files, not workspace manifests, logs or URLs.
Configured server authentication is preserved. This launcher refuses remote
binding; use the separately documented [standalone deployment](../../gui/README.md#build-for-production)
for explicitly configured advanced hosting.

## Troubleshooting

- **Standard environment:** it does not contain the GUI. Install Desktop; the
  launcher does not silently use a development checkout or another prefix.
- **Missing Node/build:** repair the Desktop installation with the source helper.
- **Startup failed/port busy:** inspect the printed `server.log` path. The
  browser is not opened for an unready server; omit `--port` for a free port.
- **Browser could not open:** open the printed loopback URL manually.
- **Existing shortcut differs:** it is not overwritten. Inspect that file/app
  before manually replacing it; another installation may own it.

For source development, continue using `npm run dev` in the repository's `gui/`
with a Desktop environment active. The raw manager/npm start command in the
[installation guide](../installation.md) is an advanced diagnostic fallback.

Directory batch input is unsupported. GUI launch does not change scientific
validation rules or qualify optional modeling backends.
