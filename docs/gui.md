# GUI guide

[User guide](README.md) · [Installation](installation.md) · [GUI Homology](gui-homology.md)

The GUI provides workspace files, two Mol* viewers, and forms for DVBfixer tools.
For server configuration rather than daily use, see the [GUI README](../gui/README.md).

## Start the local GUI

Use [Desktop](installation.md#standarddesktop-source-checkout-helper) to install
the scientific environment, Node and built GUI together. Start that installation
with `micromamba run -p /path/to/desktop npm --prefix /path/to/desktop/share/dvbfixer/gui start`.
For source development, activate Desktop before starting the server so it can
find `dvbfixer` and external tools. Node.js 22 or newer is required.
From the repository root:

```bash
micromamba activate dvbfixer-desktop
cd gui
npm ci
npm run dev
```

The development server opens the browser. Keep the terminal running while using
the application. For a built local server, see [production build](../gui/README.md#build-for-production).

## 1. Create or choose a workspace

**Library** lists workspaces. **Workspace** shows files belonging to the active
workspace. A workspace also retains tool settings, generated outputs, and Homology
workflows; switching workspaces changes that context.

Use **Import files** or drop files onto Workspace. Import PDB or CIF/mmCIF
structures and, when needed, FASTA sequences or supporting files.

With a remote server, imported files are uploaded to that server. With the default
local server, they stay on the same machine. Do not expose a development server
publicly; administrators should follow the [deployment guide](deployment.md).

## 2. Find and open files

- Category icons distinguish structures, reports, and other artifacts.
- Format badges reflect the real file path, not an editable display name.
- Hover over a row for its path and recorded provenance.
- **Show structures only** retains PDB, CIF/mmCIF, and other recognized structure
  formats. **Show all files** restores the remaining files without changing order.
- Load a structure into viewer A or B. Each viewer shows its own filename;
  use the second viewer to compare input and output.
- Open text-like artifacts in **Text Files**. Reports and logs are separate from
  coordinate files.

## 3. Run a tool

Open the **DVBfixer** panel, choose a command, and select its workspace input.
Use the [command index](commands/index.md) to decide which tool fits your task.
Set the options, run the command, and inspect its status and output.

Most tools run as managed jobs with persistent status and cancellation.
**Atom Names** uses a synchronous conversion request instead. Its optional
**Infer Protonation** checkbox is off by default and reads existing hydrogen
patterns; it does not predict pKa or rebuild hydrogens. **Dry Run** analyzes
without publishing a converted PDB. See [Atom Names](commands/atom-names.md).

## 4. Review and save results

Read the tool's report and run log, then load the resulting structure into the
other viewer. A completed job can still have scientific warnings: inspect
diagnostic findings and the affected regions before using the output.

Use a file's context menu to download it. Download a workspace from Library
as a `.tar.gz` archive. Rename and reorder files without relying on display names
to identify their format. Shift and Cmd/Ctrl support multi-selection.

Trash actions are recoverable; files move to workspace trash rather than being
immediately erased. Internal manifests and temporary helper files are hidden
from ordinary file pickers.

### Compact run storage

Workspace shows **one collapsible group per run**, rather than a flat list of
every intermediate file. Expand a group to open its outputs, summary and diagnostic
archive. The download button on the group saves all archived diagnostics at once.
Imported source files remain independent rows. The structures-only filter also
applies inside run groups.

Click a run's label to select its files; Ctrl/Cmd adds or removes groups/files,
and Shift selects a range. The arrow independently expands or collapses the run
without clearing selection. Right-click a group (including a collapsed one) for
its context menu; right-clicking an already selected group or file preserves the
current multi-selection. Trash actions apply to the selected files, not to the
whole workspace or unregistered run bookkeeping.

All managed DVBfixer commands and Homology modeling, MSA and structural fitting
use the same storage policy:

- Declared results, numbered output candidates and companion files remain
  directly accessible at their original paths. A diffusion publication bundle
  remains intact, including its digest-bound sidecars and reports.
- Additional requests, masks, motifs and intermediate structures are compressed
  into `diagnostics.tar.gz`, with an embedded file-size/SHA-256 manifest. Archive
  verification completes before the unpacked diagnostic copies are removed.
- `run-summary.json` records the status, outputs, archive inventory and any
  packaging warning. Job records and stdout/stderr logs remain available to the
  run history; runtime cache/home/scratch trees are not copied into publications.
- Failed and cancelled runs retain diagnostic evidence; their partial structures
  are not published as usable outputs. Scientifically rejected diffusion
  candidates remain labelled and are never automatically loaded as accepted ones.
- If packing fails, diagnostics remain unpacked and the summary records the
  warning. Storage quotas and access controls still apply.

The dedicated Atom Names and Altlocs APIs already publish bounded single-result
artifacts; they retain their existing contracts and appear in the same run groups.
Read-only analysis/dry runs do not create extra output artifacts. Existing run
files are grouped visually only: this does not migrate or delete historical data.

For offline inspection, download and extract `diagnostics.tar.gz` with a standard
archive utility. No automatic archive extraction is performed by the server.

## Alternate locations

Use **DVBfixer → Altlocs** to analyze a structure, choose a policy and review
exactly one checkbox per alternate residue. The embedded Mol* preview focuses
the selected residue without replacing the main viewers. Apply creates a new
single-state PDB and retains a selection report; source coordinates and occupancy
are unchanged. See [altlocs](commands/altlocs.md) for policy and connectivity rules.

## Homology modeling

For target sequences, templates, alignment editing, selected template regions,
and model generation, follow the [GUI Homology workflow](gui-homology.md).
Start with MODELLER; experimental diffusion requires separately installed runners
and has profile-specific limits.

## If something fails

1. Read the selected tool's log and diagnostic report.
2. Run `dvbfixer doctor` in the environment used to start the server.
3. If an executable was just installed, restart the server to refresh its `PATH`.
4. Check [troubleshooting](known-issues.md) and the relevant command guide.
