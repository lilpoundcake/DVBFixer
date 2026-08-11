# DVBfixer GUI

Local graphical workspace for DVBfixer, protein structures, multiple alignment,
and homology modeling. It runs in the browser with a Node-side development
backend and an optional PostgreSQL-backed mutations table.

The backend invokes `dvbfixer` from `PATH` by default. Configure an alternate
binary with `DVBFIXER_EXECUTABLE` and an optional fixed JSON argument array in
`DVBFIXER_ARGS`; subprocesses have bounded captured output and a configurable
timeout (`DVBFIXER_MAX_OUTPUT_BYTES` and `DVBFIXER_TIMEOUT_MS`).

## What it does

- **3D viewer (×2)** -- interactive [Mol*](https://molstar.org) molecular
  visualization with custom residue-type coloring. Open a second viewer to
  compare structures side-by-side; toggle camera sync via the link icon in
  the top bar.
- **Sequence viewer** -- amino acid sequence with color-coded residue types,
  click or drag to select ranges
- **Pairwise alignment** -- Needleman-Wunsch with BLOSUM62 between any two
  chains, including across two different loaded structures
- **Bidirectional selection** -- click a residue in sequence or alignment, it
  highlights in 3D (solid ball-and-stick) and vice versa
- **Elements tree** -- hierarchical view of polymers, ligands, ions, and water
  with visibility toggles + per-chain **Show Interface** button (5 Å contact
  zone, distinguishes polymer vs ligand on the same chain id)
- **Interactions table** -- computed H-bonds, ionic, disulfide, hydrophobic,
  pi-stacking, and more with chain pair filtering; auto-filters when Show
  Interface is active
- **DVBfixer panel** -- every current CLI tool, generated from argparse and
  grouped by workflow and argument section; runs persist across panel reloads,
  stream status, capture bounded logs, and can be cancelled
- **Homology workspace** -- persistent target/template projects, MAFFT,
  MUSCLE 5, or Clustal Omega MSA, manual gap editing, license-free Biopython
  structural fitting grouped by target chain (with optional Modeller SALIGN), template
  span masks, and model generation
- **Workspaces** -- each workspace owns its files, generated runs/logs,
  active A/B structures, workflow inputs/results, and Homology workflows
- **Workspace file browser** -- alternating rows, Shift/Cmd multi-selection,
  drag-handle and arrow reordering, PDB-only filtering, rename/download/trash
  context actions, and a dedicated read-only Text Files tab
- **Mutations panel** -- editable DataGrid backed by PostgreSQL for keeping
  antibody mutation sets (e.g. YTE, LS, DLE)
- **Library** -- an ordered list of workspaces; the active workspace's files
  appear in the separate Workspace panel
- **Dockable panels** -- drag panels to rearrange, "+" button to spawn any
  panel into any tabset
- **Local files** -- upload your own .pdb or .mmcif files into either viewer

Everything except DVBFixer runs and the Mutations DB runs in the browser.
Your files never leave your machine.

### Workspace data and Homology selection

For the complete target-to-model workflow, selection controls, structural
fitting behavior, template-plan schema, diagnostics, and limitations, see the
[GUI Homology Workspace guide](../docs/gui-homology.md).

Workspaces live below `structures/projects/<workspace-id>/` in a versioned
`workspace.json` plus `files/`, `runs/`, and `homology/` directories. On first
start, legacy top-level Library folders become workspaces and ungrouped files go
to `Unsorted`; legacy files remain intact as a recovery source.

Internal bookkeeping (`workspace.json`, run manifests, captured stdout/stderr,
and temporary alignment/model inputs) stays on disk but is hidden from the
workspace file list and workflow file selectors. Right-click a workspace or visible file
to move it to trash. Workspace trash is stored under
`structures/_workspace_trash/`; file trash stays in the owning workspace's
`.trash/` directory.

Library workspaces and Workspace files can be reordered from their drag
handles; the insertion marker shows the future position. The Workspace toolbar
can hide non-PDB files without changing their stored order. A file context menu
downloads the original file, while a workspace download produces a `.tar.gz`
archive. Structure metadata is stored on the workspace artifact and edited in
Info with revision-aware autosave.

The Homology Target tab parses FASTA/PIR or protein chains from a workspace
PDB/mmCIF file. **+ Add new** creates a template selector, using the active
viewer structure when available. With multiple target chains, Alignment
shows a target-chain selector and one synchronized, horizontally scrolling
alignment at a time.

| Consensus mark | Meaning |
|---|---|
| `*` | Fully conserved column |
| `:` | Strongly similar residues |
| `.` | Weakly similar residues |
| blank | Mismatch or a column containing a gap |

Click a template residue to select it, Shift-click for a range, and Ctrl/Cmd-
click to add or remove disjoint residues. Selecting a template loads it into
viewer A when necessary and, while linked, synchronizes the 3D highlight. This
selection is transient: **Use selection as modeling span** explicitly copies
it into the persisted model mask. Gap arrows move a gap without changing row
length; the toolbar adds or removes full gap columns.

## Install & run

### Minimum (viewer only — no DVBFixer, no Mutations DB)

```
git clone <repo-url>
cd DVBFixer/gui
npm install --legacy-peer-deps
npm run dev:no-db
```

Open http://localhost:5173. The Mutations tab will show a configuration
hint; everything else works.

### Full install (DVBFixer + Mutations)

DVBfixer GUI drives the package in this repository and can store antibody
mutation sets in PostgreSQL.

**1.** Create the repository's DVBfixer environment and install the package:

```
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
```

#### Modeller license

`model` and (indirectly) other DVBFixer commands rely on
[Modeller](https://salilab.org/modeller/), which requires a free academic
license. Register at <https://salilab.org/modeller/registration.html>, then:

```
micromamba activate dvbfixer
KEY=YOUR_LICENSE_KEY bash gui/scripts/set-modeller-key.sh
```

The helper finds `<env>/lib/modeller-*/modlib/modeller/config.py` and writes
`license = r'YOUR_LICENSE_KEY'` into it (backing up the previous file).
You can also pass an explicit prefix as the first arg:

```
KEY=YOUR_LICENSE_KEY bash scripts/set-modeller-key.sh /opt/conda/envs/tarantino
```

**2.** Start the app. PostgreSQL is auto-managed — if Docker is installed,
`npm run dev` spins up a local postgres container (port 5432) and sets
`DATABASE_URL` automatically. The `mutations` table is created on first
connection.

```
cd gui
npm install --legacy-peer-deps
micromamba activate dvbfixer   # so `dvbfixer` is on PATH
npm run dev
```

The first run downloads the postgres image; the container survives `Ctrl+C`
so subsequent runs are instant.

### DB controls

| Command             | What it does                                          |
|---------------------|-------------------------------------------------------|
| `npm run dev`       | Auto-starts postgres (if docker present), then vite   |
| `npm run dev:no-db` | Skip auto-postgres, just run vite                     |
| `npm run db:up`     | Start the postgres container                          |
| `npm run db:down`   | Stop & remove the container (volume persists)         |
| `npm run db:logs`   | Tail postgres logs                                    |

**Override** the auto-setup by exporting `DATABASE_URL` yourself before
`npm run dev` — the script detects an existing value and skips Docker:

```
export DATABASE_URL=postgres://my-user:pw@my-host:5432/my-db
npm run dev
```

**Override DVBFixer** if it's not on PATH (e.g. wrapped in micromamba):

```
export DVBFIXER_EXECUTABLE="micromamba"
export DVBFIXER_ARGS='["run", "-n", "tarantino", "dvbfixer"]'
```

The DVBFixer tab is usable even without the env (it will just error on
Run); the Mutations tab is usable even without DATABASE_URL or Docker
(it will show a configuration message).

## Tech stack

React 19, TypeScript 6, Vite 6, Mol*, MUI (Material UI v9, plus
`@mui/x-data-grid`), flexlayout-react, Zustand. PostgreSQL via `pg`
(loaded lazily; optional). Vite middleware backend (`server/api-plugin.ts`).

## Panels

| Panel               | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| 3D Structure        | Mol* viewer with custom SCSS skin and residue-type color theme on sticks    |
| 3D Structure (B)    | Independent secondary viewer for side-by-side comparison; camera sync toggle in top bar |
| Sequence            | Monospace amino acid grid, drag-select, independent chain selector per tab  |
| Alignment           | Pairwise Needleman-Wunsch (BLOSUM62) across chains, incl. across structure A and B; click number row to pick the same column on both sides |
| Elements            | Categorized tree (polymer/ligand/ion/water) with per-component visibility + per-chain Show Interface |
| Interactions        | Computed non-covalent + covalent contacts, filterable by type and chain pair |
| DVBfixer            | Generated CLI forms with managed queued/running state, bounded logs, restore-after-reload, and Cancel |
| Homology            | Target → templates → editable MSA and masks → Modeller project workflow |
| Mutations           | Editable DataGrid backed by PostgreSQL (`mutations` table: chain / mutation_name / mutations) |
| Library             | Ordered workspace list with rename, archive download, and recoverable trash |
| Workspace           | Active-workspace files with A/B loading, import, filtering, reordering, download, rename, and recoverable trash |
| Info                | Editable artifact metadata and structure statistics with revision-aware autosave |

Every panel can be duplicated via the "+" button on its tabset header.
Sequence panels maintain independent chain selections.

## DVBfixer commands

| Command   | What it does                                                                              |
|-----------|-------------------------------------------------------------------------------------------|
| split     | Empirical chain splitting via residue numbering / peptide-bond distance / nearest-atom gaps |
| renumber  | Renumber residues by aligning ATOM records to SEQRES; removes Kabat-style insertion codes |
| model     | Rebuild missing loops / gaps with Modeller (LoopModel + MD refinement)                     |
| prepare   | Add missing residues, heavy atoms, hydrogens via PDBFixer; supports point mutations        |
| minimize  | Energy-minimize with OpenMM (AMBER14 + GLYCAM); optional xtb / obminimize post-refine     |
| protonate | Predict per-residue pKa (PROPKA3); rename to AMBER protonation variants at target pH       |
| convert   | Convert structure naming between PDB/AMBER/GLYCAM and CHARMM conventions. Default (`--to-amber`): PDB/CHARMM → GLYCAM sugars + AMBER protonation variants (HID/HIE/HIP, ASH/GLH, LYN, CYX/CYM). With `--to-charmm`: reverse direction (GLYCAM/AMBER → CHARMM, NLN/OLS/OLT revert to ASN/SER/THR, ROH/OME caps dropped). The two flags are mutually exclusive. `--no-roh` applies only to the `--to-amber` direction. Was named `glycam` in older DVBFixer versions. |

The complete command and flag schema is generated by
`../scripts/gen_gui_spec.py` from each command's argparse parser. This keeps the
GUI synchronized with new DVBfixer releases, including repeatable and
multi-value options.

Commands run as workspace-scoped managed jobs. The panel shows queued,
running, succeeded, failed, or cancelled state; restores an active job after a
panel/page reload; streams status with polling fallback; and displays captured
stdout/stderr when the job finishes. Only one job can run in a workspace at a
time. **Cancel** stops the active process, and successful outputs are added to
Workspace and loaded into viewer A when they are structures.

## Workspace management

The **Library** contains workspaces only. Create, activate, reorder, rename,
download, or move a workspace to recoverable trash there. The separate
**Workspace** panel lists files owned by the active workspace. Import local
files from its header or the global Upload button; generated DVBfixer and
Homology outputs are registered in the same manifest automatically.

When two 3D viewers are open, the `A` / `B` selector in Workspace chooses the
destination for the next structure-file click. Internal run records and logs
remain on disk for reproducibility but stay hidden from ordinary file pickers.

`workspace.json` is the active artifact and metadata source. On first use,
older top-level and per-workspace `index.json` files are imported into workspace
manifests without deleting the legacy files; they remain recovery sources only.
The retired `/structures/index.json` scanner and `/api/library/*` mutation
routes are not served.

## Controls

| Action                       | How                                                          |
|------------------------------|--------------------------------------------------------------|
| Activate workspace           | Click a workspace in Library                                 |
| Load workspace structure     | Click a structure file in Workspace                          |
| Choose target viewer         | `A` / `B` selector in the Workspace header                   |
| Filter workspace files       | **Hide non-PDB files** / **Show non-PDB files**              |
| Reorder workspaces/files     | Drag the row handle to the insertion marker; arrows move selected files |
| Download                     | Context menu → **Download**; workspaces download as `.tar.gz` archives |
| Load your own file           | Click Upload in the top bar                                  |
| Toggle camera sync (A ↔ B)   | Link icon in the top bar                                     |
| Rotate 3D                    | Left-click drag                                              |
| Zoom                         | Scroll wheel                                                 |
| Select residue (3D)          | Click an atom                                                |
| Select residue (seq)         | Click a letter in the Sequence panel                         |
| Select range (seq)           | Click and drag across residues                               |
| Pick aligned column (both A and B) | Click a position in the number row above/below an alignment row |
| Drag-select on alignment     | Click and drag across one side                               |
| Toggle element visibility    | Eye icon in the Elements panel                               |
| Show Interface               | Network/hub icon next to a chain in Elements (zooms to 5 Å contact zone, filters Interactions) |
| Focus interaction            | Click a row in the Interactions panel                        |
| Run a DVBFixer command       | Pick the sub-tab, set the input file, click Run              |
| Cancel a DVBFixer command    | Click **Cancel** while the job is queued or running          |
| Duplicate a panel            | Click "+" on any tabset header                               |
| Clear 3D markers             | Tap empty space in 3D                                        |
| Clear everything             | Press Escape                                                 |

## Build for production

```
npm run build
```

Output goes to `dist/`. Serve it with any static file server:

```
npx serve dist
```

Note: the production build is **viewer-only**. The DVBFixer and Mutations
backends are dev-time Vite middleware in `server/api-plugin.ts` — they
don't ship in the static build. To run those in production, host a Node
server that re-uses the plugin (or port the routes).

## Notes

- Water is hidden by default when a structure loads
- The custom color theme (carbons by residue type, non-carbons CPK) applies to focus/stick representations
- Non-canonical amino acids (phosphorylated residues, modified cysteines, protonation variants, etc.) are normalized to standard codes
- Failed workflow output is kept out of the visible workspace manifest and moved to a workspace-scoped failure directory when supported
- The `postinstall` script (`scripts/fix-native-deps.mjs`) handles cross-platform rollup native bindings for macOS and Linux
