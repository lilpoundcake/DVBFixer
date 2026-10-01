# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
npm install          # Install dependencies (runs postinstall: rollup native bindings)
npm run dev          # Start the Vite development server
npm run build        # Type-check and build dist/ plus dist-server/server.js
npm start            # Run the built loopback standalone GUI + API server
npm run typecheck    # Type-check only
npm run lint         # ESLint
npm run preview      # Frontend-only Vite preview (no APIs)
```

## What This Is

Tarantino is a mostly-local browser-based protein structure viewer with an
optional Node-side dev backend (Vite middleware) for DVBFixer pipeline runs
It loads PDB/mmCIF files and provides
a dockable multi-panel workspace:

The Library lists workspaces. `server/workspace-api.ts` owns versioned
manifests below `structures/projects/<id>/`; `stores/workspaceStore.ts` is the
frontend source of truth for the active workspace. Imports and every workflow
request must carry `workspaceId` and resolve paths inside that workspace. Do not
reintroduce global `/structures/index.json` selectors in workflow panels.
Artifacts marked `hidden` are reproducibility helpers and must not appear in
the workspace file list or workflow selectors. Workspace/file deletion is recoverable:
the API moves targets to `_workspace_trash/` or the workspace's `.trash/` rather
than unlinking them.

- **3D Structure (primary + optional secondary)**: two independent Mol*
  viewers with optional camera sync between them
- **Sequence**: amino acid sequence with residue-type coloring, drag-to-select,
  per-panel chain selection; SEQRES residues missing from ATOM coords are
  rendered greyed-out + dashed-border + italic and are not interactive
- **Alignment**: pairwise Needleman-Wunsch (BLOSUM62) alignment between any
  two chains, including across two different loaded structures
- **Elements**: tree of polymers / ligands / ions / water with per-component
  visibility + per-chain **Show Interface**; clicking a row focuses the
  camera on that element
- **Interactions**: computed H-bonds, ionic, cation-pi, pi-stacking, halogen,
  hydrophobic, metal coordination, disulfide, covalent
- **DVBFixer**: form-driven UI for the DVBFixer CLI (split / renumber / model /
  prepare / minimize / protonate / convert), outputs registered as child entries
- **Library**: ordered list of workspaces; files live in the separate Workspace panel
- **Info**: stats summary at the top, single-field metadata (Name + Notes),
  and an **Equivalent chains** section that auto-groups multimeric copies
  via sequence identity (with optional manual override persisted on the
  active workspace artifact)
- **Settings**: app-wide preferences (auto-orient-on-load toggle + Alignment
  source-label toggle [File / Name]; all persisted in `localStorage`)

Selecting / hovering residues in 3D highlights them in Sequence and Alignment
(bidirectional sync). The Alignment panel routes selection per-side to whichever
viewer (primary or secondary) holds that chain.

## Tech Stack

- **React 19** with **TypeScript 6**, bundled by **Vite 6**
- **MUI (Material UI)** for UI components — no Tailwind, no shadcn
- **flexlayout-react** for dockable panels
- **Mol\*** (`molstar` npm package, used directly — not `pdbe-molstar`)
- **Zustand** for state management
- **Sass** for Mol* SCSS skin

## Architecture

### Panel System

`flexlayout-react` `Layout` + `Model` in `App.tsx`. Default layout (in the
same tabset, the first tab is the active one):
- Left column: Library, then (Info | **Settings**).
- Right column (main viewer): (3D Structure | **DVBFixer** | **Homology**),
  with (Sequence | **Alignment**) and
  (Elements | Interactions | Clashes) tabsets below.

Every tabset has a "+" button (`onRenderTabSet`) that opens a MUI Menu
listing: 3D Structure, 3D Structure (B), Sequence, Elements, Interactions,
Alignment, DVBFixer, Homology, Library, Workspace, Info, and Settings. Sequence panels keep their own chain selection.

### Data Flow

```
Workspace files / FileLoader ─┐
                              ├─→ primary plugin   → extractChains/Elements/Meta → structureStore (chains, elements, meta)
                              └─→ secondary plugin → extractChains              → structureStore (secondaryChains)
                                                ↓
            structureStore + selectionStore ────┼────────────────────────────────────┐
                                                │                                    │
                       SequenceViewer  ElementsTable  InteractionsPanel  Info        │
                                                ↑                                    │
                       AlignmentPanel  ←  primaryChains + secondaryChains            │
                                                ↑                                    │
                              MolstarViewer (primary + secondary slots)  ←  helpers ─┘
```

### Dual 3D Viewers

`MolstarViewer` takes `slot: 'primary' | 'secondary'` and each renders its own
independent `PluginUIContext`. Both plugins are kept in `structureStore`
(`plugin`, `secondaryPlugin`). Workspace files and FileLoader honor `loadTargetSlot`
('primary' | 'secondary') to choose where to load. The secondary viewer
publishes only its chains (for cross-structure Alignment) — it doesn't touch
elements / meta / Interactions.

**Tab-close cleanup** — when the user closes a 3D Structure tab,
`MolstarViewer`'s effect-cleanup disposes the plugin AND clears, for that
slot, the store's plugin/chains/fileName. Without the `fileName` clear,
the Workspace panel's A/B chip would remain stuck on whichever structure was last
    loaded into that viewer even though the viewer no longer exists.
Additionally, when the SECONDARY slot is torn down and
`loadTargetSlot === 'secondary'`, the cleanup snaps it back to
`'primary'` — otherwise the Workspace panel's A/B toggle would silently stay on
B (its UI is gated on `secondaryPlugin`, so it disappears with the
viewer), every subsequent workspace-file click would bail out with *"Open a
    '3D Structure (B)' tab first"*, and the user would be trapped unable to
    load anything into A.

Post-load, each viewer (a) hides the water component, (b) swaps the **ion**
component's default ball-and-stick representation for `spacefill` so each ion
renders as a Van der Waals sphere — `createStructureRepresentationParams(plugin,
structure, { type: 'spacefill' })` is built once and `.update()`ed into every
ion-component repr cell, then `runTask(updateTree())`. The default `physical`
size theme uses each element's VdW radius. (c) **optionally** runs
`PluginCommands.Camera.OrientAxes` to face the principal axis — gated on
`useStructureStore.autoOrientOnLoad` (default OFF, toggled in the Settings
panel, persisted in `localStorage` under key `tarantino.autoOrientOnLoad`).
When ON, camera sync is temporarily disabled during the orient so the
snapshot doesn't get mirrored to the other viewer mid-orientation.

### Camera Sync (`src/hooks/useCameraSync.ts`)

Toggled via the link icon in the AppBar (only shown when both viewers exist).
Mechanism: subscribe to each plugin's `canvas3d.didDraw`, snapshot the camera,
mirror with `dst.canvas3d.camera.setState(snap, 0)` for instant zero-anim
mirroring.

**Anti-feedback uses pending-draw flags, NOT snapshot equality**. When subA
mirrors A → B, it sets `pendingFromAToB = true`. The next `B.didDraw` consumes
the flag and skips its mirror back to A. Critical: value-based equality
(`lastAppliedToB`) DOES NOT work during a focus animation — A's snap updates
every frame, so by the time B's echo-draw fires `subB`, `lastAppliedToB` has
already been overwritten by A's newer frame, the equality fails, and the echo
back to A interrupts A's in-flight animation. The pending-flag pattern is
value-independent and survives multi-frame animations.

**Skip-mirror-from-empty-viewer**. The mirror also returns early when
`src.managers.structure.hierarchy.current.structures.length === 0`. An
empty viewer has a meaningless default camera (tiny radiusMax, origin
target); mirroring it onto a viewer that DOES have a structure freezes
the destination at a useless zoom level. Common scenario it fixes:
user opens a structure on A while B is still empty — B's idle didDraw
events would otherwise push B's default state onto A.

### Camera ownership: `manualReset: true`

`MolstarViewer` post-init calls `PluginCommands.Canvas3D.SetSettings(plugin,
{ settings: { camera: { manualReset: true } } })`. This tells Mol*'s
`canvas3d.commitScene` to NEVER call `resolveCameraReset` on its own.
Without it, every state-tree update that grows the visible bounding sphere
(adding sticks, glycam representations, etc.) would queue an
auto-reset, which fires on the NEXT draw and overwrites any in-flight
`camera.focusSphere` animation — the user sees "camera tries to jump and
snap back".

Trade-off: with `manualReset: true`, Mol*'s built-in "fit camera to scene
on first structure load" is also suppressed. We compensate in
`MolstarViewer` post-load:
- `autoOrientOnLoad === true`  → `PluginCommands.Camera.OrientAxes(plugin)`
- `autoOrientOnLoad === false` → `plugin.managers.camera.reset(undefined, 0)`

Both branches suppress camera sync during the fit so the other viewer
doesn't mirror the fit motion.

### Workspace Library and files

`ProjectLibrary.tsx` renders two active panels: Library is the ordered list
of workspaces, while Workspace shows the active manifest's visible artifacts.
Both use dedicated drag handles and persist ordering through the workspace API.
Artifact metadata lives directly in `workspace.json` and is edited through the
revisioned artifact-metadata endpoint. Legacy `structures/index.json` data is
used only by the one-time workspace migration; there is no legacy structure
tree frontend.

### Mol* Integration

Initialized via `createPluginUI` + `renderReact18` in `MolstarViewer.tsx`.
Key API patterns:

- **Load**: `plugin.builders.data.rawData()` → `parseTrajectory()` → `hierarchy.applyPreset()`
- **Query residues**: `MolScriptBuilder` → `compile()` + `QueryContext` → `StructureSelection.toLociWithSourceUnits()`
- **Select / highlight**: `plugin.managers.interactivity.lociSelects.select/deselectAll()` / `lociHighlights.highlight/clearHighlights()`
- **Read selection**: `plugin.managers.structure.selection.entries`; iterate atoms via `OrderedSet.getAt()` / `OrderedSet.size()`; read props via `StructureProperties.*(location)` on a `StructureElement.Location`
- **Custom-tagged repr** (NOT `focus.setFromLoci`): `StateTransforms.Model.StructureSelectionFromBundle` + `StateTransforms.Representation.StructureRepresentation3D` with explicit tags. See `showSticksForLoci()` in `molstar-helpers.ts`.
- **Visibility**: `setSubtreeVisibility()` from `molstar/lib/mol-plugin/behavior/static/state`
- **Interactions**: `computeInteractions()` from `molstar/lib/mol-model-props/computed/interactions`
- **Camera**: `plugin.managers.camera.focusSphere(Loci.getBoundingSphere(loci))` / `focusLoci(loci)` / `reset(undefined, 0)`

**Why we avoid `focus.setFromLoci()`**: `StructureFocusRepresentation`'s
`focus.clear()` is guarded by `if (this.state.current)` and short-circuits
when undefined; `ensureShape` unconditionally creates a `SurrSel` sub-tree
with `expandRadius` (default 5 Å, can't be 0). Workaround: push `undefined`
directly onto `plugin.managers.structure.focus.behaviors.current` and use
our own tagged repr nodes.

### Custom Tagged Representations + Color Theme

`src/lib/molstar-helpers.ts`:
- `showSticksForLoci(plugin, loci, tag, label)` — tagged `StructureSelectionFromBundle` + `StructureRepresentation3D` (ball-and-stick, **`xrayShaded: false`** for solid sticks). Replaces existing nodes with that tag.
- `deleteCellsByTag(plugin, tag)` — `StateSelection.Generators.root.subtree().withTag(tag)` → delete.
- `focusInterfaceForChain(plugin, chainId, category, radius=5)` — entity-type-aware contact MolScript (chain X residues within 5 Å of any non-self atom ∪ partner-side residues within 5 Å of chain X). Uses tag `tarantino-interface`, then `camera.focusSphere`. The `category` arg maps to a MolScript `entity.type` test so a polymer chain "A" and a ligand with chain id "A" are distinct selves.
- `showSelectionSticks(plugin, residues)` / `clearSelectionSticks(plugin)` — tag `tarantino-selection`. Called by `useSequenceSync` and `AlignmentPanel`.
- `showSurroundingsAndFocus(plugin, residues, radius=5)` / `clearSurroundings(plugin)` — tag `tarantino-focus-surr`. Builds `MS.struct.modifier.includeSurroundings({ target: residues, radius, as-whole-residues })`, renders the combined region as sticks, AND calls `camera.focusSphere` ONCE. Order: camera FIRST, then fire-and-forget sticks (so the state-tree commit can't pre-empt the camera request). Used by the Sequence panel's Zoom button and the Alignment panel's Focus button. The deterministic single camera move (combined with `manualReset: true` on the canvas3d) prevents the "tries to jump, snaps back" bug previously caused by `focus.setFromLoci` triggering its own auto-pan.
- `buildResiduesLoci(plugin, residues)` — synchronously build a Loci from a residue list without any state-tree mutation. Use when you need the loci immediately (e.g. for `camera.focusLoci`) without round-tripping through `selection.getLoci` (which is async-lagged behind the most recent `lociSelects.select`).

`src/lib/residue-color-theme.ts` registers `tarantino-residue-type`: carbons
by residue class (hydrophobic green, positive blue, negative red, polar orange,
cysteine yellow, aromatic teal, special pink), other atoms CPK. Applied to
the default focus representation (via subscription that retags
`structure-focus-target-repr` / `-surr-repr` cells) and our custom-tagged sticks.

`src/molstar-theme.scss` overrides the Mol* SCSS skin to match the app palette.
`src/lib/residue-codes.ts` has `THREE_TO_ONE` for non-canonical codes (CYX/CYM/CSO → C,
HID/HIE/HIP → H, SEP → S, MSE → M, etc.); `toCanonicalThree()` maps back.

### Sync Hooks

- `useMolstarSync` (3D → sequence): registers a `LociMarkProvider` on primary
  plugin's `lociSelects` (wrapped in try/catch, **zero re-entrant side-effects**
  so a throw can't break the next provider in the chain). Subscribes to primary
  `interaction.hover` (debounced 50 ms). **Empty-click cleanup is attached to
  BOTH viewers' `interaction.click`** — on `Loci.isEmpty` it wipes 3D state in
  both viewers (deselect / clearHighlights / focus.next(undefined) + clear /
  clearInterfaceFocus / clearSelectionSticks / `camera.reset(undefined, 0)`),
  then `setFocusedChain(null)`, `selectionStore.clearSelection()`, and
  `fireClearAll()` so the Alignment panel resets its local sel sets.
- `useSequenceSync` (sequence → primary 3D): on `_lock === 'sequence'`, async
  clearSelection → push `undefined` to focus → `await clearInterfaceFocus` →
  if non-empty `selectResiduesInViewer('select')` + `await showSelectionSticks`;
  else `clearSelectionSticks`.
- `useCameraSync` (primary ↔ secondary): see above.

`selectionStore._lock` (expires after 200 ms) prevents infinite update loops.

### Sequence Panel (`src/components/SequenceViewer.tsx`, `src/components/ChainSelector.tsx`)

`extractChains` in `molstar-helpers.ts` walks ATOM records to collect residues
per chain, then merges in **SEQRES residues missing from ATOM** by looking up
`model.sequence.byEntityKey[entityIndex]` (via `model.entities.getEntityIndex(entityId)`)
and walking its `seqId.value(i)` / `compId.value(i)` columns. Each residue gets:
- `seqId: number` — the canonical `label_seq_id` (1-based sequential per entity).
  Used by every MolScript / selection / 3D-sync path in the codebase.
- `authSeqId?: number | null` — the PDB author residue number (`auth_seq_id`,
  what the structure file literally says, e.g. 250, 251, 252...). Null for
  SEQRES-only residues (no ATOM coordinates → no authored number available).
  Displayed in the Sequence panel's "Structure" numbering mode.
- `present: boolean` — `true` if the residue has coordinates, `false` if SEQRES-only.

`SequenceViewer`, `ChainSelector`, **and** `AlignmentPanel`'s chain picker all
apply the same filter pipeline to the store's `chains`:
1. Strip water/ion residues (HOH/ZN/MG/...) from each chain's residue list.
2. Drop chains with ≤ 1 residue left.
3. Drop chains where **every** residue's `threeToOne()` is `'X'` — typically
   glycans (NAG / BMA / MAN) or other non-polypeptide chains that got assigned
   a chain id by Mol*'s polymer-classifier. Without this filter the dropdown
   shows chains whose "sequence" is just a row of X's.
4. Sort chains alphabetically via `localeCompare(..., { numeric: true, sensitivity: 'base' })`.

`SequenceViewer`'s chain-init effect **validates every fallback candidate** against
its filtered `chains` list. The store's `activeChainId` is set to `chains[0].id` of
the *unfiltered* chains, which can be a glycan filtered out by step 3 — using it
blindly would silently leave `activeChain` undefined and render a blank pane.
Fix: try `[initialChainId, globalChainId]` in order, accept the first that exists
in the filtered list, otherwise fall back to filtered `chains[0].id`.

Missing residues render with greyed text (`#b5bfcc`), dashed border, italic font,
and `cursor: not-allowed`; mouse handlers are nulled out so they can't be selected
or hovered. Tooltip reads `"<RES> PDB <auth> · #<order> (<class>) — declared in SEQRES, missing from structure"`.

**Numbering toggle.** A two-button `ToggleButtonGroup` in the toolbar (per-panel
state, default **Structure**) swaps between:
- **Structure** — show `authSeqId` (PDB author residue number, e.g. 250 …
  with insertion gaps preserved). Falls back to `seqId` for SEQRES-only residues.
- **Sequence** — show the 1-based ordinal position within the visible chain
  (1, 2, 3 …). Counts every residue including SEQRES-missing ones so the
  gutter aligns with rendered cells.
Tooltip always shows BOTH numbers (`PDB 322 · #205`) so users can cross-reference
without flipping the toggle. Selection / 3D sync continue to use the canonical
`seqId` (`label_seq_id`) regardless of toggle state.

### Info Panel + Equivalent Chains (`src/components/StructureInfo.tsx`, `src/lib/chain-grouping.ts`)

`StructureInfo.tsx` layout (top → bottom):

1. **Summary** — 4 stat cards (Chains / Residues / Atoms / Elements).
   The Chains and Residues counts are computed against
   `filterSequenceableChains(chains)` so glycan / water / ion-only chains
   don't inflate the numbers.
2. **File** — file path in a monospace `Paper` with
   `wordBreak: 'break-all'` so long DVBFixer output paths
   (`dvb_<cmd>_<ts>/<input>_<cmd>.pdb`) wrap instead of overflowing.
3. **Metadata** — three fields: `Name` (free text), **IgG Subtype**
   (singleSelect: `'' | IgG1 | IgG2 | IgG3 | IgG4 | IgA | IgM | IgE | IgD`),
   **Allotype** (free text, e.g. `G1m17,1` / `nG1m1`). The legacy fields
   (Organism / Method / Resolution) are still in `StructureMeta` +
   `index.json` for backwards compat, but their inputs are hidden.
4. **Notes** — multi-line `description` TextField.
5. **Equivalent chains** — see below.

**Equivalent chains section** (`EquivalentChainsSection` inside
`StructureInfo.tsx`, helper in `src/lib/chain-grouping.ts`):

- **Auto-detect.** `computeEquivalentChains(chains, threshold=0.95)`
  builds a sequence per chain via `chainToSequence` (from
  `src/lib/alignment.ts`), runs `alignSequences` pairwise + computes
  `trimmedIdentity` (also from `lib/alignment.ts`). Pairs with
  trimmed-identity ≥ threshold are merged via union-find (single-linkage
  clustering). Result is sorted: multi-member groups first by smallest
  chain id, singletons last.
- **Trimmed identity.** `trimmedIdentity(result)` walks the alignment,
  finds the first and last columns where **neither** sequence is a gap,
  counts `annotation === '|'` matches in that window. Internal gaps stay
  in the window and count as mismatches. This handles the FcRn case
  where chains H and I are the same protein but I has a truncated
  N-terminus — without trimming, the global identity would be much
  lower than 95% and the pair would never group.
- **Chain filter.** `filterSequenceableChains` is the same pipeline used
  by `SequenceViewer`, `ChainSelector`, and `AlignmentPanel`: strip
  water/ion residues, drop chains with ≤ 1 residue left, drop chains
  whose every residue maps to `'X'` (glycans). Re-used here so the
  groups match the chain set the user sees elsewhere.
- **Display mode.** Multi-member groups render as Paper rows with one
  filled Chip per chain id + a subtle `97.0 % over 218 aa` annotation.
  Singletons collapse into one de-emphasised `Unique: B, F, K` row.
  Footer toolbar: `[Edit groups]` + (when override exists) a `manual`
  Chip + `[Reset]` button.
- **Edit mode.** One TextField per group (comma-separated chain ids) +
  `[+ Add group]`. `[Save]` runs `validateGrouping(parsed, availableIds)`
  which canonicalises (sort + dedupe), reports `duplicates` (chain id in
  multiple groups → error) and `unknown` (chain id not present → error).
  Chains left out become implicit singletons. `[Cancel]` discards local
  edits.
- **Persistence.** `StructureInfo` saves `equivalentChains` through the
  revision-aware workspace artifact metadata endpoint. `undefined` means
  "auto-detect, don't persist"; an array (including `[]`) is a manual override,
  and `null` removes the persisted field. A successful PATCH replaces the
  active manifest with the server's incremented revision; a stale revision is
  shown as a conflict instead of overwriting newer edits.
- **Load path.** Every workspace artifact load applies the artifact's complete
  metadata snapshot through `structureMetaFromArtifact(...)`. This preserves
  manual overrides across structure switches and clears fields omitted by the
  next artifact so metadata cannot leak between structures.

### Alignment Panel (`src/components/AlignmentPanel.tsx`, `src/lib/alignment.ts`)

`alignment.ts` is pure-TS Needleman-Wunsch with BLOSUM62 (alphabet
`ARNDCQEGHILKMFPSTWYVBZX*`, affine gap penalty open -11, extend -1). Returns
aligned strings + `|` / `:` / `.` / ` ` annotation + identity / similarity /
score / length. Also exports `chainToSequence(residues)` (1-letter
sequence from a residue list — used by `AlignmentPanel` and
`chain-grouping`) and `trimmedIdentity(result)` (identity computed after
trimming leading / trailing gap columns — used by `chain-grouping` and
elsewhere when truncated termini shouldn't penalise the score).

`AlignmentPanel.tsx`: two chain pickers (A / B) grouped by source ('A' =
primary viewer chains, 'B' = secondary viewer chains); each source's chains
go through the all-X / length / alphabetic-sort filter described in the
Sequence Panel section. Drag-select per side;
mouseup commits and pushes to the corresponding plugin via
`selectResiduesInViewer` + `showSelectionSticks`. The number row above/below
each sequence is clickable — picking any column toggles residues on BOTH
sides (bilateral pick). `Focus` zooms each viewer to its selection; if
camera sync is on and viewers differ, sync is suppressed during focus.
Subscribes to `structureStore.clearAllSignal` to reset its local sel sets.

### Interactions Panel

`InteractionsPanel.tsx` uses `computeInteractions()` + `structure.interUnitBonds` /
intra-unit bond scans for disulfides and inter-chain covalents. Filterable
table with chain pair dropdowns (water excluded). Clicking a row focuses the
3D view. When `structureStore.focusedChainId` is set (via Elements'
"Show Interface"), table auto-filters and shows a banner:
`Interface: polymer chain A · 24 contacts · ↔ B (12) ↔ C (3) · [Clear]`.

### Clashes Panel (`src/components/ClashesPanel.tsx`, `src/lib/clash-detection.ts`)

Steric-clash detection. `computeClashes()` walks every atom pair within
`2 · maxVdW − minOverlap` Å via `structure.lookup3d.find`, computes VdW
overlap `(rA + rB) − distance` against a Bondi/Rowland-Taylor radius
table, and reports pairs above the threshold. Excludes: H/D, water,
1-2 and 1-3 neighbors across the **full** bond graph (`unit.bonds` PLUS
`structure.interUnitBonds`), and same-residue pairs (rotamer/ring
topology produces false-positives). The `gatherNeighbors(structure,
unit, atomIdx)` helper walks both intra and inter edges so that
glycosidic bonds between sugar non-polymer units (C1–O 1-2 and C1–C2
1-3 around the linkage angle would otherwise show as severe clashes)
and inter-chain disulfides are correctly excluded. Two tiers: `bad`
(0.4–0.9 Å), `severe` (>0.9 Å) — matches PyMOL / ChimeraX / MolProbity
convention.

Panel UI:
- Header: bad/severe count chips, severity filter (All / Bad / Severe),
  **Group** toggle (default OFF), min-overlap numeric input,
  Clear-highlight / Recompute icon buttons.
- **Flat mode** (default): one row per atom-pair clash, sorted worst-first.
- **Grouped mode**: rows collapsed by canonical `(chain,resId)` pair.
  Each group row shows worst severity (`bad × 3` when N > 1), worst
  overlap, residue endpoints, and atom count. Chevron expands the
  group into indented child rows for the per-atom-pair detail.
  Clicking a group row focuses the WORST atom-pair in the group;
  the group stays selected whenever any of its children is the active
  clash.

Row click → `showClashAndFocus(plugin, { a, b, severity })` in
`molstar-helpers.ts`:
- Renders both residues as sticks via `showSticksForLoci` under tag
  `CLASH_TAG = 'tarantino-clash'`.
- Calls `showClashLine(plugin, a, b, severity)` which uses
  `plugin.managers.structure.measurement.addDistance(lociA, lociB, …)`
  to draw a **dashed line + distance label** between the two specific
  clashing atoms — **amber** for `bad` (`0xe68a00`), **red** for
  `severe` (`0xc62828`). The returned `selection` + `representation`
  cell refs are tracked in a per-plugin `WeakMap<PluginUIContext,
  string[]>` so subsequent clicks (or `clearClashSticks`) delete them
  precisely without touching any other measurements the user might
  add through Mol*'s own UI.
- Camera focuses the residue-pair bounding sphere.

Atom-level loci is built by a new synchronous helper `buildAtomLoci`
(MolScript with `chain-test` + `residue-test` + `atom-test` on
`label_asym_id` / `label_seq_id` / `label_atom_id`), mirroring the
existing `buildResiduesLoci` pattern.

`clearClashSticks` deletes both the residue-stick cells (by
`CLASH_TAG`) and the tracked dashed-line cells. Wired into
`useMolstarSync.clearPlugin3DState` (empty-3D-click) and `App.tsx` Esc
handler.

### DVBFixer Panel + Backend

**Frontend** (`src/components/DVBFixerPanel.tsx`): MUI Tabs (one per
sub-command), input file picker sourced from active workspace artifacts, and a
flag form generated from the spec fetched at runtime from `GET /api/dvbfixer-spec`.
Per-flag controls map by `type`: bool → Checkbox, select → Select,
number/text → TextField.

**Auto-paste input** — the input picker tracks `structureStore.fileName`
(currently-loaded primary structure). Whenever the user hasn't manually picked
an input yet (`userPickedInputRef.current === false`) OR the current selection
is empty, the picker mirrors the active structure. Selecting from the dropdown
sets `userPickedInputRef.current = true` and stops the auto-mirror.

**Managed jobs and output** — Run posts to `/api/v1/workspaces/{workspaceId}/jobs`. The panel restores an
active workspace job after remount, follows state through EventSource with a
polling fallback, exposes Cancel, and shows terminal logs. On success it reloads
workspace artifacts, loads structure output into viewer A, and resets
`userPickedInputRef` so the next command uses that output. Failures leave the
viewer untouched.

**Model tab — per-chain FASTA input.** The `model` sub-tab renders a
custom section above the flag controls: one multi-line input per
polypeptide chain of the loaded primary structure (filtered via
`filterSequenceableChains`). A **Parse from PDB** button populates every
box via `chainToSequence(chain.residues)` from
`src/lib/alignment.ts` (SEQRES-aware — residues missing from ATOM coords
are still included). A **Clear** button empties all boxes. The standard
`--fasta` text field is hidden in this tab because the per-chain UI
synthesises it automatically.

**Inline missing-residue highlighting** — each per-chain box is a
`HighlightedFastaInput` (defined at the top of `DVBFixerPanel.tsx`): a
transparent-text `<textarea>` layered over a colored underlay `<div>`.
Both layers share identical font / padding / wrap rules so they
align pixel-for-pixel. After Parse from PDB, SEQRES-only residues
render greyed `#b5bfcc` + italic + fontWeight 400 (same visual as the
Sequence panel's missing residues). Highlighting is gated on
`value === parsedSequenceByChain[c.id]`: edit any character and the
whole box reverts to plain text so user-typed characters never get
mis-classified as missing. `parsedSequenceByChain` and
`presentMapByChain` snapshots are stored alongside `fastaByChain` and
wiped by Clear / input-switch.

On Run, `buildFastaContent()` assembles a valid FASTA string from the
non-empty chain boxes (60-char-wrapped lines, `>{inputBase}_{chainId}`
headers) and ships it as `fastaContent` in the request body. The
backend (the versioned managed-job route in `server/managed-jobs.ts`)
writes the content to `<outDir>/<inputBase>.fasta` and injects
`--fasta <abspath>` into the CLI args — overriding any user-typed
`--fasta` value. The materialised FASTA stays beside the output PDB so
the user can inspect / reuse it.

The per-chain UI is only ENABLED when the picker's input matches the
primary viewer's `fileName` (we need the chain list). Otherwise an
Alert tells the user to load the structure first; they can still leave
the boxes empty and let DVBFixer fall back to SEQRES from the input
PDB.

**Backend** (`server/api-routes.ts`, shared Node route composition):
- `GET /api/dvbfixer-spec` — returns `COMMANDS` from `server/dvbfixer-spec.ts`.
- `POST /api/v1/workspaces/{workspaceId}/jobs` starts a workspace-scoped DVBFixer job. `GET` on the same URL
  restores active jobs after a panel reload, job detail/events report status
  and logs, and `DELETE /api/v1/workspaces/{workspaceId}/jobs/{jobId}` requests cancellation. Successful output
  files are registered as workspace artifacts; the frontend reloads the active
  manifest and opens the primary output when appropriate.
- `PATCH /api/workspaces/:workspaceId/artifacts/:artifactId/metadata` updates
  whitelisted artifact metadata against the current manifest revision and
  returns the full incremented manifest. Optional fields accept `null` for
  removal; stale revisions return 409.

**Spec format** (`server/dvbfixer-spec.ts`):
- `FlagDef.type`: `'bool' | 'number' | 'text' | 'select'`
- `FlagDef.repeatable: true` — comma-split UI input becomes `--flag v1 --flag v2 --flag v3` (used by `--mutate`).
- `FlagDef.multi: true` — value is whitespace-split and emitted as a single `--flag v1 v2 v3` (argparse `nargs='+'`). Works with both `type: 'text'` and `type: 'select'`; the latter lets a dropdown preset like `"amber19/protein.ff19SB.xml amber19/tip3p.xml"` resolve to the right multi-arg CLI form. Used by `--ff` (minimize + protonate).
- Empty string in any `select`'s `options` is preserved as the "default" choice and the backend drops empty values before arg-building, so DVBFixer's built-in defaults apply.

### Settings (`src/components/SettingsPanel.tsx`)

App-wide preferences live in the structure store and are persisted in
`localStorage`. Two preferences today:

- **Viewer → Auto-orient on load** (`autoOrientOnLoad`, default **OFF**,
  key `tarantino.autoOrientOnLoad`). When ON, every loaded structure
  goes through `PluginCommands.Camera.OrientAxes` in `MolstarViewer`
  post-load (camera sync suppressed for the duration of the orient).
  When OFF, the structure keeps its authored orientation but is still
  fit to the viewport via `camera.reset(undefined, 0)` (necessary
  because `manualReset: true` suppresses Mol*'s built-in auto-fit).
- **Alignment → Source label** (`alignmentLabelMode`, default
  `'file'`, key `tarantino.alignmentLabelMode`). Toggles what the
  AlignmentPanel's source-labels block (above the alignment view)
  shows for each side: the **file** path (default — e.g.
  `FcRn.pdb`) or the entry's metadata **name** (the user-editable
  `name` field from the Info panel). When `'name'` and a non-empty
  name exists, falls back to the file path; otherwise displays the
  file path.

The panel renders one MUI control per preference (`Switch` for booleans,
`ToggleButtonGroup` for enums). Setters write through to `localStorage`
immediately via the matching `persistX` / `loadPersistedX` helpers at
the top of `structureStore.ts`.

`AlignmentPanel` builds its file-to-name map directly from the active
workspace's artifacts. Info edits update the same revisioned manifest, so both
primary (A) and secondary (B) source labels react without a legacy library
refresh signal.

Discoverable in the default layout's left column (paired with Info) as
of the latest layout update; also accessible via the `+` menu.

### Empty 3D Click Behavior

`useMolstarSync.attachEmptyClickCleanup` is a single unified handler attached
to both viewers' `interaction.click`. Default Mol* `clickDeselectAllOnEmpty`
requires `selectionMode === true` (default false), so Mol* doesn't auto-clear;
we check `Loci.isEmpty(event.current.loci)` ourselves (NOT `isEmptyLoci(loci)`
— Mol* fires a `StructureElement.Loci` with empty `elements`, not the
`EmptyLoci` singleton).

## Key Constraints

- **Mol* imports use deep paths** (`molstar/lib/mol-model/structure`) — no barrel export.
- **TypeScript strict**: `noUnusedLocals`, `noUnusedParameters`, `verbatimModuleSyntax`, `erasableSyntaxOnly` all on.
- **`@` alias** → `src/` (in `vite.config.ts`).
- **Workspace manifests are authoritative** for visible artifacts; hidden
  bookkeeping files must not enter file pickers or the Workspace panel.
- **The frontend never imports from `server/`** — it talks to the backend over HTTP. The DVBFixer spec is fetched at runtime from `/api/dvbfixer-spec`.
- **postinstall** `scripts/fix-native-deps.mjs` installs the right platform-specific `@rollup/rollup-*` binding.

## Workspace and Homology architecture

The authoritative user guide is `../docs/gui-homology.md`. The Library panel
is `ProjectLibrary.tsx`, backed by `workspaceStore.ts` and
`server/workspace-api.ts`. Workspace imports, renames, reorders, and
recoverable trash operations must flush/cancel pending manifest autosaves and
reload both the active manifest and workspace summaries. Helpers and run files
remain on disk with `hidden: true`. `TextFileViewer.tsx` owns read-only preview
of text-like artifacts.

`HomologyPanel.tsx` edits workflow state only: target FASTA, template/target
chain assignments, rectangular MSA rows, and masks. Build requests are handled
by `server/homology-api.ts`, which writes `template-plan.json` with absolute
source paths and invokes `dvbfixer homology --template-plan`. Scientific
fitting/mosaic/PIR logic must stay in Python (`dvbfixer.homology_plan`), not be
reimplemented in React or the Vite server.

Selection rules are click=replace, drag=range, Shift=additive extension from
the newest anchor, and Cmd/Ctrl/Option=toggle plus new anchor. Mask ranges are
zero-based and half-open. Earlier template rows win overlapping columns. The
build emits one coordinate-preserving mosaic known, never independent knowns
per painted span. Multi-chain plans share the first template structure as a
global reference; `VH` and `VL` map to PDB `H` and `L`.

## File Map

```
src/
  App.tsx                       # Layout, panel factory, "+" menu, camera sync icon, Esc handler
  main.tsx                      # MUI ThemeProvider, CssBaseline, ErrorBoundary
  theme.ts                      # MUI createTheme
  index.css, molstar-theme.scss # FlexLayout vars + Mol* SCSS skin

  components/
    MolstarViewer.tsx           # Mol* init per slot, color theme registration, post-load (hide water, ions→spacefill, OrientAxes gated on autoOrientOnLoad). Cleanup disposes the plugin and clears that slot's filename/chains.
    SequenceViewer.tsx          # Monospace residue grid, drag-select, missing-SEQRES gap rendering, validated chain init, Structure/Sequence numbering toggle
    ProjectLibrary.tsx          # Library workspace list + active Workspace files, drag handles, filtering, A/B load selector, downloads, rename/trash
    StructureInfo.tsx           # Stats and revision-aware workspace-artifact metadata, including equivalent-chain override
    ElementsTable.tsx           # Tree, visibility toggles, row-click camera focus (sync-suppressed), "Show Interface"
    InteractionsPanel.tsx       # Computed contacts, focused-chain banner
    ClashesPanel.tsx            # VdW-overlap clash table: severity filter, Group-by-residue toggle (expandable groups), row-click → residue sticks + severity-colored dashed clash line via measurement.addDistance
    AlignmentPanel.tsx          # Pairwise NW alignment, per-source plugin routing
    DVBFixerPanel.tsx           # MUI Tabs, form from /api/dvbfixer-spec, auto-pastes active fileName, auto-loads output on success
    SettingsPanel.tsx           # App-wide preferences (auto-orient-on-load, Alignment source-label mode); all localStorage-persisted
    FileLoader.tsx              # Upload button (honors loadTargetSlot)
    ProjectLibrary.tsx          # Workspace switcher/file browser: import, rename, recoverable trash, zebra rows, Shift/Cmd selection, reorder
    TextFileViewer.tsx          # Read-only text-artifact preview tab
    ChainSelector.tsx           # Chain dropdown (same filter+sort as SequenceViewer)

  hooks/
    useMolstarSync.ts           # 3D → sequence + empty-click cleanup on both viewers
    useSequenceSync.ts          # Sequence → 3D (cartoon halo + solid sticks)
    useCameraSync.ts            # Bidirectional camera mirror via canvas3d.didDraw

  stores/
    structureStore.ts           # plugin, secondaryPlugin, loadTargetSlot, chains+secondaryChains (each residue has seqId + optional authSeqId), elements, meta (incl. iggSubtype + allotype + optional equivalentChains override), focusedChainId+Category, cameraSyncEnabled, autoOrientOnLoad + alignmentLabelMode (both localStorage-persisted), clearAllSignal
    selectionStore.ts           # selected/hovered residues, _lock mechanism
    workspaceStore.ts           # workspace summaries, active manifest, debounced persistence, text-preview target

  lib/
    molstar-helpers.ts          # MolScript builders, showSticksForLoci, extractChains (label_seq_id as seqId + auth_seq_id as authSeqId + SEQRES merge + present flag), buildAtomLoci + showClashLine (severity-colored dashed line via measurement.addDistance, per-plugin ref tracking)
    clash-detection.ts          # computeClashes: VdW-overlap pairs via structure.lookup3d.find + Bondi/R&T radii; gatherNeighbors walks both unit.bonds + structure.interUnitBonds so 1-2 / 1-3 exclusions cover glycosidic + interchain disulfide bonds; severity tiers bad (0.4–0.9 Å) / severe (>0.9 Å)
    alignment.ts                # Needleman-Wunsch + BLOSUM62 + chainToSequence + trimmedIdentity (terminal-gap-aware)
    chain-grouping.ts           # computeEquivalentChains (pairwise NW + union-find), validateGrouping, filterSequenceableChains
    residue-codes.ts            # 3-to-1 letter code (incl. non-canonical)
    residue-color-theme.ts      # Mol* ColorTheme: carbons by residue class, others CPK

server/
  api-routes.ts                 # Host-neutral composition: workspaces, Homology, managed jobs, naming, security, and observability
  api-plugin.ts                 # Thin Vite development adapter over api-routes
  standalone.ts                 # Loopback-default Node HTTP + static client host
  workspace-api.ts              # Revisioned workspace/artifact CRUD, contained file serving, recoverable trash, and legacy-index migrations
  homology-api.ts               # Homology project persistence; writes CLI template plans and registers run artifacts
  dvbfixer-spec.ts              # CommandDef[] for split/renumber/model/prepare/minimize/protonate/convert (was `glycam` in older DVBFixer). renumber.--scheme options: seqres/kabat/chothia/imgt/martin/eu/aho. convert exposes --to-amber + --to-charmm + --no-roh. minimize + protonate `--ff` is a select-multi dropdown of OpenMM bundles (AMBER19/AMBER14/GLYCAM/CHARMM36 presets, empty = DVBFixer auto-pick). minimize defaults `--no-solvent` to ON for interactive development runs. prepare exposes --no-infer-conect. model exposes --num-output (top-N candidate save count; with N>1 the auxiliary _2.pdb/_3.pdb outputs stay on disk and Tarantino only auto-loads the first). Removed from UI: model --keep-workdir, minimize --dat/--padding/--platform, protonate --cys-disulfide-pka, protonate --protassign (DVBFixer now defaults it ON). Hidden flags are still valid on the CLI.

scripts/
  fix-native-deps.mjs           # postinstall
  set-modeller-key.sh           # Writes Modeller license into config.py

structures/projects/<id>/       # workspace.json plus files/, runs/, homology/, recoverable .trash/, and legacy migration recovery sources
vite.config.ts                  # React + API plugins; workspace files are served only through the contained workspace file route
```
