# GUI practical improvements handoff

Date: 2026-11-03. Branch: `feature/gui-practical-improvements`, created from
local `main`. Release: 0.10.3. Status: initial implementation delivered with the
remaining gates below; no scientific acceptance change.

## Scope and order

1. Actionable English diffusion rejection output: selected-candidate findings,
   actual measurements/limits and retained report paths, visible without verbose.
2. GUI multi-seed entry: commas/whitespace, strict validation, repeated CLI flags,
   persistence and per-candidate outcomes.
3. Navigation: Altloc before Model; Pull in Utilities; Diagnose first in Utilities;
   Utilities keeps Conect → Puppet → Pull adjacent (user follow-up);
   the glycan group keeps Convert → Transplant → Cluster;
   Cluster in the former Glycoprotein Preparation group; remove Analysis;
   accessible cascading group menus. The broader destination name still needs
   user confirmation before renaming.
4. Job Monitor: a separate tab next to Info in the same docking tabset,
   independent managed-job records, stages, honest progress, logs
   and per-job cancellation. Reuse existing API/locks; do not promise parallel
   same-workspace execution or remove resource/ownership safeguards.

The [canonical checklist](implementation-checklist.md) contains detailed
acceptance requirements. Keep unfinished items unchecked.

## Boundaries

- Python owns scientific findings and validation. Node/React only transport and
  display them. Do not change geometry thresholds or silently accept rejections.
- Keep schema-4 frozen requests and old bundle compatibility intact.
- Keep final findings separate from resolved raw/stage-1 findings; never trust
  unchecked runner report text as independent scientific acceptance evidence.
- Preserve fd-level output capture, compact retained run storage, cancellation
  evidence, case-sensitive chains and insertion codes.
- Shared layout rules are in [GUI design system](../gui-design-system.md).
- Existing untracked `FF/` directories are unrelated user work; leave untouched.

## Starting points

- `src/dvbfixer/model/diffusion_cli.py`: rejection output currently prints gates.
- `src/dvbfixer/model/diffusion/validate.py`: independent final geometry findings.
- `src/dvbfixer/model/diffusion/staged_runtime.py`: independently audited anchors.
- `gui/src/components/DVBFixerPanel.tsx`: generated tool forms and group selector.
- `gui/src/lib/managed-jobs.ts`, `gui/server/managed-jobs.ts`: shared job transport.
- `gui/src/app-panels.ts`, `gui/src/App.tsx`: dockable panel registration.

## Verification and continuation

Run narrow Python/GUI tests first, then GUI typecheck and full GUI tests.
Regenerate GUI spec only through `scripts/gen_gui_spec.py` when its source
changes. Run `git diff --check`. Scientific integration requires a healthy
environment; record blockers rather than broad skips. Any development minimize
or ZBS invocation must use `--no-solvent`. Record completed slices and remaining
gates below; no expensive sampling is required for initial reporting/UI tests.

## Progress

- Branch and scoped handoff created; earlier approved checklist changes retained.
- Independent validation now retains final residue-level geometry findings in
  summary warnings; independently audited anchor-limit violations include the
  measurement and limit. CLI renders these for rejected published candidates
  without verbose, with candidate/seed identity, gate fallback and report path.
  Rejection causes are now emitted as ERROR rather than WARNING. ToolOutput
  and Homology show explicit error/validation failures ahead of warning banners;
  report/log error lines are bold/error-colored without coloring all stderr.
  Repeated fd-summary copies are deduplicated in the banner. Historical explicit
  rejection warnings are also surfaced, without inferring scientific findings.
- Model/ZBS and Homology accept comma/whitespace-separated seeds; shared strict
  normalization feeds repeated CLI flags. Homology's former single-number API
  restriction is removed. Numeric old panel values remain compatible.
- Shared registry/spec navigation has Altloc before Model, Diagnose first in
  Utilities, Conect → Puppet → Pull and Convert → Transplant → Cluster. Analysis
  is removed and restored Analysis panel state falls back to Utilities.
- A cascading group menu is implemented. User-reviewed styling uses a compact
  outlined dropdown with a visible arrow (the trial primary fill was reverted).
  Submenus prefer right, then below/above rather than flipping left; resizing
  repositions menus and avoids a full-width trigger. The group dropdown sizes
  to its longest group label without wrapping (viewport-limited), independently
  of the compact trigger width.
- Job Monitor is registered immediately after Info in the same default docking
  tabset. It polls persisted jobs across accessible workspaces, shows lifecycle
  status/elapsed time, indeterminate progress, independent Stop and log links.
  Its description is in a separate header above a scrolling Tasks section.
  Task rows reuse shared Workspace/inspection stripe and hover colors.
- Fixed a rate-limit regression in the initial monitor: all-workspace polling
  every two seconds could exceed the API's 120/minute limit. Polling is now
  one ACL-filtered `/api/v1/jobs` request for all active jobs plus the latest
  30 completions, refreshed at most 12 times/minute. Older completions load only
  through **Load older tasks**, using a stable time/identity cursor. Active jobs
  remain visible regardless of history position; cached history is pruned for
  workspaces no longer readable. This supersedes the interim slow per-workspace
  round-robin mitigation. Polling pauses for hidden tabs and honors Retry-After. Transient
  errors retain the last known jobs. Workspace-list refreshes do not restart
  polling bursts; opening Job Monitor again focuses the existing tab. The tool
  panel no longer polls redundantly alongside a healthy SSE stream, and its
  fallback also honors Retry-After. Server rate-limit protection is unchanged.
- Workspace run/file presentation now uses independent parent-run bands and
  per-run file stripes, reset inside each run and stable across expansion.
  User review reduced the initial tint strengths, then split parent bands into
  soft blue/lavender families with white light child stripes. Further visual
  review added inter-run spacing, bold headers and stronger colored group edges
  rather than relying on barely distinguishable tints; shared selection/hover
  overrides remain unchanged. Palette and hierarchy rules are centralized in
  `ui/layout.ts` and the GUI design-system document.

## Remaining implementation/acceptance gates

- Job Monitor still needs the versioned stage/progress event contract and actual
  runner instrumentation. Current status is lifecycle-only; there are no fake
  sampling/refinement percentages. Homology synchronous jobs are not yet part of
  the managed-job list. Keep the overall Job Monitor checklist open.
- Exercise real multi-job monitor/cancellation and resource admission, reconnect,
  role changes, keyboard/touch cascading-menu behavior and narrow viewports in
  the browser. Static-render/contract tests do not close these acceptance gates.
- Exercise a real multi-seed GUI sampling workflow and per-seed result display,
  especially mixed passing/rejected candidates. Current validated-result logic
  can omit rejected seeds from a passing subset; do not claim full per-seed
  reporting from the input/argv fix alone.
- Finish reporting coverage for all failure categories, mixed candidate results
  and GUI retained-path relocation. The new CLI report path is correct at CLI
  publication time but GUI staging may relocate that bundle afterward.
- Confirm the new glycan-group name with the user; current grouping retains
  `Glycoprotein preparation` until that decision.

## Verification so far

- Python: 76 focused tests passed across validation, independent anchor audit,
  staged integration and command registry using the existing `dvbfixer` env.
- GUI: 442 tests passed across 64 files after compact-dropdown/header/row-color
  polish; typecheck passed.
- Text Files now wraps long lines/paths to the current docking-panel width,
  preserving whitespace and monospace rendering without stretching the panel.
- Latest rejection-output follow-up: all 26 validation tests passed again.
- Final pre-commit: all 76 focused Python tests re-passed, plus 5 environment/
  version-policy tests. Ruff, agent docs, OpenAPI synchronization and the three
  release metadata/CLI/GUI generated checks passed. Version metadata, current
  installation example and changelog are synchronized at 0.10.3.
- GUI spec regenerated from registry; CLI reference index regenerated after
  command-order changes. Both generated consistency checks and `git diff
  --check` pass. No diffusion sampling or full scientific fast suite executed here.
- The isolated Standard prefix has no pytest module; no package installation or
  prefix repair was attempted. Existing `dvbfixer` provided narrow-test tools;
  these tests do not establish its PDBFixer resource readiness.
