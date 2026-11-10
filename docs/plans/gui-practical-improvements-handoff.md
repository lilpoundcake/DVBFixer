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
    accessible cascading group menus. User approved **Glycans** on 2026-11-06;
    this names the subject rather than limiting the category to preparation.
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

- Job Monitor now has a version-1 stage/progress contract and a dedicated bounded
  NDJSON pipe from the managed Python process. Latest stage persists in the job
  record and is available through existing GET/list/SSE transport after reload.
  The diffusion pipeline reports preparation, combined runner execution,
  independent validation and publication at the actual operation boundaries.
  Its child boundary forwards only the validated host pipe. Maintained cc94,
  corrected/historical cc89 and Protenix wrappers now emit separate sampling and
  refinement stages; per-seed wrappers include the active seed. No sampling
  settings, process layout, checkpoint identities or request schemas changed.
  Counts, when supplied, are stage-local; current pipeline events are deliberately
  indeterminate. Terminal jobs show the last reported stage, not a readiness claim.
  Checkpoint-backed acceptance of newly instrumented wrappers and measurable update
  counts remain open. Homology modeling now participates through an additive
  persistent V1 endpoint after explicit user approval: same locks, process admission,
  Stop, persisted stages and monitor records; no new scheduler. The old synchronous
  endpoint remains compatible and request-bound. Existing installed runtime copies
  are not silently updated.
- A native Chrome controlled-CLI smoke now covers two concurrent workspace jobs,
  independent cancellation, stage display, persisted mixed-seed outcomes after
  reload/resize and retained report navigation. It does not close real CPU/GPU
  resource admission or full multi-principal browser acceptance. Source tests now
  cover ACL revocation closing SSE, queued cancellation, stale revisions,
  Retry-After, Stop during packaging and shutdown draining. Native keyboard,
  pointer-transition and narrow touch cascading-menu acceptance also passed.
  Keep the overall scientific/hardware checklist open.
- Exercise a real multi-seed GUI sampling workflow and per-seed result display,
  especially mixed passing/rejected candidates. New bundles retain a separate
  digest-linked seed inventory with final independent summaries and rejected
  inspection coordinates. Only accepted candidates remain selectable; schema-4
  results and old bundle compatibility remain unchanged. Model/ZBS, Homology and
  Job Monitor render the shared report without Node-owned geometry calculations.
- Broaden reporting coverage for execution failures before any bundle exists and
  all scientific failure categories. Mixed-result final-only reporting and GUI
  retained-path relocation now have Python, process and browser controls; native
  scientific GUI sampling remains a separate gate.
- User approved **Glycans** for Convert → Transplant → Cluster. Registry and
  generated views use that label; saved `Glycoprotein preparation` panel state
  restores to Glycans without changing command IDs or per-tool values. Verification:
  461 GUI tests and 11 command-registry tests passed, along with GUI typecheck,
  Ruff, OpenAPI consistency, maintained user-document checks and diff whitespace.

## Verification so far

- 2026-11-06 persistent Homology continuation: **481 GUI tests**, GUI typecheck/lint,
  generated OpenAPI consistency, production client/server build and standalone
  smoke passed. The final Python fast rerun passed **2202 tests**, with the same
  three optional Torch/Biotite skips and 33 slow deselections. Native Chrome
  evidence at `tmp/job-monitor-browser-yp4IDY/` includes GUI-created Homology jobs
  surviving navigation away and reload, restoration/Stop/completion, group migration,
  keyboard/pointer/touch navigation and retained values. These use controlled CLI
  outputs, not pretrained weights. Actual licensed MODELLER API smokes at
  `tmp/managed-homology-native-37x70htf/` and the final built-server-bound
  `tmp/managed-homology-native-final-4U5XKd/` passed independent output readback:
  chain A, 76 residues, 602 atoms, valid DatRecord, zero D-Cα. This is a full-template
  1UBQ workflow smoke, not gap-reconstruction or broad scientific qualification.
  Source/template snapshots and output provenance are retained; installed runtimes,
  model weights and scientific thresholds were not modified.

- 2026-11-06 continuation: 460 GUI tests, typecheck/lint, production client/server
  build and standalone smoke passed. The native Chrome/controlled-CLI evidence is
  retained under `tmp/job-monitor-browser-uAlsh2/`; reproduce after building with
  `node gui/scripts/smoke-job-monitor.mjs` (set `CHROME_EXECUTABLE` off macOS).
  This is transport/UI acceptance, not checkpoint sampling. Targeted Python
  validation/workflow/publication controls pass, including unchanged ranking,
  mixed outcomes, Unicode/case/insertion-code findings and corruption rejection.
  The first full fast run exposed five old-bundle compatibility failures because
  optional reporting accidentally required a new provenance key for old bundles.
  The reader was corrected without rewriting fixtures or weakening new-report
  bindings; all 63 focused validation/workflow/publication tests then passed.
  The final complete fast rerun passed **2201 tests**, with three existing optional
  Torch/Biotite skips and 33 slow cases deselected. Both failed and passing logs are
  retained in `tmp/job-monitor-verification/`; no fixtures or expected scientific
  hashes were rewritten to obtain that result. A further native Chrome pass is
  retained under `tmp/job-monitor-browser-gPDiq8/`, with the final built-source pass
  under `tmp/job-monitor-browser-mBrFI9/`. The isolated sampling Python has
  no pytest module; no runtime mutation was attempted to close optional controls.
  Reserved independent-finding prefixes in adapter warning text are now explicitly
  labelled as runner warnings, preventing them from impersonating Python-measured
  final geometry. Mocked cc94 and Protenix seed-loop tests verify stage ordering
  and unchanged sampling options; they are not checkpoint-backed acceptance.

- 2026-11-06 stage-transport slice: 449 GUI tests, GUI typecheck/lint, 40 focused
  Python progress/pipeline tests, Ruff and generated OpenAPI consistency passed.
  Process tests use a real child producer on the separate pipe and verify persisted
  API readback; malformed/oversized/version-mismatched events are ignored, work
  counts are validated and telemetry delivery is bounded. Python tests cover a
  closed monitor, unchanged stdout/stderr, real pipeline stage order and admission
  refusal without later-stage events. No checkpoint sampling or live browser
  multi-job acceptance was executed. Bootstrap Python lacked OpenMM; pipeline
  tests were instead run against source in the existing scientific environment.
  An actual Python emitter → Node child-pipe smoke also passed with unchanged
  stdout/stderr; retained evidence is in
  `tmp/job-progress-smoke-gnJqJH/summary.json`. This tests transport only, not a
  sampler/refiner or scientific result.

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
