# Implementation Checklist

## Installer workflow and progress follow-up (2026-11-05)

### Explicit profile attachment and numbered menus

- [x] Offer numbered actions, profiles, managers, backends, confirmations and add-ons.
- [x] Add `--extend` with separate Standard, full Desktop and `--gui-only` scopes;
  existing non-base conda environments do not need a helper ownership record.
- [x] Preview native and pip changes before confirmation, detect native/environment
  drift, recheck pip after native installation, and stop on unconfirmed pip changes.
- [x] Preserve existing Desktop/add-ons, license and workspaces. GUI-only rejects
  changes to pre-existing native libraries and keeps scientific profile ownership.
- [x] Publish full-profile ownership only after functional/Doctor/GUI checks pass
  (or explicitly report required license action), not before dependency installation.
- [x] Finish regression/documentation validation and record platform evidence below.

Native and pip installation are not an atomic transaction. Python minor/ABI
changes and conda Open Babel overlays are refused for in-place full-profile
attachment; use a new environment for incompatible stacks. A solver preview is
not scientific acceptance. Cancelled/refused previews retain diagnostic logs in
the checkout's ignored `tmp/installer-profile-preview-*/` directories, without
changing the selected environment.

Validation for this slice: 152 targeted installer/composition/documentation tests
passed, including numbered-only install/backend/extension flows, simulated
transactions through all three manager adapters, native/pip conflict/drift
refusals, interrupted GUI build, license/workspace/.env preservation, and ownership
publication ordering. Ruff, documentation/link checks, release metadata, inventory
freshness and agent-map checks passed. A separate real bootstrap/PTY process
completed the all-numeric GUI-only preview against a **simulated manager**, with
byte-identical target fixture and no ownership record; retained local evidence is
at the approved temporary `installer-numbered-smoke-g5q_cku3/` directory. This is
not a real native package solve or installed scientific/GUI acceptance run. Existing
user environments were not modified. SciPy/OpenMM-dependent full-suite and
generated-reference gates remain blocked in the available bootstrap Python as
noted below; no scientific references changed.

### Recovering a full-profile conflict without replacing the scientific stack

The reported conda Open Babel failure occurs before native/pip installation.
Keep the no-overlay invariant, but let an interactive user explicitly switch to
GUI-only, with a fresh plan and confirmation, or stop without changes. Unattended
full-profile requests still fail rather than silently narrowing scope. Requested
new scientific add-ons/license changes are skipped only after the explicit GUI
choice; an opted-in application shortcut is retained. GUI-only preserves known
editable package mode. Default extension selection shows DVBFixer installations
first and puts other environments behind a numbered entry. All discoverable
manager bases and macOS application bundles are excluded, including a foreign
Miniforge base registered with micromamba.

Validation: 161 targeted installer/composition/documentation tests passed, along
with Ruff, maintained Markdown/link checks and agent maps. Controlled fixtures
reproduce the reported Open Babel stack and cover GUI-only recovery, explicit
cancellation, retained shortcut opt-in, known editable mode, unattended refusal
without scope changes, default-list grouping and exclusion of foreign manager
bases/application bundles. Existing user environments were not modified. This is
not native installation or scientific qualification evidence; the full-profile
no-overlay guard remains enforced.

### Explicit update scope and diffusion transition

Follow-up UX fixes: the diffusion transition reuses the chosen manager, with an
explicit change-manager option. Runtime selection shows a distinct heading, full
path and retained provisioning state. Model-asset guidance lists exact cc94 files
and Zenodo before selection and in missing-asset results. Required-action summaries
do not label completed version snapshots as failures or infer network problems
from their package inventories.

After user approval, `install_weights.py` adds consented, repository-local shared
`models/` storage (Git-ignored). cc94 uses the official pinned Zenodo archive and
retains only its weights/config; Protenix uses the official v1 checkpoint URL.
Both verify existing frozen SHA256 identities before publication/reuse. Downloads
show size/terms/source, support local-file and deferred alternatives, and do not
run in previews. Automated download needs `--download-weights --accept-model-terms`.
Partial transfers are private staging, not cache hits; corrupt existing caches
are preserved/fail closed. Protenix data/Kalign remain separate manual inputs.

Launcher conflicts now identify the publication stage and provide specific
recovery. Interactive replacement is offered only after preflight and discloses
shared per-user impact; declining leaves the old launcher and action-required.
Unattended replacement remains explicit. Confirmation-time drift is rejected.
Tests use small synthetic transfers/archives and mocked provisioning, not real
model downloads, user-environment mutations or scientific qualification.

Validation: 207 targeted installer/runtime/composition/user-doc tests passed,
along with Ruff, maintained documentation links/anchors, agent maps and diff
whitespace checks. Cases include download consent and dry-run/no-network paths,
official-profile cache sharing, checksums, interrupted transfers, corrupt/linked
caches, exact-member extraction, plaintext redirect refusal, explicit manager
changes, runtime path/state presentation, completed-version-log exclusion from
required-action diagnosis, launcher replacement/decline and confirmation drift.
The generated CLI-reference check remains blocked by missing OpenMM in the
bootstrap Python; no generated command-reference pages were edited by hand.

Manual update selection now offers **only package / full installation / cancel**,
instead of conflating a negative answer to "only package?" with cancellation.
Full selection routes to the existing profile preview and keeps the final apply
confirmation separate. It does not update arbitrary unrelated packages or replace
frozen checkpoints. An explicit `--package-only` still bypasses the scope chooser.

The selected-application component menu now includes **Diffusion backend**. It
selects a host-platform maintained runtime and offers creation by name or recorded
runtime continuation by number, reusing existing provisioning checks. Runtime and
scientific prefixes cannot be identical, no Torch is installed into the selected
application, and per-user launcher discovery is explicitly disclosed rather than
claimed as exclusive per-prefix binding. Scientific add-on/license/shortcut flags
cannot silently cross into this separate operation.

Validation: 164 targeted installer/composition/documentation tests passed, plus
Ruff and maintained user documentation checks. New controlled cases exercise the
second update option's full-profile preview, isolated MPS provisioning transition,
existing runtime continuation, and byte-identical selected application files.
The provisioning calls in these cases are mocked; no user environment, checkpoint
or launcher was changed and no native sampling/qualification is inferred.

- [x] Separate main-application install/update from optional backend selection.
- [x] Add explicit helper-owned Standard/Desktop update with composition binding,
  version preview and interactive confirmation; refuse manual/runtime adoption.
- [x] Infer resume identity from the existing installation record and distinguish
  filesystem directories from environment names.
- [x] Remove installation-path prompts: new installs use an environment name in
  manager-reported writable directories; update/resume use numbered discovered
  installations. Keep explicit `--prefix` advanced-only and add `--name` automation.
- [x] Read installed versions from metadata without executing discovered Pythons;
  distinguish duplicate names, retain crash-left reservations and refuse manual
  adoption or ambiguous name-only selection.
- [x] Add named stages, live output, periodic elapsed-time notices, retained logs
  and human-readable readiness/action-required summaries.
- [x] Preserve protected MODELLER configuration across successful updates and keep
  license subprocess output out of logs.
- [x] Complete regression and documentation validation; record evidence below.

Scope does not include transactional package rollback, automatically adopting
manual environments, changing frozen backend recipes or scientific qualification.

Verification: 667 CI-selected Python contract tests passed, including the new
workflow tests in the lightweight CI lane. Ruff, strict API mypy, release metadata,
license-inventory freshness and 43-page user-document checks passed. Real
micromamba 2.9.0 name-based Standard installation and numbered existing MPS-backend
resume PTY dry runs passed without installing packages or creating prefixes or
ownership records. A prior explicit-command Desktop dry run also passed. Process tests
exercise actual child stdout/stderr retention and failure reporting; controlled
update tests cover composition binding, cancellation, license recovery after
failure/interruption, Desktop rebuilding and workspace/.env preservation.
No existing user prefix was updated, no fresh scientific install was executed,
and no additional backend was sampled or newly qualified by this follow-up.

### Existing/manual installation correction

The first selector excluded unrecorded environments, leaving existing users with
no update choice. It also inspected only wheel metadata. Update discovery now
includes legacy egg-info, editable egg-link/PEP 660 evidence and retained CLI
entry points with unknown versions. Manual entries are selectable for explicitly
confirmed **package-only** updates (`pip --no-deps`), never silent profile adoption.
Native dependencies, MODELLER configuration and workspace files are not modified;
known editable mode is retained. Automated manual updates require `--package-only`.
Frozen runtime roots/children and corrupt ownership records cannot use this route.

Latest validation: 125 focused installer/composition/documentation tests passed,
plus Ruff, Markdown/link checks, agent maps, release metadata and inventory
freshness. Tests reproduce the missing-metadata CLI, legacy/editable layouts,
manual menu selection, cancellation and dependency/license/workspace preservation.
Direct inspection of the older user prefix was permission-blocked; no alternate
access or user-environment update was attempted. The available bootstrap Python
lacks SciPy (one environment-policy smoke failed) and OpenMM (CLI-reference check
was blocked). These are retained environment blockers, not passing scientific or
generated-reference acceptance; no scientific references were changed. The older
667-test evidence above predates this correction and is not a new full-suite run.

This is the canonical completion tracker for the plans in this directory.
Current cross-plan navigation is in [Plan status and remaining work](README.md).
Status synchronized on 2026-10-20: software delivery, per-profile hardware
acceptance and wider scientific/default promotion are separate completion gates.
Check an item only after its implementation and required verification are
complete. Detailed requirements remain in these authoritative sources:

- [API roadmap](api-roadmap.md) for API implementation and deployment.
- [DDD documentation roadmap](ddd-documentation-roadmap.md) for agent maps and domain documentation.
- [Diffusion gap-reconstruction plan](diffusion-gap-reconstruction.md) for the proposed experimental modeling backend.
- [Small diffusion and Apple Silicon plan](small-diffusion-apple-silicon-experiments.md) for compact local-model research.
- [Diffusion with heterogen context plan](diffusion-heterogen-gap-reconstruction.md) for class-specific research on retaining and conditioning on small molecules, cofactors, ions/metals, glycans, PTMs, and covalent ligands during protein-gap reconstruction.
- [Diffusion runner environment plan](diffusion-runner-environments.md) for separate core, Apple MPS, and Linux/NVIDIA dependency boundaries and production protocol wrappers.
- [Continuation context](api-ddd-continuation-context.md) for historical evidence only.
- [Approved installation-profile handoff](installation-profiles-handoff.md) for Standard/Desktop, the conda-family helper and clean-install acceptance.

Status is binary: `[x]` means complete and verified; `[ ]` means remaining,
partially complete, optional but not accepted, or not yet verified.

## Installation profiles: Standard/Desktop

Added 2026-10-28. Requirements and next-session starting point:
[implementation handoff](installation-profiles-handoff.md).

- [x] Record user approval and the correction that PySCF/RESP and MDAnalysis/Plotly remain in Standard/Desktop; Desktop inherits the scientific composition.
- [x] Merge Probe/xTB work into main at `cc0716c`, verify the merged baseline, and create `feature/installation-profiles`. This is baseline verification, not installer acceptance.
- [x] Define the dependency/feature matrix, one composition source and profile/readiness contracts (WP1).
- [x] Align Standard/Desktop/add-on specifications and Python metadata, including a verified single Open Babel provider (WP2).
- [x] Implement interactive/noninteractive helper and micromamba/mamba/conda adapters without modifying existing environments implicitly (WP3).
- [x] Implement protected MODELLER setup, Desktop build and separate diffusion-runtime preparation/registration (WP3); native MPS execution remains not-tested under WP5.
- [x] Add profile-aware doctor/preflight for scientific commands, ANARCI/HMMER/data, MSA, analysis and included PySCF RESP (WP4).
- [x] Execute feasible Linux WP5 acceptance: final six manager/profile installs, slow structural regressions, extended formats/chemistry/topologies/analysis, GROMACS and separate QM.
- [x] Synchronize the Linux WP6 candidate: generated Standard alias, Desktop inheritance, explicit add-ons/compatibility input, CI/manual/docs/agent metadata.
- [x] Execute native Apple micromamba 2.9.0 Standard/Desktop licensed 14-step smokes, Standard extended 27-step acceptance, installed Desktop's 405 tests and a real Chrome workspace/scientific job.
- [x] Extend native Apple acceptance with mamba 2.10.0 profiles/installed GUI, MSA alternatives plus PSI4 extended checks, separate PSI4 RESP, and one digest-verified public cc94 MPS installation bundle (not general scientific qualification).
- [x] Verify a pinned manual native Probe source build, exact-version doctor, real electron/nuclear diagnostics and the absent-Probe fallback; keep the managed conda add-on packaging restriction explicit.
- [x] Extend native cc94 MPS installation evidence with five predeclared public model seeds and a public Homology template-plan case, including actual generated-only CPU refinement and digest/geometry-verified publication; retain broader-context qualification as pending.
- [x] Recheck Linux `3309bd8` with new installed Standard (25 extended checks), managed Protenix/cc94 runtimes and 15 independently validated private CUDA candidates, including Protenix 152/512/513 size controls; retain failed publication and interrupted 9CQD evidence inside the project.
- [x] Resolve project CephFS publication with an exclusive manifest-last commit: recover 15 independently revalidated candidates into 11 readable bundles and pass fresh public Protenix batch/resume/ZBS/Homology on the project filesystem.
- [x] Complete the supplied full-FASTA 9CQD seeds 101–105 with retained raw continuation for seed 101, three accepted and two explicitly rejected bundles, independent movement/geometry/trace audits and unchanged originals.
- [x] Final Linux installed regression: new Desktop with 407 GUI tests and real Chrome workflow/viewport, 30 extended scientific checks (MSA/Probe/GROMACS/PSI4 included via explicit tooling), all 33 slow tests and a complete installed suite (1831 passed; two Torch skips passed separately). Execute cc94 batch/resume/ZBS/Homology and Protenix repair-before-full-context model/Homology.
- [x] Return to native Mac at `8f0cd33`: new intact Standard/Desktop, 1797 fast tests, 28/14 scientific smokes, installed GUI407/Chrome, native runner/publication87 and public cc94 MPS batch/resume/ZBS/Homology. Preserve the incomplete old-prefix attempts.
- [x] Execute all 33 native slow cases (32 together, optional PSI4 separately configured), plus the three Torch/Biotite fast skips in isolated test environments; preserve the original summary and dependency-failure attempts.
- [x] Run a separately prepared protein-only 5J13 baseline with original heavy coordinates retained: authoritative seven-gap/42-residue target, five executed MPS seeds, one accepted and four geometry-rejected outcomes. Do not count template-built sidechains as diffusion-generated atoms.
- [x] Extend predeclared native cc94 MPS case evidence with 1FR2/1A3N and a distinct public-altloc-selected 1DQJ derivative (5/5 each, zero fixed-heavy movement); preserve initial 1DQJ refusal and protein-only/observed-sequence limits instead of claiming broad scientific qualification.
- [x] Provide an explicit catalog-composed `ambertools-compat` route retaining all Standard science: native new-prefix helper/doctor, `pip check` and 28 extended checks pass with metadata-derived Python/NumPy/Biopython constraints and pdb2pqr. Default composition and separate runtimes remain unchanged; Linux acceptance of this optional route is not inferred.
- [ ] Resolve exact frozen seed-7 hardware-repeat provenance: current v8 repeats are deterministic/valid, but the unchanged historical v6 expected PDB hash is not reproduced; neither current-version repeats nor an arbitrary v6-label checkout closes that gate.
- [ ] Verify clean Standard/Desktop installs and add-ons with real functional smokes for each claimed manager/platform; retain not-tested states (WP5).
- [ ] Synchronize manual installation, documentation, CI and agent status, then finish the proven GUI/MSA migration without removing Standard scientific capabilities (WP6).

WP1–WP4 and the feasible Linux acceptance/migration are implemented. Native Apple
micromamba/mamba Standard/Desktop, browser, MSA and optional PSI4 acceptance passed.
Native conda package collisions, unavailable Probe packaging and broader cc94
MPS scientific qualification remain open and block whole-plan completion.
The exact `05c16c9` historical hardware source was also replayed twice on Linux:
both valid results are identical, but still differ from the old expected hash.
Keep that failed Linux forensic gate and upstream AmberTools auxiliary `pip check`
metadata distinct from the available native Mac work; do not regenerate a
favorable historical baseline after execution.
Continue using the [Apple acceptance handoff](installation-profiles-apple-acceptance.md).
Actual clean-install
evidence are tracked in [profile acceptance](../installation-profile-acceptance.md).
xTB remains removed under ADR 0012; Probe stays optional. Old core/full/extras
proposals are not the approved installation interface.

## Convenient Desktop GUI launch

Added 2026-11-03. Normal Desktop startup must not require users to remember
installation paths or compose `micromamba run -p ... npm --prefix ... start`.

- [x] Add `dvbfixer gui` to launch the installed GUI from the active Desktop
  environment. Resolve GUI assets and Node from that environment, not the
  current working directory or an unrelated tools prefix. Standard-only installs
  must provide an actionable explanation when GUI assets are unavailable.
- [x] Let the Desktop installer create a macOS application shortcut and a Linux
  desktop launcher that remember their installation prefix, start the local
  server and open the browser without prior shell activation. Make shortcut
  creation explicit and do not overwrite unrelated user applications/launchers.
- [x] Preserve loopback-first hosting, authentication and existing shutdown/
  child-process cleanup rules. Handle an already-running instance, occupied
  ports, startup failure and missing runtimes with clear diagnostics; do not
  open a browser before readiness or silently expose the service remotely.
- [x] Make the one-command/shortcut route the primary installation quick start;
  retain the manager/npm command only as an advanced diagnostic fallback.
- [ ] Test prefix discovery, launch outside the repository, Standard-only
  refusal, startup/readiness and shutdown, plus real installed Desktop command
  and shortcut smokes on macOS/Linux before claiming platform support.

Implementation slice 2026-11-05 (`feature/desktop-gui-launcher`): CLI-only
registry entry (not a remotely managed scientific command), prefix-bound Node/
Python/build lookup, inherited readiness pipe, private authenticated local
instance control, `--stop`, opt-in non-overwriting shortcuts and revised quick
starts are implemented. Unit/process tests cover reuse, failed/timeout startup,
signal cleanup, malformed control requests, isolation and helper ordering.
The installed-package/built-server macOS smoke passed local/authenticated arms
and direct execution of the generated app's executable. Native macOS
LaunchServices bundle activation and explicit stop also passed with browser
opening suppressed; evidence is retained under
`tmp/desktop-launcher-native-mcuzvs1z/evidence-v5/`. It used an isolated
wheel install and explicitly linked local Node/dependency fixtures, not a fresh
complete scientific Desktop solve. The first smoke's public-health auth
assumption failed and remains retained; the corrected test uses `/api/session`.
Manual Finder/default-browser UX acceptance and native Linux desktop/installed
command smokes remain open. This is not a fresh Standard/Desktop scientific
installer solve. Reproduce with
`python scripts/smoke_desktop_launcher.py --prefix PREFIX --output NEW_DIR --shortcut`.
Local verification: 614 tests in the CI-selected Python contract set, 443 GUI
tests, GUI production build/server smoke, lint/typecheck, strict API mypy and
generated CLI/GUI consistency passed. These do not substitute for the remaining
native Linux or full scientific installer acceptance gates.

## Geometry-gated anchor relaxation

Added 2026-10-27. Detailed execution requirements and evidence boundaries:
[headless handoff](geometry-gated-anchor-relaxation-handoff.md).

- [x] Record the 9CQD uniform/B-weighted positional-restraint pilot and visual-review exports.
- [x] Record user-approved geometry-triggered stages, raw-input restart, default-enabled anchor repair and explicit notifications.
- [x] Add initial pure B-weighting/fallback, explicit movement-bound measurement, staged-controller and scientific anchor-force helpers with narrow regression tests.
- [x] Add an independent research audit with explicit per-anchor limits, raw-input bindings, old strict-result preservation, outer-junction and observed-record checks; verify no-op native 1UBQ controls.
- [x] Replay all eight corrected-cc89 pilot cases with metadata/B provenance checks and shared-minimization accounting: five raw no-ops, three generated-only repairs, no observed-anchor motion.
- [x] Add bounded typed numerical-failure recovery and independent raw-based failure auditing; reject unclassified failures and contradictory/misbound evidence in regression tests.
- [x] Implement pure B-provenance/weighting/fallback policy and minimal problematic-anchor selection.
- [x] Extend the scientific refiner with explicit observed-anchor mobility and source-centred positional restraints.
- [x] Implement generated-only stage 1 and default-enabled, opt-out stage 2 restarted from raw coordinates.
- [x] Version validation/provenance and preserve identity, nonmobile coordinates, geometry and all-L invariants.
- [x] Thread policy, opt-out, diagnostics and resume compatibility through CLI, runners, ZBS, Homology and GUI.
- [ ] Evaluate the new chain-normalized B formula against uniform and calibrate per-anchor movement bounds on a broader cohort.
- [x] Complete regression and workflow verification plus native 9CQD refinement replay; update user docs and generated schemas for the experimental integrated path (not broader scientific/default promotion).

The experimental software path is integrated; see
[native replay and custody verification](../research/staged-anchor-integration-v1.md).
The remaining broader-cohort item is scientific follow-up: numerical defaults
are approved engineering limits, not calibrated experimental uncertainty.

### Actionable diffusion rejection output

Added 2026-11-03 following GUI review of the 9CQD cc94 seeds 7, 42 and
101. This is a reporting improvement, not a change to scientific acceptance
criteria. Currently the CLI prints hard-gate identifiers while concrete
residue-level findings remain in the refinement JSON report.

- [ ] Print human-readable final rejection findings in standard output, in
  English, without requiring `--verbose` or manual JSON inspection. Include
  the case-sensitive chain/residue identity, relevant atoms or metric,
  measured value, units and allowed range/limit. Example output:

  ```text
  Diffusion candidate REJECTED.
  W/ASP28: N–CA–C bond angle = 79.1°; allowed range: 95–125°.
  W/VAL29: backbone RMSD = 0.6001 Å; maximum allowed: 0.5000 Å.
  ```

  These illustrate separate retained candidates; do not combine measurements
  from different seeds or stages into one candidate's report.
- [ ] Carry authoritative Python-owned findings from refinement/independent
  validation to the CLI; do not duplicate geometry calculations in Node or
  substitute hard-gate identifiers for concrete findings. Preserve identifiers
  for machine-readable reports and provide an explicit fallback for failures
  without a residue-level finding.
- [ ] Separate remaining final rejection causes from intermediate findings
  already resolved by refinement. Identify candidate/seed and selected stage;
  keep raw and stage-1 details available in the full report.
- [ ] Print an actionable retained report path and explicitly state that a
  rejected structure is saved for inspection, not accepted. Avoid pointing
  users only to a temporary path removed during GUI run-storage compaction.
- [ ] Ensure GUI run output and captured logs show the same concrete English
  findings through the shared CLI/fd-level diagnostic capture. Clearly
  distinguish process completion from scientific candidate rejection; any
  exit-code/job-status compatibility change requires an explicit decision.
- [ ] Add regression coverage for bond angles, peptide planarity,
  Ramachandran outliers and per-anchor movement limits, including multiple
  candidates/stages, insertion codes and case-sensitive chain IDs. Verify
  standard output without `--verbose`, retained log/GUI delivery, no resolved
  intermediate findings presented as final failures, and no empty rejection
  banner for accepted candidates. Keep validation thresholds unchanged.

### Return to the real-complex findings

Evidence: [pre-anchor benchmark](../research/protein-complex-pre-anchor-benchmark-v1.md).
Execution and acceptance criteria are in the
[follow-up work packages](geometry-gated-anchor-relaxation-handoff.md#real-complex-follow-up-work-packages).
These are open follow-ups, not completed by merging the implementation.

- [ ] **1BRS D:64–65:** run an explicitly declared independent/sparse-context test with incomplete GLN D:58; preserve absent OE1/NE2 as unobserved rather than silently repairing the experimental source. Keep the original single-gap input refusal as a distinct result.
- [ ] **1GC1 G:397–409:** test resampling or fixed-anchor conformational search for the residual generated SER G:398 Ramachandran outlier. Retain the original failed seed; do not loosen geometry gates or release observed anchors without independent eligibility.
- [ ] **1HE8 A:322–352 and A:528–619:** profile the 240 s generated-only timeouts, then preregister a separate resource-budget experiment. Distinguish incomplete computation from geometry failure and preserve all raw evidence.
- [ ] **Live raw/final integration:** rerun a declared native MPS workflow after the raw-diagnostic fix, retaining complete sampler resource/probe provenance through refinement, publication and resume. Archived coordinate recovery does not close this acceptance item.

## API Phase 0: Decisions And Baseline

### Accepted Decisions

- [x] Accept an ADR making Node/TypeScript the HTTP, workspace, artifact, and job shell.
- [x] Accept an ADR making Python the sole owner of scientific naming policy.
- [x] Accept an ADR for host-neutral route composition rather than a duplicate Python web service.
- [x] Record `dvbfixer atom-names` as the dedicated naming command.
- [x] Record naming conversion as synchronous and long-running workflows as asynchronous jobs.
- [x] Record workspace artifact IDs as the only V1 server-side input reference.
- [x] Record the V1 policy rejecting multiple `MODEL` blocks.

### Open Maintainer Decisions

- [x] Select an authentication model.
- [x] Define principal-to-workspace authorization semantics.
- [x] Decide whether V1 needs direct upload/download conversion.
- [x] Decide whether HTTP should accept a JSON override artifact as well as inline overrides.
- [x] Decide whether `tleap-reduce` must honor the force-field naming option.

### Baseline Evidence

- [x] Add regression coverage for naming blockers 1-6 and 9 from the API roadmap.
- [x] Complete end-to-end nucleic-acid naming coverage.
- [x] Capture and document any still-needed representative golden inputs.
- [x] Complete ADR and test evidence distinguishing current behavior from intended V1 behavior.

## API Phase 1: Naming Application Service

### Application Boundary

- [x] Add typed naming request, result, diagnostic, target, profile, override, identity, and report vocabulary.
- [x] Implement `convert_force_field_naming` as PDB text in and transformed text plus report out.
- [x] Keep file I/O outside the pure service.
- [x] Retain the atomic in-place compatibility adapter for existing pipeline callers.
- [x] Report changed atoms, residues, variants, and atom names separately.
- [x] Keep naming conversion distinct from canonical `rename` and chemical parameterization.

### V1 Semantics And Correctness

- [x] Accept PDB text and support AMBER and CHARMM with the GROMACS profile.
- [x] Preserve exact, case-sensitive chain IDs and insertion codes.
- [x] Use typed three-field override identities and restrict legacy two-field identities to blank insertion codes.
- [x] Apply explicit override, explicit source variant, then no inferred change precedence.
- [x] Do not infer HID/HIE/HIP from canonical HIS alone.
- [x] Validate before returning or committing output and leave sources untouched on failure.
- [x] Reject multiple `MODEL` blocks.
- [x] Detect source and target atom-name collisions per residue.
- [x] Update correlated `ANISOU` and parseable full `TER` records.
- [x] Preserve unrelated records and non-name atom fields.
- [x] Support four-character CHARMM names and LYN/LSN hydrogen pairs idempotently.
- [x] Remove ambiguous two-field override fallback from the public service.
- [x] Add end-to-end DNA and RNA naming tests.
- [x] Document unvalidated terminal 5-prime/3-prime hydroxyl cases as unsupported.
- [x] Resolve `tleap-reduce` naming parity.
- [x] Document `standard` as suppressing shifts, not reversing GROMACS names.
- [x] Do not claim naming compatibility proves template or bond compatibility.

### Test Matrix

- [x] Verify every `GROMACS_AMBER_ATOM_RENAMES` mapping.
- [x] Verify AMBER and CHARMM caps and termini.
- [x] Verify exact coexistence of residues such as `H:82` and `H:82A`.
- [x] Verify one-residue and multi-chain proteins.
- [x] Verify byte-identical repeat conversion across all supported mappings.
- [x] Complete malformed-line coverage.
- [x] Complete `ATOM`, `HETATM`, `ANISOU`, `TER`, and untouched-record matrix coverage.
- [x] Verify empty/no-op conversion and report counts.

Evidence: `tests/test_force_field_naming.py` checks independent expected AMBER
spellings against the complete mapping inventory, both coordinate record types,
caps/termini, variant directions/profiles, DNA/RNA, repeated conversion, malformed
identity fields, correlated records, line endings, and report counts. This
verifies naming transformations, not chemical template or bond compatibility.

## API Phase 2: Dedicated CLI Adapter

- [x] Add `dvbfixer atom-names` with one input and one explicit output.
- [x] Add `--target-ff {amber,amber19sb,charmm}` and `--profile {gromacs,standard}`
  (native profile exposed publicly in 0.9.1; GROMACS remains the default).
- [x] Support variant overrides from JSON, `--dry-run`, and `--report-json`.
- [x] Call the pure naming service and keep the command distinct from `rename` and `convert`.
- [x] Register the command and regenerate CLI and GUI command references.
- [x] Add subprocess-level tests and deterministic conversion.
- [x] Explicitly verify source preservation at CLI level.
- [x] Explicitly validate the report schema at CLI level.
- [x] Explicitly verify non-zero exit and no destination on unsafe input.
- [x] Explicitly verify the selected CIF boundary behavior.

## API Phase 3: Versioned Naming API

### Core Vertical Slice

- [x] Add a V1 naming-conversion route with TypeBox runtime schemas and OpenAPI 3.1.
- [x] Resolve inputs by workspace artifact ID.
- [x] Invoke the Python command instead of duplicating naming rules in Node.
- [x] Validate bounded command reports and output digests.
- [x] Reload the current manifest before publication.
- [x] Register exactly one output with reproducible provenance and none for dry runs or failures.

### Request And Artifact Transaction

- [x] Authenticate the caller.
- [x] Authorize the principal before workspace path resolution.
- [x] Validate the request body and resolve `inputArtifactId` through the manifest.
- [x] Explicitly verify accepted source artifact types.
- [x] Enforce source artifact size limits before execution.
- [x] Verify operation-directory privacy and cleanup behavior.
- [x] Parse the report, validate output existence/digest, and publish against the current manifest revision.
- [x] Return the output artifact and report.
- [x] Define and verify retained failure-log behavior.
- [x] Add bounded retention or pruning for failed operation directories.

Evidence: `gui/server/naming-retention.test.ts` and `naming-api.test.ts` verify
age/count pruning, workspace and symlink isolation, retained diagnostics,
configuration validation, and cleanup on subsequent successful execution.

### Provenance, Observability, And Tests

- [x] Record source, target, profile, overrides, report, and digest provenance.
- [x] Explicitly record command/service and DVBFixer package versions.
- [x] Add server-owned request IDs and structured API access logs.
- [x] Add bounded process-local service metrics.
- [x] Add durable audit events and retention.
- [x] Generate a client from OpenAPI.
- [x] Verify stable validation and error envelopes across V1 routes.
- [x] Add authentication and authorization tests after implementation.
- [x] Explicitly test path containment, source/request size limits, and naming timeouts.
- [x] Test artifact lookup, dry runs, one-artifact publication, concurrent manifests, and failures.
- [x] Validate OpenAPI examples against runtime schemas.

## API Phase 4: Standalone Server

### Local Standalone Host

- [x] Extract host-neutral route composition and keep Vite as a thin adapter.
- [x] Add a bundled Node entry point serving the built GUI and API.
- [x] Bind to loopback by default and require explicit insecure remote acknowledgment.
- [x] Validate host, port, roots, and shutdown configuration.
- [x] Reject lexical and symlink-resolved static/data root overlap.
- [x] Return JSON for unknown API paths and use SPA fallback only for HTML requests.
- [x] Stop intake, reject new work, terminate children, then drain HTTP.
- [x] Bound the emitted-server startup and shutdown smoke test.
- [x] Complete final review and verification of the standalone slice.
- [x] Commit the standalone slice.

### Public Deployment Hardening

- [x] Add authentication and principal-to-workspace authorization.
- [x] Add restrictive CORS.
- [x] Add upload and workspace limits plus process-wide child concurrency control.
- [x] Bound the process-wide child admission queue and reject overload before spawning.
- [x] Add bounded process-local request-rate limits before authentication.
- [x] Add OS-level CPU, memory, and filesystem quotas.
  - [x] Add a fail-closed systemd/cgroup-v2 deployment profile and bounded-storage preflight.
  - [x] Run privileged kernel-enforcement acceptance tests on the Ubuntu 24.04 reference deployment host in CI.
- [x] Document current local deployment and security boundaries.
- [x] Document the supported single-host public-deployment profile after security controls exist.
- [x] Add storage-level manifest locking or compare-and-swap.
- [x] Replace process-local locks and event subscribers for multi-instance operation.
- [x] Complete reference public-deployment smoke tests through HTTPS on the CI VM.

## API Phase 5: General DVBFixer API

### Route And Contract Policy

- [x] Version managed-job routes under `/api/v1`.
- [x] Retire the duplicate synchronous generic command route.
- [x] Classify every operation as a deterministic transform or managed workflow.
- [x] Keep minimization, modeling, and parameterization out of synchronous APIs.
- [x] Preserve cancellation for long-running operations.
- [x] Generate clients and contract tests from OpenAPI.
- [x] Define API compatibility, versioning, and deprecation policy.
- [x] Add durable scheduling, state, and event delivery before horizontal scaling.

### Operation Rollout

- [x] Expose `diagnose`.
- [x] Expose `conect`.
- [x] Expose `renumber`.
- [x] Version or migrate managed `prepare`, `minimize`, `model`, and `zbs`.

## GUI tool navigation and grouping polish

Added 2026-11-03. User-requested navigation changes; scientific command
behavior stays unchanged.

- [ ] Move the **Altloc** tab immediately before **Model** in tool navigation.
- [ ] Move **Pull** to **Utilities**, with adjacent order **Conect → Puppet → Pull**.
- [ ] Move **Cluster** from **Analysis** to the current **Glycoprotein
  Preparation** group, with order **Convert → Transplant → Cluster**.
  Choose a broader, descriptive name for the destination
  group that fits its complete tool inventory, including Cluster; confirm the
  final name before implementation rather than silently selecting one.
- [ ] Move **Diagnose** to the first position in **Utilities**, and remove the
  now-obsolete **Analysis** group. Preserve tool identities, saved tool-panel
  state and existing workspace/job references while changing presentation.
- [ ] Add cascading tool submenus to the utility-group dropdown: hovering a
  group reveals its tools and selecting a tool opens it directly. Support
  keyboard navigation/focus and click/touch activation as equivalent routes;
  keep submenus reachable across pointer transitions and within viewport
  bounds, with clear active-tool/group indication.
- [ ] Keep group names, ordering and tool membership consistent across the
  navigation registry and other generated/shared views; do not hand-edit
  generated specifications. Follow the shared GUI design system.
- [ ] Add navigation regressions for exact ordering/membership, removal of
  Analysis, retained panel state, submenu selection, pointer transitions,
  keyboard and touch/click operation. Run GUI typecheck/tests and inspect the
  resulting navigation at desktop and narrow viewport sizes.

## GUI diffusion multi-seed input

Added 2026-11-03. The CLI supports repeated integer options such as
`--diffusion-seed 1 --diffusion-seed 17 --diffusion-seed 148`, but the GUI does
not make entering multiple seeds clear. Fix usability and end-to-end transport,
not just the field label.

- [ ] Provide an explicitly labelled **Seeds** control for diffusion, with
  visible help and examples. Accept comma-separated and whitespace-separated
  integers (`1, 17, 148` or `1 17 148`, including mixed separators) and normalize
  them to an ordered integer list. Show the default seed 7 when unspecified.
- [ ] Validate entries before submission: reject negative, fractional,
  malformed, duplicate and non-safely-representable values with actionable
  feedback; preserve input order and match Python's unique non-negative seed
  contract. Never silently discard invalid entries or pass a whole list as
  one integer argument.
- [ ] Trace all GUI diffusion entry points, including Model and Homology where
  applicable, through client/API schemas, CLI argument construction and saved
  panel state. Serialize the normalized list as repeated `--diffusion-seed`
  arguments and preserve it in job provenance and reload/resume behavior.
  Keep single-seed inputs and existing persisted values compatible; do not
  edit generated GUI specifications by hand.
- [ ] Explain that each seed requests a candidate, not necessarily a separate
  concurrently scheduled job. Show requested seeds and per-candidate
  completion/validation outcomes; integrate with Job Monitor without implying
  that seed count guarantees parallel execution. Keep existing whole-run
  timeout semantics visible for multi-seed requests.
- [ ] Add parser, UI and API/argv regression tests for both separator styles,
  mixed separators, empty/default and single-seed inputs, invalid/duplicate
  values, persistence and exact repeated-option transport. Verify a real
  multi-seed GUI workflow retains a separately identified result/outcome for
  every requested seed, including rejected candidates; run GUI typecheck/tests
  and relevant backend tests before marking complete.

## GUI Job Monitor

Added 2026-11-03. Planned dedicated **Job Monitor** tab for tracking task
execution, with a design that supports multiple simultaneous jobs rather than
one global active-task indicator. This is planned work, not a claim that all
workflow types currently permit concurrent execution in the same workspace.

- [ ] Add a Job Monitor tab next to **Info** in the same docking tabset,
  listing queued, running and recently finished jobs.
  Show job identity, tool/workflow, owning workspace, start time, elapsed time,
  current execution stage and terminal outcome. Keep process completion and
  scientific acceptance/rejection distinct.
- [ ] Show per-job stage progression and progress using authoritative runtime
  events. Display completed/total work and percentages only when measurable;
  otherwise show an indeterminate indicator, stage name and elapsed time.
  Distinguish sampling, refinement, validation and publication where supported;
  do not invent progress or expose a partial output as ready.
- [ ] Define a shared, versioned job-stage/progress event contract and connect
  supported Python workflows and external runners to the existing managed-job
  transport. Reuse persisted job state and event delivery rather than adding
  a separate GUI-owned scheduler or inferring scientific stages from fragile
  free-text matching. Keep logs available for detailed inspection.
- [ ] Support independent selection, progress and cancellation controls for
  multiple jobs. Review existing workspace locks, concurrency limits and
  CPU/GPU/memory admission before enabling parallel execution; queue conflicting
  jobs and explain why they are waiting. Preserve workspace/artifact isolation
  and do not remove safety locks simply to permit concurrency.
- [ ] Provide an accessible per-job **Stop** action through the existing
  cancellation API. Show cancellation requested/stopping until termination is
  confirmed; terminate tracked subprocesses, handle completion/cancellation
  races, retain failure/cancellation diagnostics and never publish incomplete
  results as accepted. Stopping one job must not stop unrelated jobs.
- [ ] Restore monitor state after tab changes, reloads and event reconnects
  from authoritative persisted records. Define cross-workspace visibility with
  ownership/ACL enforcement; switching workspaces must not misattribute running
  jobs or leak another principal's job data.
- [ ] Link each job to its logs, actionable validation findings, retained run
  reports and completed outputs. Follow the shared GUI design system and keep
  monitoring usable while a different task or workspace is open.
- [ ] Test concurrent jobs, queued/resource-blocked jobs, stage/progress events,
  indeterminate work, reload/reconnect, authorization, cancellation of one job
  without affecting others, process cleanup and rejection versus completion.
  Run GUI typecheck/tests plus backend contract tests and a real multi-job
  start/monitor/stop smoke before marking this feature complete.

## API-Wide Quality Gates

- [x] Use one runtime-schema source for every V1 route and OpenAPI description.
- [x] Ensure no public route accepts unrestricted filesystem paths.
- [x] Keep scientific naming policy implemented once in Python.
- [x] Record reproducible inputs, options, and versions for every output artifact.
- [x] Enforce workspace ownership before path resolution.
- [x] Return structured reports for every exposed operation.
- [x] Require durable workflow state before horizontal scaling.
- [x] Adopt one common V1 error envelope and explicit compatibility policy.
- [x] Test status mappings for 400, 401, 403, 404, 409, 413, 415, 422, and 500.
- [x] Stabilize documented naming error codes.
- [x] Run release verification: GUI checks, focused Python tests, generated-file checks, agent-doc checks, and the non-slow Python suite.
- [x] Push the completed implementation and require every configured GitHub Actions CI check on that commit to pass.

Verification on 2026-09-21: the [post-push CI run for `cf56aaa`](https://github.com/lilpoundcake/DVBFixer/actions/runs/35661988827)
passed all four configured jobs. The scientific lane reported 552 passed and
3 skipped (no Modeller license); GUI reported 232 passed. Locally, the
non-slow Python suite reported 533 passed and 3 skipped, and the Python 3.11
fast-lane selection reported 220 passed. Ruff, mypy, TypeScript, generated
CLI/GUI/OpenAPI files, agent-doc validation, and `git diff --check` passed.
Python 3.13 changes argparse help formatting, so CI pins its scientific lane
to Python 3.11. This verification does not include licensed Modeller
integrations or target-host public-deployment acceptance.

Reference deployment acceptance on 2026-09-22: the
[Ubuntu 24.04 CI run for `b72e477`](https://github.com/lilpoundcake/DVBFixer/actions/runs/35695278505)
passed its new privileged deployment job. The systemd service's mandatory
preflight accepted a dedicated ext4 data mount, read-only root, private tmpfs,
CPUQuota=200%, MemoryMax=16 GiB, MemorySwapMax=0, and TasksMax=512. Scaled
probes recorded `cpu.stat nr_throttled` 1→31, `memory.events oom_kill` 0→1,
and `pids.events max` 0→2; writes stopped at the data and temporary filesystem
limits, and stopping the unit removed a session-detached child. HTTPS smoke
covered health, required bearer auth, CORS rejection, and workspace creation;
MAFFT, tleap, and Reduce produced output under a contained unit. These results
validate the Ubuntu reference profile, not an operator's future production host;
the installation-specific acceptance in `docs/deployment.md` remains required.

The versioned workflow API uses a durable workspace run lock, cross-process
cancel markers, and persisted job-record polling for SSE on a shared local
filesystem. `docs/adr/0009-versioned-workflows-and-compatibility.md` defines
the V1 error and deprecation policy. Horizontal scaling across multiple hosts
remains explicitly unsupported until an external scheduler/state store exists;
these checked workflow-state items cover one-host multi-process operation.

## Diffusion Gap Reconstruction

- [x] Complete Phase 0: policy, license inventory, corpus, leakage metadata, and predeclared thresholds.
- [x] Complete Phase 1: backend-neutral contract, masks, Kabsch/reinjection primitives, fake runner, validation, provenance, and atomic publication.
- [x] Complete Phase 2: pinned RFdiffusion v1 Linux/NVIDIA benchmark adapter and minimal digest-pinned container recipe; final container build/GPU smoke awaits a host with an OCI runtime.
- [x] Complete Phase 3: select Protenix v1 through the all-atom constrained-sampler ablation and 231-group confirmatory benchmark; retain Boltz-2 as proxy-only.
- [x] Record the confirmatory result: Protenix 220/231, MODELLER 169/231, noninferiority demonstrated, and `candidate-preferred` decision.
- [x] Deliver the Phase 4 experimental `model --backend diffusion` implementation while retaining
  MODELLER as the default and without automatic fallback. Public CLI, production
  wrappers, additive Doctor preflight, fail-closed batch/GUI admission,
  record-preserving candidate materialization, and portable protocol failure
  handling are implemented and tested. Public cc89 CUDA and corrected Apple
  integration each pass 5/5 within their bounded scope.
- [ ] Complete the selected Protenix production-wrapper 231-case replay and final
  OCI image build/GPU smoke; profile-lock pinning and public GPU smokes already exist.
- [x] Deliver mosaic-first Homology with authoritative template coordinates,
  coverage ownership, CLI/GUI routing and bounded public integration evidence.
- [ ] Complete the larger independent Homology/context/edit evaluation before
  claiming broader scientific acceptance.
- [ ] Complete Phase 6: separate evidence-backed production decision.
- [x] Execute the macOS Apple Silicon CPU/MPS portability evaluation separately
  from Linux/NVIDIA acceptance and remote-runner evidence; all 231 frozen cases
  reached terminal outcomes with fallback disabled, with 213/231 selected passes.
- [x] After native GPU gates pass, define the external NVIDIA Docker runner with a digest-pinned base and digest-verified mounted checkpoint.
- [ ] Build and GPU-smoke the final image when an OCI runtime is available; defer registry signing, multi-arch, embedded weights, and orchestration until deployment creates those requirements.
- [x] Complete compact-model Linux scientific selection, the 231-case
  Protpardelle descriptive follow-up, and the frozen Apple Silicon CPU/MPS
  portability matrix on the physical Mac.
- [x] Complete Protpardelle-1c artifact audit, Linux CPU/CUDA operator smoke,
  initial three-residue withheld-coordinate pilot, and frozen six-case v5 pilot;
  generated-only CPU refinement raises the latter from `2/6` to `6/6` hard-gate
  passes, while training-membership leakage remains unresolved.
- [x] Complete the preregistered 24-case compact-model extension: raw
  Protpardelle passes `13/24`, refinement alone passes `21/24`, and validation-
  first raw/refined selection passes `22/24` versus Protenix `23/24` and MODELLER
  `21/24`; keep the result exploratory because training leakage is unresolved.
- [x] Add and exercise the optional Protpardelle post-update callback. Corrected
  internal-frame reinjection is byte-identical across two same-seed `9dvi` runs
  and exact for all 569 represented fixed atoms, but does not improve raw geometry;
  OpenMM CPU refinement passes while remaining measurably nondeterministic.
- [x] Run corrected per-step reinjection over the frozen 24-case compact follow-up.
  Exact callback control passes all operational checks, but validation-first
  coverage is `21/24` versus `22/24` for final-only projection, so per-step
  reinjection is not selected as the compact baseline.
- [x] Diagnose OpenMM's CUDA PTX failure and add an isolated CUDA-12.2-compatible
  refinement environment with digest-verified ff19SB data. Two byte-identical
  `9dvi` runs and the 440-residue `9ina` run pass every CUDA hard gate. This historical
  kernel diagnostic preceded the completed 231-case follow-up recorded below.
- [x] Freeze the full 231-case Protpardelle descriptive-follow-up manifest before
  cohort-wide inference, with 128 gap-5 and 103 gap-10 cases, CUDA refinement,
  validation-first selection, digest-tracked code/environment inputs, and no
  replacement of failures.
- [x] Preserve the stopped v5 full-run manifest after its internal CUDA
  fixed-particle tolerance rejected `35zr`; pass the v6 standalone gate and
  freeze a successor manifest before a clean namespace restart without changing
  sampling, physical refinement, validation, selection, or cohort membership.
- [x] Complete the clean v6 231-case follow-up: Protpardelle raw passes `101/231`,
  CUDA-refined and validation-first selected output pass `217/231`; retain Linux
  evidence as immutable comparator data for the Apple portability study.
- [x] Complete the staged Apple Silicon handoff: explicit CPU/MPS device
  selection, fallback-disabled operator smoke, same-seed repeat, six-case gate,
  24-case portability pilot, 100-case operational soak, and full frozen 231-case
  run. Evidence remains descriptive because training membership is unresolved.
- [ ] Benchmark at most two adapted pretrained compact candidates before training a purpose-built local gap model.
- [ ] Train and freeze a leakage-controlled gap-only model only after the public-structure dataset and sequence-cluster splits are immutable.
- [ ] Require withheld-coordinate noninferiority and Apple memory/latency gates before any local experimental backend decision.

Detailed requirements and acceptance gates are in the
[diffusion gap-reconstruction plan](diffusion-gap-reconstruction.md).

### CUDA expansion follow-up (2026-10-19)

- [x] Share profile-specific gap-length admission across CLI, mosaic construction,
  and scientific scope; preserve separate public and research capabilities.
- [x] Independently audit withheld accuracy of all 50 full-input clean-protein
  CUDA outcomes, including rejected candidates and native-quality exclusions.
- [x] Validate adaptive queries and previous source context; do not report a
  context expansion when the policy selects no additional atoms.
- [x] Exercise deletion flanks 1–5 and reject invalid Homology requests before
  writing the normalized candidate input.
- [x] Freeze and execute generated-only CUDA refinement for all six rich cc89
  frontier outputs; retain 0/6 acceptance and two guarded chirality exceptions.
- [x] Diagnose the original catastrophic cc89 failure with frozen tensor/axis
  controls and implement an explicit CUDA-only neutral denoiser chain-phase mode.
- [x] Confirm all six repaired candidates match causal controls byte-for-byte;
  unchanged generated-only CUDA refinement accepts 6/6 at 0.17–0.25 Å RMSD.
- [x] Complete the 13-source, 30-mask known-sequence follow-up: 60 cc89 raw and
  refinement outcomes, 30 MODELLER comparisons, four public Homology/Protenix
  controls, exact masked-template checks, and independent artifact revalidation.
- [x] Record the negative broad-scope result without changing thresholds or
  replacing failures: cc89 37/60 after refinement, including 0/16 CDR outcomes;
  retain three guarded chirality failures. Public Protenix controls pass 4/4.
- [x] Audit actual first-attention chain/position features, raw prefit tensors,
  frame conversions and native oracle round trips; distinguish packing from
  false peptide adjacency and from whole-complex sampling misuse.
- [x] Implement explicit CUDA source-physical position features and local-input
  correction; recheck all 16 original CDR mask/seed pairs with unchanged
  refinement: 16/16 pass, median backbone RMSD 1.152 Å. Preserve the original
  negative holdout and label the correction as diagnostic, not independent acceptance.
- [x] Complete and independently audit all 44 corrected general-protein outcomes:
  44/44 pass after unchanged refinement, retaining all 37 original passes and
  recovering seven failures. Combined corrected acceptance is 60/60, zero D-Cα,
  and exact outside-gap preservation. Record paired raw/final transitions;
  treat reference RMSD as secondary rather than an acceptance criterion.
- [x] Freeze and independently revalidate the first source-disjoint corrected-mode
  check: 9WW8, gaps 5/10, seeds 7/19, 4/4 after unchanged refinement. Keep the
  source denominator explicit and freeze a separate expanded CDR/general cohort.
- [x] Complete the frozen 8EK6/9Z0U/9LDN prospective cohort and independent audit:
  20/20 pass after unchanged refinement. Combined prospective acceptance is
  24/24 across four new sources, zero detectable D-Cα and exact outside-gap
  preservation; keep separate from the development recheck.
- [ ] Complete prospective long-gap noninferiority and independent review.
- [x] Expose corrected cc89 as the public Linux/NVIDIA `protpardelle-1c-cuda`
  profile for one internal 5–10-residue gap, automatic fixed-policy local context,
  isolated CUDA refinement, CLI/GUI routing, exact bundle validation and resume.
  Public model/Homology integration accepts 5/5 candidates across four workflows.
- [ ] Complete wider independent cc89 CUDA acceptance beyond the corrected
  development recheck, bounded four-source prospective check and public integration.
- [ ] Complete larger independently reviewed context/CDR and Homology edit cohorts.

Evidence: [CUDA expansion completion audit](../research/cuda-expansion-completion.md).
Follow-up: [cc89 chain-phase diagnosis](../research/cc89-chain-phase-root-cause.md).
Expanded follow-up: [known-sequence structural holdout](../research/cc89-known-sequence-holdout.md).
Input correction: [representation audit](../research/cc89-input-representation-audit.md).
Apple follow-up completed in `0e2a8c8`: corrected research pilot 8/8 and versioned
`protpardelle-1c-mps-local-v1` public integration 5/5, including multi-seed publication
and exact frozen-request replay; historical MPS semantics remain separate.
- [ ] Complete the larger corrected-Apple memory/latency/operational cohort before
  wider acceptance. See [Apple public evidence](../research/cc89-public-apple-acceptance.md).

## DDD Documentation Roadmap

### Phases 0-2: Foundation And Domain Contexts

- [x] Define source precedence, statuses, references, delegation, and ownership rules.
- [x] Use dependency-free TOML maps and validate syntax, links, IDs, and evidence in CI.
- [x] Map force-field naming callers, boundaries, contracts, invariants, and known gaps.
- [x] Add Force-Field Naming, Structure Identity, and Parameterization contexts.
- [x] Add their task, contract, invariant, and change-group records.
- [x] Give every `dvbfixer.domain` export an active meaning or explicit vocabulary-only status.

### Phases 3-4: Artifacts And Pipelines

- [x] Document `.dat`, structure identity, connectivity, workspace, Homology, and topology artifacts.
- [x] Record owners, producers, consumers, keys, mutation policy, atomicity, cleanup, and failure behavior.
- [x] Add Structure Preparation, Topology Generation, Diagnostics, and Workflow Execution contexts.
- [x] Document boundaries without copying algorithms and reconcile cross-context contradictions.

### Phases 5-6: Tasks And Validation

- [x] Index high-risk `.dat`, rebuild, CONECT, topology, workspace, API, and fixture changes.
- [x] Add new task records only when repeated high-risk work justifies them.
- [x] Validate TOML, metadata, paths, Markdown, linked records, evidence, sections, and indexed pages.
- [x] Add symbol checks, stale-commit reporting, generated views, or PR checks only when justified.
- [x] Keep a general Python/TypeScript AST index out of the initial rollout.

### Ongoing Maintenance

- [x] Use `docs/agent/README.md` as the canonical agent-map entry point.
- [x] Maintain shared maps, contexts, validator, and tests through one coordinating editor.
- [x] Separate implemented behavior from vocabulary, proposals, research, and gaps.
- [x] Surface unresolved contradictions instead of choosing silently.
- [x] Add context pages only when structured records cannot explain the concern.
- [x] Run agent-doc validation and `git diff --check` after every map change.
