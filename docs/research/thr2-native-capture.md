# Original-input THR2: bounded fresh native capture

## Prospective protocol and authorization

The user approved this separate bounded original-input/platform capture and
autonomous continuation on 2026-11-12. This supersedes the earlier authorization
boundary in [whole-path research](carbon-H-whole-path.md), not its findings.
Documentation scope: this note (protocol, results, limits) and the canonical
[handoff](../plans/pdbfixer-reconstruction-handoff.md); validate both links and
evidence before publishing. A follow-up link in the earlier whole-path note marks
its authorization boundary as superseded without altering historical findings.

Before scientific execution, `research_thr2_native_capture.py plan` binds the
original comparator plan, all 160 historical source hashes, original THR2
prepared PDB/Dat/log/output, OpenMM Python source and observer dependencies.
Execution resolves scientific imports from the preserved historical tree.

Exactly **three fresh native CLI executions**: unseeded, H seed 1 and H seed 2,
once each, automatic platform, `--no-solvent --atom-naming standard`, otherwise
historical defaults. Seed only Python random during `Modeller.addHydrogens`;
restore it for seeded arms, preserve its natural evolution in the unseeded arm.
Retain RNG states and before/after H coordinates. No initializer, force,
restraint, minimization tolerance/iteration or safeguard changes.

Observe every actual `LocalEnergyMinimizer.minimize`, including internal H
placement: preserve topology/bonds, System/Integrator/State XML, checkpoint,
actual platform properties, arguments and before/returned coordinates. No native
reporter or energy/force readback is added. Boundary serialization is itself
instrumentation; this is not proof of an identical uninstrumented trajectory.

After **all three native executions**, bind captured artifacts, then replay each
captured call independently **once on captured-checkpoint platform/properties and
once on Reference from XML State** with a read-only native callback reporter.
Never chain replay endpoints into the next phase. Native call parameters remain
unchanged. Record starts, callback arrays/objectives/gradients and returned state;
report OpenCL initial/endpoint differences, never demand cross-platform identity.

Bounds: at most eight native calls per label, 4,096 callbacks per replay,
1,800 seconds per job and 14,400 seconds total. All enrolled labels are attempted;
partial capture, failures, timeouts, missing calls and budget nonexecutions remain
evidence. No retry, adaptive extension, extra seeds or success filtering.

## Interpretation gates

The original failing THR2 run has no saved pre-minimization System/State/checkpoint.
These new captures cannot recover or replace that stochastic trajectory. Failure
to reproduce inversion is a valid result. Reference behavior cannot establish
OpenCL causality; callback observations are not every line-search evaluation.
Native error enums are not exposed by this observer. Production and upstream
issue publication remain outside this authorization.

Private payloads belong under `tmp/pdbfixer-research/thr2-native-capture-run-01/`;
only scripts, tests, bounded aggregate evidence and notes may be tracked.

The [frozen plan](thr2-native-capture-plan.json) is a byte-identical public copy
of `tmp/pdbfixer-research/thr2-native-capture-plan-01.json`; it contains hashes and
protocol metadata, not native coordinates, forces or checkpoint bytes.

## Results: fresh native inversion captured

All **3 native executions / 9 native calls / 18 independent replays** completed.
Every native call selected **OpenCL / Apple M3 Pro / single precision**. H calls
used native **1.0 / 50**, full phases native **10.0 / 1,000**. All outcomes were
retained; no additional seeds, repeats or minimizations were run.

| Fresh native label | Before whole phase 1 Cβ | After phase 1 | After phase 2 / final PDB |
|---|---|---|---|
| unseeded | R, +0.198389 Å³ | **S, −2.472333 Å³** | **S**, −2.444469 / −2.446819 Å³ |
| seed 1 | R, +0.198389 Å³ | R, +2.463277 Å³ | R |
| seed 2 | R, +0.198389 Å³ | R, +2.484431 Å³ | R |

Volumes use `(CB, CA, OG1, CG2)`; independent RDKit CIP labels agree. H placement
did not invert heavy stereochemistry. All 42 Cα centres remained positive at
native boundaries, every replay callback and returned state. The starting target
Cα/Cβ volumes were already small, **0.813111 / 0.198389 Å³**: this is not a
demonstration of inversion from an initially well-formed tetrahedral frame.

The three whole-phase-1 Systems and Integrators are byte-identical. All **327
heavy-atom coordinates and identities** match exactly between those three starts;
the **315 H coordinates differ** (seed 1/2 vs unseeded H RMSD **1.062439 /
1.049278 Å**). This isolates the saved-start difference to H coordinates, not
changed heavy geometry or force parameters. It does not isolate Python RNG from
OpenCL numerical variability in generating those H starts. Across each H build,
heavy coordinates moved by up to **1.282947×10⁻⁶ Å** relative to its input; do not
claim byte-exact original-observation preservation.

### Accepted-iterate inversion, not terminal leakage in this cohort

All **18/18 replay initial coordinates** match their saved native starts exactly.
No replay endpoint is bit-identical to its original native endpoint. OpenCL
maximum endpoint differences span **0.000192–0.111995 Å**; Reference controls are
not numerical reproductions of OpenCL. Nevertheless, both numerical modes
preserve the per-label stereo outcome, including unseeded R→S.

The unseeded whole-phase-1 inversion appears between **callback ordinals 2 and 3**:

| Mode | Cβ volume before → after (Å³) | Objective change (kJ/mol) |
|---|---|---|
| captured-checkpoint OpenCL | +0.026375 → −0.174341 | −72.616150 |
| Reference XML | +0.026369 → −0.174343 | −72.616189 |

Across **1,114 callbacks**, same-platform independent physical readback has zero
measured mobile-gradient and objective discrepancy. Every returned state equals
its last reported coordinates: **0/18 terminal mismatches**, unlike the earlier
THR24 Reference cohort. All 12 full-phase replay endpoints satisfy the native
gradient criterion; six H endpoints do not (their native 50-update budgets were
exhausted). Passing that numerical criterion is not stereochemical qualification.

Thus a fresh native R→S event is now captured, with matching accepted-iterate
inversion in both saved-start controls. **Terminal leakage is not the explanation
of this fresh event.** It does not invalidate the distinct THR24 terminal-leakage
findings or recover the original November 6 trajectory.

### Saved-only force decomposition

The separate `audit_thr2_force_terms.py` binds the completed capture and analyzes
**all 18 paths / 1,150 initial, callback and returned nodes**, without running an
optimizer. Only force-group measurement metadata changes, verified by exact
restoration to the original serialized force fields. Positive projection increases
the Cβ signed volume along its Cartesian unit gradient; this is not a work integral
or an L-BFGS direction.

At the unseeded whole start, OpenCL total Cβ projection is **−898.181
kJ/mol/nm**; the main contributions are harmonic angles **−1,414.322**, bonds
**+792.185**, nonbonded **−308.273** and periodic torsions **+19.672**. At the last
positive-volume callback, angles contribute **−1,576.020**, bonds **+915.105**,
and total **−918.015**. The Cα guard contributes zero there: its own centre
remains L and is above its activation boundary. This is not a Cβ chirality wall.

By contrast, the seed-1/2 whole starts have positive total Cβ projection
**+4,025.314 / +3,952.327**, with angle contributions **+4,469.238 /
+4,456.165**. Reference agrees on these signs. Identical heavy geometry and
System parameters mean the initial angle-projection differences involve the
different H geometry; they do not establish that a particular H or angle term
alone causes the entire later trajectory.

Measured force-group sum discrepancies across both platforms are at most
**0.000223 / 0.000264 kJ/mol/nm** for Cα/Cβ projections and **0.005157 kJ/mol**
for energy. Preserve those single-precision reduction differences rather than
asserting exact additive equality.

## Evidence and continuation boundary

- [Capture/replay audit](thr2-native-capture-evidence.json): aggregate geometry,
  exact-start/endpoint differences, gradient contracts and sign crossings.
- [Saved-force audit](thr2-force-terms-evidence.json): all-path bounds and
  initial/first-crossing-pair/last/returned force terms, no coordinate arrays.
- [Matched-start audit](thr2-start-lineage-evidence.json): complete atom/bond
  identity, byte-identical Systems/Integrators/parameters and heavy/H differences.
- Private captures and native callback/force arrays remain below the declared
  project-local `tmp/` root. The first audit is retained; the expanded independent
  audit was rerun and produced byte-identical private and public reports.
  Saved-force and matched-start reports also reproduce byte-identically.

The bounded sampling enrollment is finished, not extended until another inversion
appears. A new intervention to assign causality to a specific H placement/angle
would be a **separately frozen experimental stage**, not another repetition of
this cohort. No production fix, broad force-field conclusion, upstream issue or
recovered-original-trajectory claim follows from these results.

## Verification

Research checkpoints `0f6814d` (prospective protocol) and `2eb97eb` (results) are
pushed on `research/pdbfixer-reconstruction`.

- Full fast suite: **3,074 passed**, three optional Torch/Biotite skips,
  33 slow deselected. Two matched-start helper tests added after that collection
  are included in the final **17-passed focused suite**.
- The first full attempt retained **six failures / 3,066 passes** because the
  test shell lacked `dvbfixer` on `PATH`. With the interpreter's bin directory
  prepended, all six failing CLI cases and the full suite passed; no test skips,
  code changes or scientific-capture retries were used to hide that failure.
- Ruff (research/tests and production), strict-target mypy, generated CLI/GUI
  checks, agent-documentation checks and changed-document links passed.
- Expanded native audit, force decomposition and matched-start audit each
  reproduce byte-identically; original native artifact hashes remain unchanged.
- Public reports contain no coordinate arrays, native force/gradient arrays,
  checkpoints or RNG states. All scientific payloads remain private in-project.

Retained test logs: `tmp/pdbfixer-research/thr2-fast-suite-01.log` (initial PATH
failure), `thr2-fast-suite-02.log` and `thr2-fast-suite-02.xml` (corrected full run).
