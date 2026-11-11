# Mac reconstruction feasibility: execution log

Date: 2026-11-11. **RFD3 and Boltz-2 completed native MPS checkpoint smokes.
Independent audit passed this RFD3 case without repair; the Boltz proposal moves
observations and fails a junction check after source restoration. Protenix v1
private MPS sampling completed; its clash was repaired generated-only.** This is separate from earlier CUDA
evidence and Protpardelle MPS experiments; no broad model qualification is claimed.

Design: [Mac reconstruction portfolio handoff](../plans/mac-reconstruction-portfolios-handoff.md).
Evidence root: `tmp/mac-reconstruction-2026-11-11/`; owned runtime root:
`.artifacts/mac-reconstruction-2026-11-11/`. These local trees are not guaranteed
to transfer through Git. Existing engines/environments and historical results
were not overwritten.

## Current results (supersede the historical first-attempt pause below)

The user requested autonomous continuation. `resource-policy-amendment-v2.json`
replaces the overly conservative standalone swap-growth stop with `pressure-v2`:
8 GiB RSS, 6 GiB MPS driver, 1800 seconds, critical kernel pressure or <1 GiB
available sustained five seconds. Missing pressure telemetry fails closed. Swap
is retained telemetry, not a stop condition by itself. No MPS watermark override,
generic unsupported-operation CPU fallback or killing user applications occurred.

| Seed-7 arm | Wall time | Peak aggregate RSS | Independent result |
|---|---:|---:|---|
| RFD3 attempt02 | 50.05 s | 2,634,924,032 bytes | Exact 72-residue/574-atom inventory; observed drift after fit <4e-14 Å; generated held-out RMSD 1.012 Å; shared geometry and junction checks pass without refinement |
| Boltz-2 attempt02 | 37.59 s | 3,889,594,368 bytes | Exact inventory; native observed RMS/max drift 1.026/4.482 Å after observed-only fit; generated held-out RMSD 1.131 Å; source-restored A9–A10 C–N 1.482 Å fails existing junction connectivity gate |
| Protenix v1 attempt02 | 56.87 s | 2,532,851,712 bytes | 200 verified reinjection callbacks, exact overwrite each step; saved native arrays recovered without resampling after an extra-OXT materialization refusal; generated held-out RMSD 1.099 Å; LYS11–GLU34 clash before repair |

The audit fits only observed coordinates; held-out coordinates enter scoring
afterwards, never fitting, proposal or repair. Native, aligned and restored arrays,
unaltered original prediction files, masks and independent reports remain under
`rfd3-audit/`, `boltz-audit/` and the attempt02 directories. `aligned-native.pdb`
is a source-metadata materialization, not an untouched native model export: the
original export and full-precision native arrays remain authoritative for drift.
The evaluation mask uses the existing request builder but **does not claim another
backend executed the Protpardelle profile or produced a Protpardelle trace**.

RFD3 telemetry recorded 77,273 model-module calls with MPS outputs and peak MPS
driver allocation 2,512,093,184 bytes. Its requested schedule is 200 timesteps;
the observed diffusion-module call count is **199**, with 398 inner encoders for
two recycles. Do not relabel that as 200 measured updates. This tiny case is not
an efficacy benchmark or evidence for large/joint/multichain requests. A valid
junction is not rebuilt merely to exercise a repair path.

Boltz remains a template-conditioned proposal, not strict source-preserving
inpainting. Its native geometry pass cannot excuse observed drift, and restoring
the source cannot conceal the newly invalid junction. No seed or threshold was
changed to obtain these results.

Protenix assets completed exact pre-existing v1 checksum verification. The
maintained callback patch was applied **only to the owned official clone**, with
a patch receipt. Native Kalign 3.6.0 (`edd04a358dd8ea1e9379dd74053768d28c7b8613`)
passed a real two-sequence alignment. Private MPS reinjection passed a synthetic
rotated/translated control. Torch 2.7.1 emits a built-in MPS SVD CPU-fallback
warning even when generic fallback is disabled; the private implementation
therefore explicitly computes only the 3x3 covariance SVD on CPU with NumPy.
Coordinate transforms and model computation stay on MPS. This is declared mixed
rigid-fit implementation, not a fully MPS-native SVD claim. Existing CUDA code is
unchanged. OpenMM and MDAnalysis were added only to the new Protenix environment
to satisfy imports of the reused mapping/materialization/validation helpers.

Latest focused engineering verification: eight input/resource/audit tests passed.
The audit reference-coordinate control is explicitly not diffusion evidence.

### Protenix execution and common generated-only refinement

Attempt01 passed model inference but failed in upstream confidence reporting:
`torch.ones(..., dtype=float).to(mps)` requested unsupported float64. The separate
maintained `mps-confidence-dtype.patch` preserves the PAE tensor dtype and device;
it applies only to the owned clone. The failure and original attempt return code
remain visible. Attempt02 completed 200 steps and confidence, saving full native
arrays and per-step pre/post-overwrite drift, but its first PDB materialization
refused one unrequested upstream A72 OXT. A separate CPU feature-axis readback
verified the saved atom names/residue names, selected the existing canonical
requested-atom projection, and materialized **without resampling or reloading the
checkpoint**. The extra OXT stays in the original 575-atom arrays; the requested
projection contains 574 atoms. The original attempt02 exit remains failure, with
materialization-only recovery recorded separately. Subsequent entries include
the explicit projection and retain failures in `finally` cleanup.

The post-reinjection Protenix candidate passed inventory/source preservation but
had one severe LYS11 NZ–GLU34 OE1 overlap (0.99 Å). Its fitted/reference score is
after declared reinjection, not an untouched unconditioned proposal. The raw
pre-overwrite coordinates and per-step drift remain retained separately.

Boltz and Protenix were independently repaired with the same generated-only
Reference refiner: all 28 new atoms mobile, all observations fixed, 1000 iteration
cap, zero restarts, seed 7, existing geometry-based scoring, no hidden reference
coordinates. Both final assessments passed geometry and source custody. No
observed anchor was released; this is not anchor-limit calibration. Clean RFD3
junctions were not rebuilt. These are labeled derivatives, not retroactive raw
passes or broad comparative efficacy evidence.

### Prospective single-case six-seed development replication

`development-panel/plan.json` freezes 15 new serial jobs: seeds
42/101/102/103/104 × RFD3/Boltz/Protenix, with unchanged per-arm sampler settings,
input/code hashes, no retries and the same independent evaluation gates. The
three original seed-7 pilots remain retrospective controls with all technical
attempts retained. Reserved seeds 211/223/227/229/233 are untouched. This small
single-case replication is not the later multi-case conditioning or anchor
calibration matrix and includes no ranking optimization or public capability
change. Completion is established only by retained job results/completion index,
not directory existence.

All 15 jobs subsequently completed with independent audit readback. Including
original seed 7, source-restored raw geometry is RFD3 6/6, Boltz 0/6, Protenix
4/6; these numbers are **not an efficacy ranking**. Boltz's max observed drift
spans 3.877–6.957 Å. Protenix's two rejected raw cases are seeds 7 and 103 (clashes).

**RFD3's nominal six-seed ensemble is not six distinct stochastic proposals.**
All exported generated coordinates were identical. Two separately declared,
non-mutating model-entry controls at seeds 7/42 confirmed the actual final
consumer mask still has exactly 546 fixed / 28 diffused atoms and all sequence
fixed. CPU RNG states differ, but MPS RNG states and first noisy arrays are
identical; final full-precision arrays differ only by 9.54e-7 Å. The reset/seed
transport cause has not been established. Keep the original cohort unchanged
and do not infer structural stability from this lack of diversity.

A separately declared `private-mps-seed-at-model-v1` arm explicitly seeds MPS
immediately before model forward, after checkpoint/preprocessing. Seeds 7/42
now have different MPS RNG and first noisy arrays, with final-coordinate max
difference 0.354 Å. Raw seed 7 fails geometry and seed 42 passes; held-out RMSDs
are 1.033/1.023 Å. This is a **new seed contract**, not sampler/mask/weight changes
or a relabeling of the untouched upstream cohort. The prospectively allowed
101–104 extension proceeds only after these seed-transport controls passed.
Reserved seeds remain untouched. No raw failure is silently replaced by a retry.

The explicit-MPS-seed extension completed at all six development seeds, with six
distinct first-noisy array hashes. Raw geometry passes 4/6 (7/101 have LYS11
clashes); all six pass after the same generated-only repair policy, skipping
already-clean cases. Boltz and Protenix also pass 6/6 after their separately
retained source-restored/generated-only derivatives. This one tiny case does not
establish relative model quality, broad reconstruction readiness or calibrated
anchor motion. The nominal upstream RFD3 cohort remains unchanged and separately
flagged for failed seed transport.

### Experimental backbone/surface diagnostics on retained 1HE8

New private primitives `backbone_intersections.py` and `loop_surface_probe.py`
implement separate spatially indexed 3D backbone-segment distance alarms and an
explicit straight-chord/centroid-fan loop-surface probe. They do not implement
P3, repair topology, change old acceptance or claim an open-path knot invariant.
`screen_reconstruction_backbones.py` binds the source/request and honors shared
peptide-connectivity evidence, TER, chain case and insertion codes. Noncanonical
records and unconnected boundaries are reported. Explicit CONECT links are
retained; SSBOND and sidechain/disulfide surface topology remain untested.

On the historical selected patchwork, all five retained relaxation trials and
the final structure have zero ≤0.75 Å nonbonded backbone-segment alarms, but the
chosen artificial surfaces show the user-suspected generated pair:

- gap 2 (A436–456) surface is pierced by generated C492–N493;
- gap 3 (A489–497) surface is pierced by generated N446–CA446;
- gap 2 also has an observed N458–CA458 piercing.

These persist in the retained endpoint/trial structures. This supports examining
the visually identified pair, **not proof of impossible threading, absence of
chain passage along unsaved optimizer steps, or a qualified selection veto**.
An artificial closure can intersect legitimate paths; the observed A458 hit
illustrates why alarms must not become blanket fatal gates. No held-out false-
positive calibration or complete disulfide-topology analysis has executed. The
first surface-report serialization failed on a NumPy bool; its empty artifact
is preserved, and the corrected report has a separate attempt02 filename.

Fifteen focused segment/projection/surface controls pass: synthetic intersections,
projection-only crossings, parallel/endpoint cases, rigid invariance, chain case,
insertion codes, TER/connectivity, closure/end sensitivity and bounded-work refusal.

### Hash-bound visual constraints and private compatible-group search

`compatible_portfolio.py` adds atomic groups, manual locks/allow-lists, pairwise
pruning and mandatory complete-assembly assessment. Small constrained spaces are
enumerated; larger spaces use a deterministic width-4/64-evaluation search, with
homogeneous controls when permitted by the same manual mask. Ranking keeps
admission failures and actionable geometry counts lexicographically separate.
Joint groups cannot be split into donors from different runs. This pure planner
does not materialize coordinates, grant scientific acceptance or fabricate a
backend trace; its scientific adapter must verify custody and whole assemblies.

The new 1HE8 adapter verified the original five donor coordinates, trace/evidence
hashes and visual review bindings, and applied exactly the user's 36/50 approved
donors. Its constrained space is 75,000 combinations. It independently assessed
64 complete assemblies; the best still has 26 actionable findings and zero
admission failures. Optimality is not established, and this is not an accepted
assembly. The gap2/gap3 direct-clash pre-screen does not treat surface alarms as
qualified vetoes. All other pairs are checked by whole-assembly assessment rather
than silently assumed compatible. No old donor or prior selection was overwritten.

The selected assembly proceeds to one shared refinement with **all 1712 new atoms
mobile simultaneously**, no observed motion and no issue-local freezing. The
first launch failed before any refinement because the adapter initially copied
the source as `source.pdb` while its retained request expects `1HE8.pdb`; the
byte-identical request-name alias and correction receipt are separately retained.
This private experiment is not yet public portfolio/GUI integration. Six atomic-
group search tests and the existing portfolio regressions pass; the combined
focused suite initially had 66 passing tests; the latest combined portable
portfolio/diagnostic/input/mapping/refinement suite passes **79 tests**.

The plain all-generated CPU control completed but did **not** fix this assembly:
26 initial actionable findings became 27, with no observed-record admission
failures. The environment requested two CPU threads, but the shared refiner
explicitly supplies `Threads=1`; the actual control is single-thread CPU. It is
separate from the Reference attempt, not a replacement or a bit-identity
comparison. Reference subsequently completed with 26 remaining findings, also
rejected. A distinct existing-geometry-preconditioning control ran with
bonded/soft-steric/full-force-field iteration caps
100/75/75/1000, zero restarts and the same all-generated/fixed-observation mask.
This reuses the existing conditioning algorithm, not P3 or a new closure solver.

The preconditioned CPU derivative reduced the 26 initial findings to seven:
detailed Ramachandran findings at A349/A581, amide/junction findings at A895–896
and their aggregate gates. Observations remain exact; this is still rejected.
All three completed refinement endpoints have zero segment/surface alarms under
the experimental diagnostics, which does not excuse their geometry failures.
A separately declared four-case derivative search changes retained donors only
for residual-issue gaps 0/4/6 (three deterministic user-allowed recipes), plus one
explicitly **unreviewed** richer-context RFD3 joint group for gaps 2/3. Every case
still minimizes all 1712 generated atoms jointly, with observations fixed,
CPU Threads=1, seed 7, zero restarts and caps 100/75/75/1000. Results belong under
`1he8-shared-refinement-derivative-search/`; no completion is inferred from its
directory or partial files. No anchor has been released or calibrated here.

All four derivatives subsequently completed and passed independent readback of
exact observed records and global L-chirality. **None passed geometry.** The
three user-allowed donor swaps leave 33/85/23 actionable findings, with zero
segment alarms but 1/15/5 experimental surface piercings. The separately unreviewed
RFD3-context joint group leaves six findings, zero segment alarms and zero surface
piercings: detailed A896 backbone angle, A581 Ramachandran and A895–896 junction
C–N 1.5127 Å, plus aggregate gates. These are retained in
`independent-final-audits.json` and each derivative's independent surface report.
The baseline cc94-only conditioned assembly's seven findings and all four
derivative failures remain distinct; six findings is not an accepted assembly,
nor a demonstrated superiority result for the new group.

### Two new complete-reference pilots and common repair derivatives

Source/reference admission passed for explicit projections 1FR2 A3–85 with
hidden A24–31 and 1A3N A1–59 with hidden A30–37. Original bytes, unchanged
occupancies and source/native geometry reports are retained in
`reference-panel-preparation/`. These are cropped fragments, not their complete
native complexes or unseen-family qualification. RFD3/Boltz reader inventories
and Protenix hidden-template masks were checked before loading weights.
All six predeclared seed-7 jobs completed serially without a guard stop.

| Case | Backend | Raw held-out RMSD (Å) | Max native observed drift (Å) | Source-restored geometry |
|---|---|---:|---:|---|
| 1FR2 | RFD3 explicit MPS seed | 6.625 | <1e-12 | pass |
| 1FR2 | Boltz stock proposal | 1.196 | 11.945 | pass |
| 1FR2 | Protenix reinjection | 3.514 | <1e-12 | junction rejection |
| 1A3N | RFD3 explicit MPS seed | 0.966 | <1e-12 | clash rejection |
| 1A3N | Boltz stock proposal | 1.326 | 5.410 | junction rejection |
| 1A3N | Protenix reinjection | 1.072 | <1e-12 | pass |

Native drift is measured after an observed-only rigid fit, before restoration.
Protenix here means its post-reinjection output, not an unconditioned proposal.
Held-out coordinates enter scoring only, never fitting or refinement. A clean
1FR2 RFD3 geometry result at 6.625 Å illustrates why geometry acceptance is not
native accuracy; Boltz's restored pass does not establish hard conditioning.

Three raw failures received the same separately declared generated-only
Reference policy (seed 7, 1000 iteration cap, zero restarts, no preconditioning,
all observations fixed). All three passed final geometry and source custody:
1FR2 Protenix refined RMSD 3.556 Å; 1A3N RFD3 0.885 Å; 1A3N Boltz 0.896 Å.
Raw passes were not rebuilt and raw failures remain unchanged. Evidence is in
`reference-pilot-generated-only-repair/`, including independent audit readback.
This adds case evidence, not broad ranking, positive anchor-motion calibration
or public capability authorization.

The frozen Boltz CUDA callback patch was dry-checked against the owned MPS fork
and refused at `diffusionv2.py:299`; it was **not applied**. Fork-compatible strict
per-update conditioning requires its own maintained patch, exact feature-axis
audit and separate prospective arm. The stock proposal cohort stays unchanged.

The separate `deploy/boltz-2/mps-fork-per-step-callback.patch` subsequently
dry-checked and applied only to owned fork revision
`1431d573c3a02a2cb6b5c3b775e7128c7a98d539`. Its source snapshots, original hashes,
applied diff and prospective seed-7 pilot declaration are retained in
`boltz-strict-hook-source-review/`. Nine isolated real-sampler/fake-score and
padded-frame controls pass on CPU/MPS, including rejection of steering, invalid
shape/dtype, nonfinite returns and repeated callback steps; the disabled and
identity callback paths also yield identical seeded control outputs. These are mechanism
controls, not checkpoint acceptance.

The private `mps_research_smoke.py` entry validates real token traversal, names,
elements, atom-to-token mapping and unpadded prefix, repeating the feature audit
at the actual model input. It retains full native arrays and per-update drift;
only the 3x3 fit SVD is explicitly on CPU. OpenMM 8.6.1 was added to the newly
owned Boltz environment for canonical PDB materialization imports, not to an old
environment. The first 1UBQ seed-7 checkpoint pilot is separately declared and
guarded under pressure-v2. Its first attempt stopped before model loading because
Lightning's tensor-only loader refused OmegaConf checkpoint metadata; a second
attempt used explicit `weights_only=False` only after the pinned checkpoint
checksum passed. No global unsafe-load override was introduced.
**Fixation is post-update:** random augmentation and noise still move coordinates
inside each iteration, so this is not an assertion that every noisy model input
is byte-fixed or that Boltz has acquired a qualified inpainting capability.

The second attempt completed all 200 exact MPS overwrite callbacks but exited
through an overly strict returned-device/shape guard. Its original returned
prediction was not saved; the retained full-precision last pre-overwrite array
and declared exact float32 overwrite permit a separately labeled reconstruction
of the final callback state **without resampling**. This is not an untouched
returned-model export. That callback-state derivative independently passes
inventory, geometry and source custody (held-out RMSD 0.940 Å). The original
attempt exit remains failure. A third technical attempt now retains returned
device/shape/coordinates, checks MPS at model input and each callback, allows
Lightning's CPU result transport, and verifies the returned fixed coordinates.
It is not an independent new seed or a retroactive acceptance of attempt02.

Attempt03 completed in 28.41 seconds, peak aggregate RSS 3,459,072,000 bytes,
without a guard stop. Actual final-consumer tensors and all 200 callbacks were
MPS FP32; Lightning returned a CPU tensor of shape `(1, 576, 3)` (574 real atoms,
two padding slots). The returned real coordinates exactly match the explicitly
reconstructed attempt02 callback state (zero max full-precision difference).
Independent audit of the actual returned candidate passes inventory, fixed-source
records and raw geometry at 0.940 Å held-out RMSD, without repair. The extra
technical attempts do not increase the independent sample count beyond one.

Two further seed-7 fixation controls were prospectively declared for the same
admitted 1FR2/1A3N cases and original source masks, under
`boltz-reference-post-update-fixation/`. They retain the same FP32 sampler budget
and no steering; no reference coordinates enter model conditioning. They test
the coordinate intervention separately from the unchanged stock proposal arm,
not broad multi-seed efficacy or qualification.

Both reference fixation controls completed and independently passed raw geometry,
source custody and exact inventory, with <1.6e-14 Å fitted observed drift.
Held-out RMSDs are **4.132 Å on 1FR2** and **1.418 Å on 1A3N**, versus stock
proposal scores 1.196/1.326 Å before restoration. Thus the intervention removes
observed drift in these cases but does not establish improved accuracy; 1FR2 is
substantially worse in this single-seed comparison. The 1A3N stock junction failure
and repaired derivative remain separate from the new raw pass. No multi-seed
causal/efficacy claim follows from these two case controls.
The stock CLI and private checkpoint entry also differ in orchestration;
byte-equivalence of every intermediate feature/cache state was not established.
Treat these as source/mask/model/budget-matched mechanism pilots, not a proven
pure-intervention causal estimate. `boltz-fixation-comparison.json` retains this
limitation alongside the raw comparison values. Independent global L-chirality
checks passed for all three fixation candidates and all three new repaired
reference candidates.

### Existing MODELLER junction route: source review before new closure solvers

`src/dvbfixer/model/modeller_run.py::run_modeller` is not a drop-in repair for an
already complete diffusion proposal. It returns the input unchanged when its
alignment contains no gaps (lines 665–673). `_PinnedLoopModel` restricts only
`select_loop_atoms`, deliberately not initial automodel CG (lines 562–587), so
that production route does not promise byte-fixed observations during every
stage. Its fallback `Selection(self)` when no loops are found must not be reused
as permission to move an entire observed complex in this experiment.

Consequently no ordinary `model` invocation has been counted as an executed
junction-repair comparator. A prospective private control would need an explicit
generated-atom selection on the complete proposal, retained original/source
records, the existing chirality restraints and independent source-custody audit,
followed by the same all-generated shared refinement. That controlled adapter
and actual scientific comparison remain pending. No CCD/KIC implementation or
P3 work is justified merely by this source review.

A prospective private complete-proposal CG control was then declared on the
six-finding RFD3-context assembly, with all 1712 generated atoms selected and no
observed mobility. The first attempt failed before optimization: `complete_pdb`
renumbered the two chains globally and added one terminal OXT, so its raw identity
inventory did not match. Original failure and inventory forensics remain retained.
The second separately labeled attempt uses an explicit input-file-order residue
bijection per case-sensitive chain, checks every residue name and original atom
coordinate before optimization, and reports any added terminal OXT as a temporary
fixed topology atom excluded from output. It reuses the existing stereochemical
and Gaussian chirality infrastructure, not the production no-gap shortcut or a
sequence-based identity guess. Only after these checks may it run CG and the
same mandatory shared generated-only refinement. Completion and scientific
assessment are still required; no positive MODELLER result is inferred here.

Attempt02 subsequently executed CG followed by mandatory all-generated shared
refinement. The MODELLER intermediate has ten actionable findings; the shared
CPU result has five (detailed amides A895–896/A901–902 and Ramachandran A979,
plus aggregate gates). Independent final audit verifies exact observed records,
global L-chirality and zero segment/surface alarms. This remains a geometry
rejection, not a successful junction comparator or permission to waive amide
gates. The original six-finding pre-MODELLER derivative, ten-finding intermediate
and five-finding final are all retained separately. The temporary fixed terminal
OXT and explicit numbering bijection are recorded in `identity-transport.json`.
No observed motion, new closure algorithm or public capability was introduced.

### Bounded completion of the two reference development panels

`reference-development-extension/plan.json` prospectively freezes exactly
**30 serial calls**: seeds 42/101/102/103/104 × the same admitted 1FR2/1A3N
projections × RFD3/Boltz/Protenix. Original seed-7 outputs are reused, not
resampled. These complete the already approved six-development-seed panels;
validation seeds remain untouched. No anchor sweep, new crop, new reference
selection or geometry threshold is introduced. Commands, per-seed request/source
copies and code/input hashes are retained. The optional Boltz hook now exists in
the owned clone but is disabled in its stock CLI arm; the changed source diff is
explicitly retained, and the original pre-hook pilot remains unchanged.

The new orchestration has a no-checkpoint declaration test for distinct outputs,
correct seed arguments and Protenix request/source fidelity. The independent
summary script waits for the completion index and reports every declared seed in
execution/geometry denominators, including technical/audit failures. Heavy-atom
and backbone held-out RMSDs are scored independently over **all auditable outputs**,
not merely accepted or best seeds; missing output errors are not imputed as zero.
These experiments are still development on two cropped proteins, not a frozen
validation-family cohort or a calibrated positive anchor-motion envelope.

Latest combined portable test run: **82 passed**. Separate owned Boltz-runtime
mechanism controls: **nine passed**. Ruff, CLI-reference synchronization and
agent-document maps pass; the complete Python/GUI suites were not rerun and no
GUI/public capability changed. The extended panel runs asynchronously; partial
results are not completion evidence and are not polled for progress.

The extension subsequently completed **all 30 calls**, with no technical/audit
failure. Together with retained seed 7 this yields 36 auditable outputs. Both
RFD3 reference panels have six distinct first-noisy array hashes. The independent
all-seed summary is `reference-six-seed-summary-attempt02.json`:

| Case | Backend | Source-restored geometry passes / 6 | Native source + geometry passes / 6 | Median hidden heavy / backbone RMSD (Å) |
|---|---|---:|---:|---:|
| 1A3N | RFD3 | 4 | 4 | 0.973 / 0.384 |
| 1A3N | Boltz stock proposal | 1 | 0 | 1.435 / 0.880 |
| 1A3N | Protenix reinjection | 3 | 3 | 1.143 / 0.318 |
| 1FR2 | RFD3 | 6 | 6 | 6.130 / 5.582 |
| 1FR2 | Boltz stock proposal | 2 | 0 | 3.559 / 2.657 |
| 1FR2 | Protenix reinjection | 4 | 4 | 4.450 / 3.880 |

All medians include every auditable seed, not selected successes. “Native” for
Protenix means after declared reinjection; stock Boltz source preservation fails
even when its native geometry is clean. The first summary's field
`raw_source_and_geometry_passes` actually assessed restored PDBs; its values were
not changed, but a separately retained label amendment splits the two meanings.
The original JSON remains preserved; use attempt02 for unambiguous accounting.

A common Reference generated-only arm then repaired the 13 newly failed extension
outputs, reusing the three existing seed-7 repairs. All 16 repair derivatives
passed geometry, exact-source custody and independent global L-chirality; the
20 clean proposals were not rebuilt. Final source-restored/repair medians are:
1A3N RFD3/Boltz/Protenix **0.917/1.090/0.958 Å**, 1FR2
**6.130/3.562/4.486 Å**, with final geometry/custody 6/6 in every cell. Evidence is
`reference-extension-generated-only-repair/` and
`reference-six-seed-final-derivatives-summary.json`. Raw failures, stock-native
drift and clean-but-inaccurate 1FR2 proposals remain visible. This shows that
fixed-observation repair suffices for this small reference panel, **not** that
observed motion is calibrated, all large/joint requests work or any model has
general accuracy superiority. No positive motion envelope can be established
merely from these zero-motion successes; the 1HE8 five-finding rejection is not
retroactively accepted or used to relax thresholds.

### Next bounded controls: longer hidden masks and minimal anchor feasibility

Two 16-residue hidden intervals were declared before further sampling on the
same unmodified projections: 1FR2 A24–39 and 1A3N A30–45. Source/native-reference
admission and all three engines' input inventories/template masks pass. These
are changes in task length on the same proteins, **not additional families**.
Their seed-7-only pilot is bounded to six serial calls in
`longer-reference-pilot/`; reference coordinates remain scoring-only, and no
threshold, occupancy or source geometry was changed. Original eight-residue
panels remain immutable.

After the fixed-observation donor and MODELLER routes still rejected 1HE8, one
separate necessity-triggered **engineering feasibility** arm was declared in
`1he8-two-anchor-engineering-feasibility/`. It selects only observed direct gap
flanks A895/A902, each independently required to have a residual amide/junction
finding; the entire 1712-atom generated mask remains mobile. All other observed
records and all observed metadata stay fixed. The restraint references are the
original source, uniform k=10000 kJ/mol/nm², no B-factor weighting, CPU Threads=1,
seed 7, zero restarts and phase caps 100/75/75/1000.

The predeclared per-anchor envelope is `(BB RMSD, Cα, heavy max) =
(0.10, 0.10, 0.50) Å`, the conservative candidate engineering point from the
handoff—not an inferred physical uncertainty or approved production limit. No
limits will be loosened to accept this case. Geometry, outside-mask custody,
global L-chirality and each movement metric are reported separately; even an
engineering pass cannot establish calibrated motion or authorize public
publication. The RFD3 joint group is still not user-visually approved. Earlier
zero-motion results and historical movement thresholds remain untouched.

All six longer-mask seed-7 pilots completed and were audited. On 1FR2, all three
source-restored raw proposals pass geometry, with hidden-heavy RMSDs
RFD3/Boltz/Protenix **3.392/1.255/0.968 Å**. On 1A3N all three raw proposals fail
geometry; the unchanged common generated-only repair passes all three, with
final RMSDs **1.400/1.259/2.442 Å** (raw **1.434/3.317/2.507 Å**). All final
outputs independently pass global L-chirality. No observation moved during these
reference controls. Boltz's original max observed drift is 11.727/6.713 Å before
restoration, so its geometry results are not native source preservation.
Evidence: `longer-reference-final-summary.json` and the raw/repair directories.

Full-mask RMSDs for eight- and sixteen-residue tasks score different atoms. A
separate post-hoc check therefore fits both outputs on the same shared observed
atoms and scores only the original hidden core. For 1FR2 the short→long core
errors are RFD3 **6.625→4.696**, Boltz **1.336→1.333**, Protenix
**3.514→0.848 Å**; for 1A3N **0.966→1.012**, **1.462→2.355**,
**1.072→1.726 Å**. These are different conditioning masks/generated scopes, not
matched-noise or pure causal comparisons. They support mask sensitivity, not a
rule to hide measured residues or enlarge production gaps. Original frozen
requests remain unchanged; `nested-mask-common-core-comparison.json` retains
the common fit mask, atom counts and audit bindings.

The minimal 1HE8 anchor feasibility arm completed but is **rejected on both
geometry and movement**. It leaves the A901–902 amide defect (plus its aggregate
gate). Per-anchor `(BB RMSD, Cα, heavy max)` movements are A895
**(0.221, 0.043, 0.361) Å** and A902 **(0.211, 0.106, 0.513) Å**, violating the
predeclared `(0.10, 0.10, 0.50)` envelope. Independent audit confirms exact
outside-mask records/observed metadata, 18 explicitly movable observed atoms,
global L-chirality and zero segment/surface alarms. These do not excuse either
failure. No threshold is loosened, and the result is not a calibrated positive
motion limit or accepted prepared structure.

### Coupled 1HE8 context-group seed extension

Five additional RFD3 context-group development calls (42/101/102/103/104) were
declared to complement the retained seed-7 pair. Every coupled output owns gaps
2/3 atomically; all conditioning is original observed context, and the other
eight cc94 donors enter only after prediction. These are unreviewed controls,
not covered by old user visual allow-lists. No additional families, reference
coordinates, public traces or new refinement are supplied in this sampling arm.

The first series failed all five calls before checkpoint/model forward: broad
string-prefix substitution rewrote `…pair-context-input` as `…seed-N/run-input`.
No structures were generated. The original failures remain; a separate amendment
and attempt02 use semantic input/output/seed argument rewriting, backed by two
regression tests. The controller now stops at the first technical failure rather
than advancing through further seeds. The retry is technical, not five extra
scientific samples. Completion must be established by its retained index, not
partial directory contents.

Before inspecting the completed retry results, a separate bounded follow-up
budget was frozen in `1he8-context-refinement-budget-declaration.json`: at most
two of the five new groups may receive shared refinement. Independent native
and source custody, native observation drift ≤0.001 Å, complete-assembly
admission and global L-chirality are required. Rank complete raw assemblies by
actionable geometry count, then seed; do not select by confidence, reference
error or unqualified surface veto. All five scientific outcomes remain retained.
Every selected derivative minimizes all 1712 generated atoms jointly, all
observations fixed, CPU Threads=1, seed7, zero restarts, caps 100/75/75/1000.
This is a bounded group experiment, not exhaustive optimization over all donor
combinations. Seed7 remains the existing comparator and is not rerun. No result
is claimed until completion and independent final audit.

Latest combined portable regressions: **85 passed**; nine separate Boltz MPS/CPU
hook controls pass. New context argument and repair/ranking scripts pass Ruff;
the agent-document maps remain valid. Full Python/GUI qualification is unchanged.

The corrected coupled panel subsequently completed all five calls: raw complete-
assembly actionable findings are seeds **42/101/102/103/104 = 27/27/30/28/31**.
Each has a distinct first-noisy hash, global L-chirality and native observed drift
<1.3e-13 Å. All five raw full assemblies have zero segment/surface alarms under
the unqualified diagnostics. Per the unchanged declared count/seed rule, seeds
42 and 101 enter the two allowed all-generated shared refinements; no raw full
assembly is accepted merely because it was selected for repair. Complete
admission independently rematerializes both gaps together from native exports,
checks custody hashes, and preserves the same other eight donors.

A visualization-only six-model raw ensemble (retained seed7 plus these five)
is exported under `1he8-rfd3-joint-context-visualization/`. MODEL 1–6 correspond
to seeds 7/42/101/102/103/104. Each model includes both joint gaps plus the same
fixed observed context, 1337 atoms; the other eight generated gaps are omitted.
Noncontiguous source segments have explicit TER boundaries. Independent Gemmi
readback verifies every exported coordinate against its retained raw input.
This is **not** a structure acceptance/publication bundle or a user visual
approval. Both loops must be taken from the same MODEL if later manually
reviewed; do not independently mix members of these joint proposals. Original
cc94 review bindings and the source-label recovery for historical seed7 remain
separate and unchanged.

Latest combined portable regressions including the ensemble readback controls:
**87 passed**. Ruff passes for the coupled-panel argument, bounded refinement
and visualization scripts. Broader qualification and public capabilities remain
unchanged.

### Completed selected-group repairs and a small-ring diagnostic coverage gap

Both newly selected coupled-group repairs completed and independently failed
geometry: seeds 42/101 leave **21/23 actionable findings**, versus six in the
retained seed7 coupled comparator. Observations remain byte-exact and global
L-chirality passes; existing backbone/loop-surface screens show zero alarms.
Thus the lower raw count (27 versus seed7's 30) did not predict the better shared
repair endpoint. This bounded count-ranking heuristic is not qualified as a
general post-refinement optimizer. All rejected outcomes remain retained.

Residual localization finds the dominant new cluster between retained cc94 gaps
7/8, not the newly sampled gaps 2/3: TYR972 aromatic-ring contacts with SER1044,
including an approximately 3.04 Å SER CA–CB bond. Direct coordinate inspection
finds the six TYR ring edges approximately 1.74 Å in these failed endpoints
(the retained seed7 endpoint is approximately 1.41 Å). These are already
geometry-rejected structures, not hidden defects in accepted output.

An explicit topology-bound investigation verifies every ordered TYR972 ring
edge and the SER1044 CA–CB bond in the OpenMM topology, then detects the latter
piercing the ring's chosen surface in both failed endpoints. The retained seed7,
MODELLER-final and rejected two-anchor endpoints have no such piercing in that
declared ring. **The selected raw TYR ring is nonplanar beyond the conservative
0.10 Å diagnostic scope**, so the new ring probe reports it undefined; an earlier
arbitrary fan had zero hits, which must not be treated as proof of an unthreaded
initial state. These endpoints do not establish the actual unsaved optimization
trajectory or causal timing of any passage event. This finding is separate from
the historical gap2/gap3 loop-surface suspicion.

`covalent_ring_probe.py` adds a separate private small-ring surface diagnostic,
requiring closed covalent edges and tested path bonds rather than duplicating
residue/ring templates. Neighboring covalent paths are excluded; missing/invalid
topology, nonfinite coordinates and budgets fail closed, while nonplanar and
degenerate surfaces remain explicitly undefined. The engineering planarity
limit is not a protein-acceptance rule. Five synthetic controls cover crossings,
near contacts, bonded-neighbor exclusion, topology binding, rigid invariance,
chain case/insertion codes, undefined surfaces and bounded-work refusal.

`screen_covalent_rings.py` derives a deterministic RDKit graph cycle basis from
the supplied canonical-residue OpenMM topology (no bond perception, aromatic
assignment or force-field reimplementation), restricting scope to intra-residue
3–8-member cycles and canonical/peptide paths. Unsupported atoms/cycles and
nonpeptide inter-residue links are explicitly reported, not silently treated as
tested disulfide/macrocycle topology. Native 1UBQ is an initial no-piercing
engineering control; broader nine-case retained/native diagnostics are declared
under `small-ring-diagnostic-controls/`. This remains a nonfatal prototype,
without legitimate-knot, full-disulfide or minimization-trajectory qualification.
No P3/topology repair or new public fatal gate was implemented.

All nine small-ring diagnostic controls completed. The three native projections
contain 7/14/15 small cycles, zero piercings and three undefined proline surfaces
each. This is 27 defined native surfaces without hits, **not a general false-
positive-rate estimate** or legitimate-knot/disulfide qualification. Each full
1HE8 structure has 202 candidate cycles: the selected raw assembly has 50
undefined surfaces (46 PRO, three TRP, one TYR); refined endpoints have 48
(46 PRO, two TRP). Only the failed seeds42/101 endpoints show a defined-surface
piercing, one each, both the verified TYR972/SER1044 pair. Undefined surfaces are
not counted as unthreaded passes. Complete reports and method/input hashes remain
in `small-ring-diagnostic-controls/`; no structure or earlier acceptance policy
was changed by this diagnostic study.

The latest combined portable suite passes **94 tests**, with nine separate owned
Boltz-runtime hook controls also passing. Ruff and agent maps pass. All currently
declared bounded sampling/refinement/diagnostic jobs are complete; broader
conditioning, positive movement-envelope calibration, detector qualification
and any public CLI/Homology/GUI integration remain separate unfinished gates.
No geometry-accepted whole 1HE8 assembly has been obtained in this work, no
historical outcome was overwritten, and no commit/push was requested or made.
The six RFD3 joint alternatives remain explicitly unreviewed; they must not be
merged into the user's manually approved donor pool without new visual review.

### Actual missing 1HE8 loops: local RFD3 conditioning controls

Three prospectively declared seed-7 controls completed with the explicit MPS seed
contract, unchanged 200-timestep/two-recycle sampler, source-only observations and
exact consumer masks:

| Arm | Model residues | Fixed / generated atoms | Whole-assembly actionable findings | Backbone segment alarms |
|---|---:|---:|---:|---:|
| Joint A436–456/A489–497, contiguous A428–505 | 78 | 410 / 201 | 32 | 53 |
| Separate eight-flank proposals, assembled together | 37 + 25 | 141/137 + 126/64 | 32 | 49 |
| Joint pair plus established context-span union | 168 | 1136 / 201 | 30 | 0 |

These are actual unresolved residues with no held-out native coordinates. The
other eight gaps come from the same retained visual-constraint assembly **only
after sampling**, never as observed model conditioning; their original generated
ownership remains explicit. Every model-native observation fit has <1.3e-13 Å
max drift and every assembled raw candidate is globally L. All three whole
assemblies still fail geometry. Pair-attributed detailed findings are 5/5/4;
whole-complex totals also include unrelated retained gaps and aggregate gates.

The richer context includes original request gap2/gap3 context spans plus the
complete intervening span. Its ten noncontiguous source segments have explicit
`/0` model cuts. These are artificial model entities, **not preservation of the
native whole-chain topology**. The short independent arms have less intervening
context than the joint arm, so this is not a conditioning-matched causal estimate
of joint generation. The context control ran in 80.30 seconds, peak aggregate RSS
2,893,824,000 bytes, without a guard stop. Its source/frame audit is retained.

Short-flank proposals have real geometric near-crossings with omitted observed
residues (e.g. A379 versus generated A439), not just projected crossings. They
also pierce the artificial gap2 surface 7/5 times; the richer-context candidate
has no such surface hits. This supports context sensitivity in this **one case**,
not a validated topology guarantee, general success or a new fatal gate.

The richer-context materializer initially refused an upstream identity defect:
`set_indices()` stored `src_component` using chain-ID string width, collapsing
A1030–A1037 into `A103` in exported metadata. The new maintained source-component
width patch applies only to the owned RFD3 clone. A CPU-only pinned-consumer
readback regenerated all 168 identities, checked all 160 intact prior mappings
(the original 161-entry dictionary also contains the collapsed `A103` alias),
and verified identical built coordinates, masks, sequence and atom axis. The
saved native output was then materialized with an **explicit hash-bound mapping
recovery**, without resampling. Observed coordinates independently verify that
recovery. The original truncated metadata remains untouched; equivalence of all
downstream scientific features beyond the built input is not claimed. No silent
positional/sequence guess or public compatibility claim is made.

Portable local-projection/mapping controls and explicit isolated native MPS
controls were added. In particular, long source IDs and MPS seed-repeat/difference
tests pass in the owned RFD3 environment; private Protenix reinjection controls
pass in its own environment. Reference-derived fake predictions in unit tests are
not checkpoint evidence.

## Authorized first slice

The user's subsequent “приступай” authorizes starting staged execution, not an
unbounded sweep. `authorization.json` records a small seed-7 smoke, one ML process
at a time, 8 GiB aggregate process RSS / 6 GiB MPS driver guards, 1800-second
timeout, disk headroom, no automatic CPU fallback and no sequence uploads to MSA
services. Development seeds remain `7, 42, 101, 102, 103, 104`; validation seeds
remain reserved. No six-seed scientific cohort, conditioning sweep, GUI selector,
automatic compatibility search or anchor calibration has executed in this slice.

## Initial source/runtime review (historical G1 checkpoint)

| Route | Pinned source | What passed | Not established |
|---|---|---|---|
| RFD3 | `fnachon/foundry` `7a9c246fbf140a743d14a35248f5039e7723a23a` | Isolated Python 3.12 environment, CLI help, model-package import, finite basic MPS tensor operation; dependency check | Completed checkpoint forward or reconstruction quality |
| Boltz-2 | `Novel-Therapeutics/boltz-community` `1431d573c3a02a2cb6b5c3b775e7128c7a98d539` | Isolated environment, CLI help, model-package import, finite basic MPS tensor operation; dependency check | Completed checkpoint forward, template drift or strict source-preserving inpainting |
| Protenix v1 | official `bytedance/Protenix` `85767b811c40ed46e73a9b39519cf6bfca8701ba` | Separate reviewed dependency environment; v1 model class imports with `LAYERNORM_TYPE=torch`; finite basic MPS operation; callback patch dry-check applies | MPS inference entry point, real checkpoint execution and reinjection acceptance |

The initially considered `DeepFoldProtein/Protenix-MPS` fork at
`66e519b1dcee4bc85fa8c96d6f4d6c24153be7ea` reports version 0.7.3 and registers old
v0.5 Base checkpoints, not the requested v1. It was **not** installed/executed as
a substitute for v1. The official v1 clone matches the existing CUDA source pin,
but its runner still selects CPU without CUDA. New private MPS work is required;
the frozen CUDA adapter was not changed.

RFD3 and Boltz environments resolved Torch 2.15.0 / NumPy 2.4.6; Protenix uses
upstream Torch 2.7.1 / NumPy 2.4.1. Exact dependency freezes and `uv pip check`
receipts are retained. Protenix uses an explicitly derived source-PYTHONPATH
inference environment, not a claim that unconditional upstream CUDA dependency
metadata is compatible with macOS. CUDA kernels and unused scikit-learn-extra
are omitted with documented reasons. Missing imports are not hidden by broad
exception handling. An initial default FastLayerNorm import attempted CUDA
compilation; the corrected attempt uses upstream's explicit Torch selector.

License review retained source-license hashes: RFD3 BSD-3-Clause expressly covers
its referenced weights; Boltz code and weight card identify MIT; Protenix v1 is
Apache-2.0. Restricted Protenix v2 weights are excluded. This is an artifact-level
engineering review, not permission to redistribute an assembled environment.

## Prepared input and actual RFD3 consumer audit

The new smoke is an explicitly labeled **1UBQ A1–72 projection**, withholding
A10–13 coordinates. It is not the original full-entry/frozen request. Selected
source occupancies and identities are preserved. There are 546 observed heavy
atoms and 28 canonical generated atoms, with the complete 72-residue sequence.
Held-out reference coordinates are stored separately, not supplied to prediction.

For RFD3, unknown-coordinate residues are represented by sequence-derived
canonical atom identities with zero dummy coordinates. The actual pinned
`DesignInputSpecification.build()` consumer now confirms:

- 574 atoms; exactly 546 coordinate-fixed and 28 diffused;
- diffused residues exactly 10, 11, 12, 13;
- every atom has fixed sequence;
- every diffused coordinate is zero at the built input.

This admission check found two preparation problems **before sampling**:

1. High dummy atom serials followed by original serials failed Biotite's strictly
   increasing-serial requirement. A new adapter input serializes monotonically;
   source/reference coordinates and residue identities are unchanged.
2. Selecting only the fixed flanks with a contig string left all atoms fixed
   under this consumer's defaults. The new explicit dictionary includes
   `"A10-13": ""`, clearing coordinate fixation for the withheld region.
   A hypothetical all-fixed run would not be reconstruction and was refused.

All input attempts and failure receipts remain retained. The current input is
`smoke-input-03/`; the consumer proof is `rfd3-input-preflight.json` and
`rfd3-built-input-attempt04.npz`. Initial import used the wrong class name and
is retained as a harness error, not an upstream/model failure. AtomWorks reports
long bonds at dummy/observed boundaries; these are visible input diagnostics,
not grounds to claim a completed reconstruction or a geometry pass.

Boltz's prepared input supplies the observed-only structure as a template and
uses `msa: empty`, without contacting a remote MSA service. This is initially a
**template-conditioned proposal arm**, not established strict inpainting. The
actual upstream YAML/template parser passed after chemical asset provisioning:
72 target residues and 574 atoms. A detailed final feature-to-source atom mapping
and coordinate-restoration audit is not yet completed; parser success alone is
not that proof.

## First RFD3 checkpoint attempt: stopped by host-resource guard

`rfd3-smoke-declaration.json` froze seed 7, one candidate, 200 steps, two recycles,
eager Torch kernels and no compilation/CUDA graphs. The upstream 2,690,316,669-byte
checkpoint has locally computed SHA-256
`9b3f85923e0d51e9453e15cdd2f8c666e7ce096a60577f57d11bbc54ae6d67c1`.
The upstream registry supplies no expected checksum; do not call this an
independently upstream-digest-verified artifact.

The attempt stopped after 14.01 seconds because host swap grew 643,694,592 bytes,
exceeding the declared 512 MiB guard. Peak observed aggregate process RSS was
1,570,406,400 bytes. The retained device trace has **zero actual model MPS forward
calls**. Input prevalidation ran, but no completed generated candidate exists.

This is a **technical resource interruption before a model forward**, not a
scientific rejected seed, proof of unsupported MPS, or an established model memory
requirement. Other host activity and memory state were not controlled. The host
subsequently had about 4.68 GiB available and 5.99 GiB swap in use. No user
applications were terminated and no safety limits were weakened. Retain this
attempt; a future retry requires a new explicit attempt/resource declaration.

## Asset transfer correction and remaining work

The first Boltz checkpoint download was conservatively refused because the
downloader incorrectly used an HTTP Xet ETag as an expected file SHA-256. The
retained bytes have SHA-256
`090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`,
which matches both upstream LFS metadata and the existing
`deploy/boltz-2/README.md` pin. No unverified checkpoint was loaded.

A separate corrected transfer uses immutable Hugging Face revision
`6fdef46d763fee7fbb83ca5501ccceff43b85607` and the existing LFS SHA-256 pins for
structure weights, molecule data and affinity weights required by the stock CLI.
The original refused file/receipt remain preserved. Completion must be read from
the transfer result, not inferred from partial-file/directory existence.

The corrected transfer **completed** for all three pinned Boltz files with exact
expected byte lengths and SHA-256. The molecule archive was extracted only after
checking member paths/types and disk headroom, using Python's data extraction
filter: 45,228 members, 1,820,698,819 unpacked bytes. The first refused download
receipt remains unchanged; the corrected transfer is separately identified.

### First Boltz checkpoint attempt

`boltz-smoke-declaration.json` records seed 7, one sample, 200 diffusion steps,
one recycle, explicit MPS, Torch kernels and no remote MSA. Native CLI template
preprocessing completed and Lightning selected MPS. The resource guard stopped
the process after 17.12 seconds due to 544,407,552 bytes of host swap growth.
Peak observed aggregate RSS was 1,285,537,792 bytes; the device trace again has
**zero actual model MPS forward calls**. There is no completed candidate or
scientific quality result. A shutdown semaphore warning remains in the log.

The two interrupted attempts do not isolate the cause of swap growth or establish
the models' peak requirements. No memory-watermark override, application killing,
CPU fallback, automatic retry or weaker safety budget was used. Pause further
checkpoint launches until host headroom can be established. Protenix checkpoint/
assets and the private MPS inference entry point remain pending. Next work is
bounded G1/G2 completion, then independent source/sequence/geometry audits, not
automatic cohort expansion.

## Engineering verification

Two new input-preparation regressions check hidden-coordinate exclusion,
monotonic dummy serials, preserved observations and refusal to overwrite inputs.
Together with the focused existing joint/component tests: **22 passed**. Ruff
passed for the new preparation/guard scripts and tests. All three new environment
dependency checks passed. These are engineering checks, not checkpoint-backed
scientific acceptance. The complete Python/GUI suites were not rerun for this slice.
