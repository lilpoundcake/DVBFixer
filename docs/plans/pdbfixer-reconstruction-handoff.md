# PDBFixer reconstruction: research plan and live handoff

## Authorization and baseline

**O–H implementation authorization:** user approved a separate implementation
and qualification branch after `347efe6`. Work continues on
`fix/hydrogen-postflight-contract` under the
[frozen postflight proposal](hydrogen-postflight-contract.md). Native forces,
minimizer settings and existing regression assertions must remain unchanged.
This is not authorization to merge/release without qualification or to activate
the separate empirical-library/reconstruction research routes.

**Implementation progress:** candidate v1 (`94589d9`) passed the narrow 54-record
distance replay but failed its full suite (3,144 passed/two failed): source-inventory
binding and a real retained HG–OXT contact of 0.718741 Å despite valid O–H length.
All raw failure artifacts and code are now preserved in-project. Prospective v2
amendment `693f54b` adds the existing >1.0 Å nonparent-heavy contact predicate for
O-bound H, conservative radial admission and typed fail-closed/diagnostic handling.
All 67 v2 records are accounted for: original54 pass, eight safe contact-graph
outcomes, four explicit blocked-toy refusals, and one real failed SER state repaired
to O–H0.970000/HG–OXT1.788592 Å. Independent replay is exact; focused94 pass.
See [implementation and validation](hydrogen-postflight-contract.md). No force
field, native optimizer, existing SER assertion or historical evidence changed.
Candidate-v2 full fast validation is **3,174 passed**, three optional skips and
33 slow deselections; focused94, Ruff, strict mypy and generated/map checks pass.
The original SER test file's complete digest is unchanged. v1's failed source,
PDB/Dat, log and XML remain retained; the old lost 1.245813 Å trajectory is still
not recovered. Implementation is isolated on the fix branch, not merged/released.

**Real-input follow-up:** prospective `5d99286` ran the declared six-source ×
two-FF cohort once. Eleven preparations completed; all twelve saved postflight
states replay exactly, with no new refusals and heavy/valid-H/identity/bonds
preserved. Twelve previously unguarded stretched O–H were repaired, all on
forcefield-free fallback, not native AMBER/CHARMM. Complete-cohort acceptance is
still false: existing GLYCAM-named hinge + CHARMM falls back after NLN template
failure, creates duplicate UYB H1 before postflight and fails later naming.
Its partial PDB is invalid evidence of success and `.dat` is absent. No native
retry. Frozen auditor initially failed on repeated chain-ID blocks; separate
offline transport amendment `6080fe0` preserves exact nm/graph/residue ownership
without changing frozen artifacts. See the
[full results and remaining coverage gates](../research/oh-native-qualification.md).
Do not claim native CHARMM glycan or direct ligand-OH repair qualification.
Final engineering validation: 3,184 fast tests pass (three optional skips,
33 slow deselected), ten observer/transport tests pass, and independent offline
audit readback is byte-identical. Production science is unchanged from `3ffa1f4`.

**2026-11-12 continuation:** the user approved a separate bounded original-input
1CRN/THR2 native capture and instructed autonomous continuation until genuine
input is needed. The prospectively frozen three-run protocol and separate
saved-state OpenCL/Reference observations are documented in
[THR2 native capture](../research/thr2-native-capture.md). Historical findings and
uncaptured original trajectory remain unchanged; production integration and
upstream issue publication are not authorized by this approval.

**Bounded THR2 capture completed:** frozen protocol checkpoint `0f6814d` precedes
all science. All 3 native labels and 18 saved-call replays completed; unseeded
native THR2 Cβ changed **R→S in whole phase 1**, seed 1/2 remained R. Both saved
OpenCL and Reference controls invert between callbacks 2/3 for the unseeded
start, with ~−72.616 kJ/mol objective change. All 1,114 callbacks independently
match physical mobile gradients/objectives; **0/18 returned/last-callback
mismatches**, unlike the separate THR24 terminal leakage. No negative Cα found.
The three native whole starts have byte-identical Systems/Integrators and all
327 heavy coordinates; 315 H coordinates differ. Initial Cβ force projection
is −898 vs +4,025/+3,952 kJ/mol/nm, dominated by opposite angle projections.
Saved-only all18/all1150 force decomposition retains numerical group-sum errors.
This captures a fresh native event, not the lost original November 6 trajectory.
Details, evidence and the no-adaptive-extension boundary:
[THR2 native capture](../research/thr2-native-capture.md). Production unchanged.
Protocol/results checkpoints `0f6814d` / `2eb97eb` are pushed. Correct-PATH full
fast suite passed **3,074**, three optional skips and 33 slow deselections; final
focused suite passed **17**, including two later-added lineage controls. Retain
the initial six missing-CLI PATH failures, not just the corrected pass. Ruff,
strict mypy, generated specs, agent docs and independent byte-identical saved-only
audits passed. Capture enrollment is complete; no further sampling is implied.

**Next separately frozen mechanism stage:** after the user's renewed instruction
to continue, saved-only per-angle attribution on all three THR2 starts identifies
HB-associated angle projection −2,729 vs +2,979/+2,980 kJ/mol/nm, with exact same
heavy-only +909 contribution. Independent analytic/native force agreement is
within 4.1e-12. The prospective [H-swap protocol](../research/thr2-hswap.md) enrolls
all six ordered donor/recipient pairs × HB/complementary-H masks plus three
unchanged starts, once in each numerical mode (30 paths/60 native calls). No new
H generation, seeds, physics, repair or production behavior; retain all outcomes.

**H-swap completed:** protocol `de12422` precedes all 30 paths/60 native phases.
All 4,097 callbacks independently match native physical gradients/objectives;
zero D-Cα/terminal mismatches, three exact old Reference controls. Outcomes follow
HB source in the full 3×3 HB/remainder grid in both numerical modes: unseeded HB
→ S, either seeded HB → R. The 30 paths contain **nine unique starts**, not 30
independent samples; reciprocal masks duplicate six hybrids and Reference phase
endpoints match exactly. Maximum OpenCL duplicate endpoint delta is 0.095902 Å.
Saved RNG reconstruction accounts for all315 H/all945 draws per fresh capture:
the unseeded perturbation already places HB in a four-ligand-sign conflict before
H minimization (−1.571 Å³ vs +0.198 heavy-centre volume), which persists through
H-only/repair before full heavy inversion. This is finite-context mechanism
evidence, not broad qualification, a production fix or the lost original run.

**Current genuine access blocker:** full fast suite retained **3,082 passed / one
failed** (three optional skips, 33 slow deselected). Unchanged SER126 HG regression
emitted OG–HG 1.245813 Å vs existing <1.2 assertion. Do not rerun until green,
weaken thresholds, skip or patch production. Log/XML are retained in-project;
read access to the original failing pytest workdir was denied as external.
Ask for access to the exact directory documented in
[H-swap continuation](../research/thr2-hswap.md) before copying/inspecting the raw
failure. Its cause is not yet established. Do not bypass the sandbox restriction.

**Access/loss update and renewed continuation:** user authorized read/copy, but
the tool still denied access and their own `cp` reported source path missing.
Treat the original failed SER-HG structure as unavailable; log/XML remain, and
fresh outputs must not replace it. The prospective
[SER-HG capture protocol](../research/ser-hg-capture.md) uses the same tracked
synthetic input in six fixed fresh preparations (unseeded + H seeds1–5), with all
H/minimizer/production-repair boundaries retained inside the project. Production
and assertions unchanged. A deterministic graph at the logged 1.245813 Å O–H
length shows the existing detector/test contract gap (<1.2 vs distance-only
1.552 Å), but does not establish missing predicates or the lost-run cause.

**SER-HG stage completed, protocol `89aaf27`:** all six fixed fresh prepares
completed (automatic OpenCL) and final OG–HG is 0.970200 Å in all six. The logged
original 1.245813 Å endpoint was **not reproduced**; do not add seeds or erase the
original 3,082-pass/one-failure full-suite result. All native intermediates instead
have HG near OXT (0.000371–0.091563 Å) and OG–HG 1.669–1.748 Å; actual production
repair catches HG by parent distance and restores it. All artifacts are private
in-project from the start. Static detector/test gap is deterministic.
All six native Systems are identical: correct OG–HG bond (0.96 Å), mobile HG,
pinned OG/OXT, graph distance5, no pair exception, HG charge +0.4474/OXT −0.8132,
zero mixed LJ epsilon. Saved-only all24-state force audit supports local Coulomb
collapse without an index swap; none of the returned H states passes the native
gradient criterion. Exact termination reason/terminal leakage remains unobserved.
Do not perpetuate the source-comment index-swap explanation as proven for these
fresh Systems. Focused research tests: 30 passed; production and thresholds
unchanged. [SER-HG results and implementation gate](../research/ser-hg-capture.md)
define the next genuine decision: authorize a separate behavior-changing H
postcondition/repair-contract implementation and qualification branch, or retain
research-only scope. No general cutoff or native FF change is qualified here.

User approved source research, controlled experiments and evaluation of a custom
geometry-constrained atom rebuilder on 2026-11-02. Branch:
`research/pdbfixer-reconstruction`, based on local `main` at `4a328a6`.
Production backends, thresholds and safeguards must remain unchanged.

Routine research-note/evidence/handoff maintenance is part of this ongoing study;
avoid repeated approval questions for each checkpoint. Discuss substantial new
scope or production changes separately.

## Current authorization: installed PyRosetta for research only (2026-11-03)

The user explicitly authorizes using the rotamer library already available in the
local PyRosetta installation for this research, **not production**. Continue local
experiments through PyRosetta under the applicable installed licenses; do not
block that research merely because the original Dunbrack server is unreachable.
This supersedes earlier statements that an upstream archive is needed before
the empirical-library research arm can proceed. It does not change prior results.

Keep bundled tables out of tracked files, packages and public artifacts; do not
redistribute or silently promote them into a production dependency. Record the
installed PyRosetta revision, selected library/configuration and relevant hashes
for reproducibility. Rosetta-bundled data are not automatically identical to the
original upstream release. User authorization is a research-scope decision, not
a replacement for third-party license terms or proof of redistribution rights.

**Mandatory production gate:** the user will try to obtain the original tables
through the appropriate license. Before ANY production integration, public
capability, packaging or redistribution of this rotamer-backed route, stop and
explicitly ask the user for the original tables, source/version provenance and
documented license/permissions covering the intended use and distribution.
Verify those materials and obtain explicit production approval. Do not assume
they have arrived, carry over research-only permission, or relax the scientific
qualification gates. If the materials are absent or unsuitable, production stays
blocked; local authorized research need not wait for them.

### Official-source retry (2026-11-03, 16:54 UTC)

The official host is reachable again: overview, license and example download
pages returned HTTP 200. Saved source HTML, retrieval URLs/times and SHA-256
records: `tmp/pdbfixer-research/dunbrack-official-retry-03/provenance.json`.
No original library archive has been downloaded and no application submitted.

- Overview: <https://dunbrack.fccc.edu/lab/bbdep2010> explicitly states CC BY 4.0
  since July 25, 2019, with free academic/industrial use, modification and
  redistribution with attribution.
- Linked terms/application: <https://dunbrack.fccc.edu/lab/bbdep2010_lic> instead
  contains an Open Data Commons Attribution text. It says a download link is
  emailed only after the application is reviewed, with real applicant details.
- Example download links lead to that application; the old
  `/bbdep2010/license/index.html` redirects there too. Do not confuse the example
  page or underlying training dataset with a downloaded rotamer-library archive.

Record the discrepancy, rather than silently choosing a license or asserting
that a paid agreement is necessary. The user's original-table acquisition and
the mandatory production approval/license check remain pending. Local installed
PyRosetta research remains authorized and does not depend on that application.

## Latest checkpoint: exact Reference replay identifies rejected-trial Context leakage

Snapshot-force attribution was committed/pushed as `f82eabc`. The
[whole native-path replay](../research/carbon-H-whole-path.md) and
[independent evidence](../research/carbon-H-whole-path-evidence.json) retain all 24
paths/48 phases and 3,991 callbacks: original two seeds/three conditioning states/
four arms. All starts and both original nm endpoints exact. Native physical
gradient/objective closure at every callback, max force-component error 0.0.
No new seed/force/H/phase/iteration change, no retry/early exit, no PDB or repair.
Max callback count 188 < frozen 1,000; these are computational replays, not samples.

17/48 returned Contexts differ from last reported iterate and all 17 fail the
source-defined native gradient criterion. Five returned energies higher, 12 lower;
do not equate a lower-energy rejected trial with convergence. All three negative
CA boundary snapshots have positive last reported CA:
seed 1/beta-only/methyl phase 1 +0.021875→−0.002076 Å³, ΔE +357.415 kJ/mol;
seed 2/raw/CH2 phase 2 +0.020362→−0.006118, ΔE +943.606;
seed 2/raw/combined phase 2 +0.017283→−0.006974, ΔE +1093.030. Last reported states
also not converged; never silently restore/publish them as qualified candidates.

[Pinned source metadata](../research/openmm-minimizer-source-36a30cb.json): runtime
declares 8.5.2.dev-36a30cb, revision 36a30cbca54e727b216b606f3c011b67201eb8b4,
Darwin arm64. Exact upstream C++/L-BFGS files with notices private/checksummed;
public URLs/hashes only. This is NOT attestation of an unmodified compiled binary.
Reference evaluate() writes trial coordinates into Context; L-BFGS failed search
restores its own x/g buffers; Reference ignores return code and does not restore
Context to accepted x. Unconstrained check then exits. Source-consistent rejected-
trial leak explains these controlled returns. Exact native error enum unobserved;
do NOT assert MAXIMUMLINESEARCH specifically or extrapolate to other builds/platforms.

Preserve adverse accepted states too: 16 target-CA negative measurements in eight
phase-1 paths; 137 target-beta negative measurements in 16 phase-1 paths (raw 8,
beta-only 8; full-generated 0). Min beta volume −19.101650 Å³ is badly distorted,
not 137 valid independent D residues or minima. All saved beta endpoints positive.
Thus not every adverse sign is a terminal leak and no finite CA guard protects CB.

Plan SHA `748b13a82a83656dd99677084ac19b1b08e490c69f066adf8c5a0dfa8e0e8bea`;
private `carbon-H-whole-path-{plan,run,audit}-01*`, independent audit02.
Eleven focused replay/audit/source-provenance tests and Ruff pass. Full fast suite
through this checkpoint: **3,066 passed**, three optional skips/33 slow deselected.
Log: `tmp/pdbfixer-research/carbon-H-whole-path-fast-suite-01.log`.
Production/defaults/thresholds unchanged.

**Authorization/access boundary:** this fixed Reference replay question is
complete. Original unforced 1CRN THR2 pre-min System/State/checkpoint was never
retained (original comparator folder has input/prepared/output PDBs/logs only).
For that ORIGINAL cause, ask user for an original failing checkpoint or explicit
approval of a separately declared bounded fresh capture on original input/platform.
Do NOT silently expand seeds/platforms until inversion, claim Reference proof
transfers to the uncaptured original event, change production or publish an
upstream issue without permission. Empirical-table production provenance/license
approval remains separate and closed; no original tables presumed received.

## Preceding checkpoint: native guard is present; adverse snapshots are not stationary minima

Complete whole results were committed/pushed as `cac26cc`. The
[saved-only attribution](../research/carbon-H-whole-attribution.md) and
[independent evidence](../research/carbon-H-whole-attribution-evidence.json)
retain all 24 THR24 paths/72 states/216 focal measurements: both original seeds,
raw/beta-only/full-generated conditioning, all four arms, target CA/beta and
context GLN16 CA. No minimization/integration/new H/force ablation/seed or phase
search. Fresh original and group-cloned Reference readback with exact group
restoration. All native energy/force sums close; max absolute force/projection
error 1.907e−6 kJ/mol/nm, energy 1.490e−8 kJ/mol; independent finite-difference
Cartesian gradient error 9.154e−11 nm². Complete force arrays private/checksummed.

All 144 CA measurements have the native guard, nine active; analytic guard
derivative matches observed projection. Raw/beta-only target vref 0.025934 Å³,
qmin 0.2 => activation floor 0.005186748 Å³, NOT 0.5 Å³ diagnostic threshold.
Full-generated target floor 0.511023 and GLN16 floor 0.505681 Å³. Frozen parameters,
not proposed settings. No missing guard or unit/sign mismatch established.

Seed 2/raw negative phase-2 target CA at −0.006118/−0.006974 Å³ (CH2/combined) has
positive guard projections **+11,554,131.899/+12,438,578.149 kJ/mol/nm**; native
totals +11,553,861.831/+12,438,326.631. Selected whole Systems have no constraints
or frozen particles. Huge nonzero forces mean NOT converged stationary D minima;
do not claim FF intrinsically prefers D here. Does not determine iteration cap vs
line-search failure or an unrecorded crossing. No extra relaxation is attempted.
Seed 1/beta-only/methyl-zero also has negative phase-1 CA −0.002076 with restoring
guard +7,205,585.536, then positive-flat final +0.007662 Å³. All 408 whole states
retain exactly these three negative-CA snapshots and no negative THR-beta.

At raw seed 2 final snapshots, angle projections ~−1109/−1163/−1130/−1135,
torsion +817/+840/+829/+850, nonbonded ~−170–173; restraint ~+159–163. Unit-volume
gradient projection is NOT an actual quasi-Newton step, integrated work, pair
causality or recovered trajectory. Original unforced THR2 event stays separate.

Plan SHA `81a5588e2cc1bac6471841cc31f38b2e9d6b670e57ed3580ea02e0ffbaf78d86`;
private `carbon-H-whole-attribution-{plan,run,audit}-01*`, independent audit02.
Nine new focused tests plus seven existing volume/guard tests/Ruff pass;
independent publications byte-identical. Full fast suite through attribution:
3,055 passed, three optional skips/33 slow deselected. Retained log:
`tmp/pdbfixer-research/carbon-H-whole-attribution-fast-suite-01.log`.

**Resume:** separately freeze observational native-callback replay of the same 24
paths, requiring exact starts and both frozen endpoints. Original two Reference/
tolerance 10/max 1,000 phases/tenfold positional-weight reduction only; no new seed,
force/H, reflection, retry or budget/iteration tuning. Reporter records original
coordinates/gradients/objective and never requests early exit. Computational
replays, not new samples. Establish exact transport before termination/crossing
claims; production/defaults/thresholds and empirical-table production/license
approval remain unchanged. Continue until genuine user input is required.

## Preceding checkpoint: complete whole factorial gains are limited and retain negative stress signs

Baseline/auditor checkpoint was committed/pushed as `12a3572`. The
[complete whole results](../research/carbon-H-whole.md#complete-whole-stage-results)
and [independent evidence](../research/carbon-H-whole-evidence.json) retain all 160
entries: 136 completed, 24 inherited source refusals; no retries/blockers/timeouts/
budget nonexecutions. All 408 states audited. All 34 unchanged native starts/both
endpoints and 14 known historical whole intersections exactly reproduce, all phase
Systems/anchors unchanged. Independent publications byte-identical (audit SHA
`a92343851a6996f17e6498ad9b367391d95c32f895ab4c44e09e852c02ed6101`).

On 24 1SHG stress contexts/arm, actual-frame raw/rounded masked heavy gate passes
unchanged/methyl/CH2/combined **10/11/10/13**. Final half-Å³ CA counts **8/7/9/5**,
negative CA **0/0/1/1**, all-acute HA **1/3/2/0**; no final THR-beta shape flags.
Earlier THR24-only methyl **10/24 /11/24 /12/24** remains separate unchanged
historical evidence, not replaced by this broader factorial. Combined gains three
gate passes without losing old passes but still fails 11/24; no general repair or
native-accuracy qualification. CH2-only loses THR24/seed 2/full-generated because
GLN16 globally flattens to 0.483534 Å³ even though its masked target passes.

Both negative endpoints are THR24 CA in THR24/seed 2/raw. Initial already-collapsed
volume 0.025934 Å³ in all arms; phase-1 volumes 0.903253/0.848707/0.910136/0.852397,
phase 2 **+0.009691/+0.017269/−0.006118/−0.006974 Å³**. Rounded CH2/combined remain
negative −0.004679/−0.005540. Preserve signs, but these are all near-planar stress
outcomes, not normal L→well-formed D or original unforced THR2 Cβ causality.
These are pre-emergency native phases, no public candidate/PDB and production
zero-D safeguards unchanged. Combined retains other CA volumes 0.004414/0.005554
at THR24 and ~0.0652 at THR32. Final-zero inward HA is not heavy acceptance.

The worst H-stage combined CG2–HG23 bond 1.855795 Å (THR24/seed 1/raw) is corrected
BEFORE whole minimization by the unchanged production repair bridge: exactly
HG21/HG22/HG23 change; bond 1.090000 Å at whole start, 1.079736/1.083349 after phases.
Do NOT credit whole minimization for that bridge correction or erase H-stage harm.
Heavy coordinates are not repaired by transport and remain badly flattened.

Ten prepared contexts/arm (five proteins/two seeds) have no final negative/half-Å³/
all-acute findings; masked gate N/A, not zero passes. Minimum final CA
0.500677/2.157882/2.156509/2.157992 Å³. Known ASP49/1TIM positives are exact
historical intersections, not independent samples. All original-heavy movement,
raw/rounded continuous volumes and H strain retained. Max stress H-angle
42.240/47.023/44.166/21.130°, prepared 28.629/11.852/11.915/11.674°. Observations
move; no immutable-repair/native-RMSD/safety certificate. Full fast suite: 3,046
passes, three optional skips/33 slow deselected; 12 focused whole tests/Ruff pass.

**Resume:** freeze separately saved-only force attribution on ALL THR24 source
conditioning contexts, both old seeds × raw/beta-only/full-generated × four arms
and all three saved states (24 paths/72 nodes). Original phase Systems and native
Reference force/energy readback; independent signed-volume gradient checks and
cloned-force grouping with exact parameter restoration. No new minimization/H,
force ablation, strength/seed/phase tuning or favorable-endpoint selection.
Snapshot forces alone cannot recover an unrecorded crossing or unique cause.
Original THR2 cause and empirical-table production/license approval stay open;
continue autonomously until genuine user access, decision or authorization needed.

## Preceding checkpoint: all 34 independent whole baselines audited

The all-context whole protocol was committed/pushed as `581de3d`. The
[baseline checkpoint](../research/carbon-H-whole.md#independent-baseline-checkpoint)
and [complete baseline evidence](../research/carbon-H-whole-baseline-evidence.json)
retain all 40 entries: 34 completed, six inherited source refusals, no retries/
blockers/timeouts/budget nonexecutions. All original nm starts and all native
phase force/anchor parameters agree with the source captures. All 24 legacy
historical initial/phase1/phase2 Å snapshots reproduce exactly at their original
precision; no old nm arrays are fabricated. Independent baseline publications
regenerate byte-identically. The earlier private baseline audit remains retained;
metadata-only addition of auditor hash does not resample or rewrite coordinates.

All 136 whole continuations are now running under unchanged frozen plan
`6dac11974efad016509db0d3482b7e3f0bcbe4f46072d3e383f72bfa858a5a11`.
Private `carbon-H-whole-baseline-01/`, baseline-audit01/02 and `carbon-H-whole-run-01/`.
Continuation outcomes NOT evaluated at this checkpoint. Saved-only auditor is
prepared for exact 34 native nm controls, 14 historical whole intersections,
existing masked stress-target gates, continuous raw/rounded global geometry,
all H findings and original-heavy movement. Prepared-view masked gate is N/A,
not a false zero-pass scientific result. Twelve focused worker/readback tests
pass, including no-minimization audit and unsafe recovery/identity/frame refusals.
Full fast suite through this implementation: **3,046 passed**, three optional
Torch/Biotite skips and 33 slow deselected; retained log
`tmp/pdbfixer-research/carbon-H-whole-fast-suite-01.log`.

**Resume:** wait for all-context whole run completion, then independently audit
every arm and source refusal. Require exact unchanged starts/both new native nm
endpoints, retain every adverse state and old Å precision distinction. No new
seed/force/phase/retry, emergency reflection, observation restore or public PDB.
Production/defaults/thresholds unchanged. Original THR2 cause, broad qualification
and empirical-table production/license approval remain open; continue autonomously.

## Preceding checkpoint: all-context whole protocol and exact transport admitted

Complete H-stage factorial results were committed/pushed as `9a192b9`. The
[whole protocol](../research/carbon-H-whole.md) and
[preflight](../research/carbon-H-whole-preflight.json) retain all 40 original
contexts, 34 independent native baselines and all 136 continuations plus 24
inherited refused arm entries. All 34 original native whole starts preflight
exactly in nm after unchanged H-label repair/object-index transport; mapped bond
graphs agree. No whole native-phase outcomes evaluated at this checkpoint.

Original full Systems/States/integrators/anchors, Reference/no-solvent/10 kJ/mol/nm/
max 1,000 iterations per phase; only original phase-2 positional weights divided ten.
Baseline native Simulation calls precede corresponding direct continuations.
Require all 34 unchanged starts and both independent native nm endpoints exact,
all phase Systems equal. Old legacy historical endpoints are Å-only snapshots:
compare actual new baseline nm→Å output honestly in Å, never divide old Å back
and call it original nm. Existing masked target checks remain separate from
global volume descriptors; all original-heavy movement/H strain/adverse sites
and source strata are required, no success-selected subset.

Whole plan SHA `6dac11974efad016509db0d3482b7e3f0bcbe4f46072d3e383f72bfa858a5a11`;
private `carbon-H-whole-plan-01.json`. Budgets 1,800 s/cell, 7,200 s baseline batch,
21,600 s continuation batch, no retry. Four focused baseline/precision/transport
tests and Ruff pass. Production/defaults/safety net unchanged, no public PDB.

**Resume:** execute 34 independent native baselines first; any blocker prevents
arm continuations. Then execute/audit ALL 136 arms and retain 24 refusals. Do not
turn narrow H-stage positives into heavy-geometry qualification. This remains two
phases BEFORE emergency repair, no reflection/observation restore. Original
unforced THR2 cause and empirical-table production/license approval remain open;
continue until genuine user access, decision or authorization is needed.

## Preceding checkpoint: complete carbon-H factorial retains mixed H-stage results

Admission/implementation was committed/pushed as `e316ce3`. The
[fixed factorial](../research/carbon-H-factorial.md) and
[independent evidence](../research/carbon-H-factorial-evidence.json) retain all 160 entries:
136 completed, 24 inherited source refusals; no retries/blockers/timeouts/budget
nonexecutions. All 34 unchanged native nm endpoints and all 14 known historical
intersections reproduce exactly. All 6,800 native mobile gradients independently
match (max component error 0.0); 7,072 initial/callback/final nodes retain H-bond
and H-angle strain at every node, complete HA findings and fixed observations.
Two independent private audits agree byte-for-byte.

Keep strata separate. Across 24 1SHG stress contexts per arm, any all-acute callback
cells are unchanged/methyl/CH2/combined **3/4/4/2**, final all-acute HA counts
**1/2/2/0**. Across ten prepared-view contexts per arm (five proteins/two seeds),
both descriptors are **1/0/0/0**. All matched denominators complete; six source/seed
refusals stay in all four arms. Reused positive cases are NOT independent samples.

Stress methyl-only creates a new LYS18 final inward HA in THR24/seed1/beta-only.
Known GLN16 final failure in THR24/seed1/full-generated persists under either
individual construction; CH2-only also yields final GLN16 in seed2/full-generated.
Combined has zero final all-acute H but retains THR32 seed2 raw/beta-only transient
findings at 10/11 callbacks. Zero endpoint count is NOT universal path repair.

Severe adverse combined bond: THR24 CG2–HG23 length **1.855795 Å**, deviation
**0.765795 Å**, THR24/seed1/raw. Matched unchanged maximum deviation **0.090327 Å**.
That fixed heavy input already has Cα volume 0.025934 Å³, explicitly a poor-heavy
stress context, NOT healthy production admission. Methyl-only worst deviation
0.726748 Å (THR24/seed2/raw). Maximum final H-angle deviations are
74.305/71.818/74.336/61.842° in stress and 71.533/72.281/71.256/43.884° in prepared
views. Final <1 Å nonbonded-H descriptor counts are zero but miss stretched bonds.
No geometry/convergence/whole-system safety qualification follows.

Execution plan SHA
`d136dbceb4016651f7f76e665ed18a33ab2e60ef7eade97f7e1e9670fc017b48`;
private `carbon-H-factorial-{plan,run,audit}-01*`, second independent audit02.
Audit SHA `1ecf63570b9f00e0f7f491a6fc8a08e14588d8048762ad5dffabd9a2fe0f943b`.
Twenty-nine focused tests and targeted Ruff pass. Full fast suite through initial
audit tests: 3,032 passed, three optional skips/33 slow deselected; later kept-H
and published-replay-binding tests separately pass. Production unchanged.

**Resume:** separately freeze ALL 136 whole continuations and all 24 inherited
refusals, not combined-only or success-selected endpoints. Establish independent
native baselines from original full System/State/integrator; exact unchanged nm
starts and both new native endpoints required. Historical legacy Å-only endpoints
must be compared honestly at Å level, never divided back and called original nm.
Keep all stress/prepared strata, original-heavy movement, continuous raw/rounded
volumes, H strain, masked heavy-target checks where already defined and every
adverse site. Original two Reference/no-solvent/10 kJ/mol/nm/max1000 phases and
tenfold native positional-weight reduction only, BEFORE emergency recovery.
No new force/seed/phase, accepted candidate/PDB or production/default change.
Original THR2 cause and empirical-table production/license approval stay open;
continue autonomously until genuine user access, decision or authorization needed.

## Preceding checkpoint: complete carbon-H factorial admission, before execution

The downstream checkpoint is committed/pushed as `f32595f`; metadata-only
[factorial enrollment](../research/carbon-H-factorial-enrollment.json) was then
committed/pushed as `74d6bfa` BEFORE new adapter/group evaluation. The
[protocol](../research/carbon-H-factorial.md) and
[admission evidence](../research/carbon-H-factorial-admission-evidence.json)
retain all 40 contexts: ALL 24 original 1SHG reconstruction-stress starts and
ALL 16 fixed eight-input/two-seed capture entries. All 34 native starts pass
transport and CH2/CH3 construction admission; six inherited source/seed refusals
remain unchanged. Four arms imply 160 entries, 136 possible replays and 24
inherited refused arm entries. No factorial outcomes evaluated at this checkpoint.

Legacy labels exactly match retained full-topology identities; actual native
HarmonicBondForce graph equals saved topology bonds. Original nm starts/endpoints
match retained Å snapshots exactly. Modern V2 topology/ownership records match
native graph/masses. Full case-sensitive identities, pre-H observations and
elements are exact; positive masses equal singly bonded newly added H only.
No guessed aliases, measured/kept H rewrite, constraints/periodic/virtual sites.
Resource/transport failure blocks rather than becoming a geometry exclusion.

Fixed-zero methyl and heavy-compatible CH2 constructors preserve original radii
and observations. Combined ownership is disjoint and reads both constructions
from original coordinates, no sequential energy/phase search. All initial heavy
shape findings and exclusions remain. Across native starts: 1,472 methyl and
3,326 methylene group instances, not independent protein samples. Poor-heavy
1SHG stress conditions remain explicit, not healthy production admissions.

Enrollment SHA `ec69ad151110223d96725ea85ed8a108d139806f42e7eded9362b3cc9c428ba7`;
admission SHA `6474696f4b4b63bad7ce78b55d8bcaf8842aca6a57bf96a3fc364b332a51394b`.
Private admission V1 is retained as an earlier transport-only report; V2 adds
continuous initial-heavy findings and explicit replay labels without arm execution.
Twenty focused tests and Ruff pass; independent admission regeneration is exact.

**Resume:** freeze complete factorial execution implementation/plan and run all
four arms (unchanged/all_methyl_zero/all_methylene/combined_zero), original
Reference/1.0/50 physics, no new H/seed/force/phase choice or retries. Budgets
900 seconds/cell, 14,400 seconds total. Require all 34 unchanged nm endpoints
exact, all native callback/mobile-gradient readbacks and adverse global H geometry.
Known graph-zero protein-view and 1TIM CH2 intersections are replay controls,
not independent new evidence; do not conflate all-methyl 1SHG with the old
THR24-only intervention. Keep six inherited refusals in all four arms and report
stress/natural strata separately. Production/public repair, original unforced
THR2 cause and empirical-table production/license approval remain unresolved.
Continue autonomously until genuine user access, decision or authorization needed.

## Preceding checkpoint: local CH2 rescue transfers whole-phase flattening to ALA151

The H-stage checkpoint was committed/pushed as `fb4fd4e`. The
[complete whole continuation](../research/ASP49-methylene-whole.md) and
[evidence](../research/ASP49-methylene-whole-evidence.json) execute all 18 cells /
54 states, no retries/blockers/timeouts/budget nonexecutions. All six unchanged
starts AND both corresponding prior methyl-source-arm whole endpoints exactly
reproduce in nm, every atom. Native H-label repair/object-index bridge and every
phase System/anchor remain exact. Original Reference/no-solvent/tolerance 10/
max 1,000 iterations per phase, tenfold native phase-2 weight reduction. No
emergency reflection, observation restore, new seed, public candidate/PDB or
production/default/threshold change. These are pre-emergency native phases only.

Local two-H placement prevents ASP49 flattening in all six contexts: final CA
volumes ~2.609–2.613 Å³. For seed 1/unchanged methyl source, prior ASP49 volume
2.575429→1.065215→0.500677 becomes 2.575429→2.588732→2.612701 with local CH2
(all-carbon 2.610218 final). But seed 2/graph-zero source transfers the failure
to ALA151: original heavy volume 2.315224 in all arms, local phase 1 0.460936 and
phase 2 **0.445384 Å³**, versus matched unchanged **2.466134** and all-carbon
**2.507008**. Final local ALA151 HA remains inward (cos −0.982127), N/C/CB angles
84.318/91.124/85.969°. No final all-acute HA, so that descriptor misses the damage.

Global final half-Å³ CA flags: unchanged/local/all-carbon **1/1/0**. Equal aggregate
counts hide the changed site; unchanged flag is known ASP49 selected-source seed 1
0.498283, with other native source near 0.500677 just above threshold. New ALA151
0.445 is NOT merely borderline classification. No negative CA or THR-beta shape
finding in any arm. All raw/rounded centre volumes and both sites are published.

All-carbon minimum final CA across six contexts: 2.162403 Å³, max H-angle deviation
11.915° versus 28.629/25.182 unchanged/local. This remains ONE protein with
correlated arms and 620 changed H, not general/complete geometry qualification.
Original 1,870 heavy observations move (new ALA151 context RMS
0.285218/0.287578/0.284049 Å); no coordinate restoration or native accuracy score.

Whole plan SHA256 `1d9871f16b8b53aade2f9415141675833707f3de857a9865bf81f88ae6e8a238`;
private `ASP49-methylene-whole-{plan,run,audit}-01*`. Five focused tests/Ruff pass;
independent publications byte-identical. Full fast suite through initial
downstream tests: 3,004 passed, three optional skips/33 slow deselected; later
end-to-end saved-only audit test separately passes. Production stays unchanged.

**Resume:** freeze a source-diverse carbon-multi-H FACTORIAL cohort before any
new evaluation: all existing native starts from the eight-input convenience
enrollment plus ALL 24 original 1SHG reconstruction-stress contexts, not a
success-selected replacement/seed search. Compare unchanged / all methyl phase 0 /
all carbon-CH2 / their fixed-zero combination using the two existing constructors.
Retain all six inherited source/seed refusals, known positive replays and bad-heavy
stress contexts explicitly; unchanged endpoints must match. No geometry cleanup,
energy-selected phase search, new force/seed/iteration profile or production
capability. Original unforced THR2 Cβ causality remains unproven; do not expand
seed/platform search until a desired inversion appears. Empirical-table
production provenance/license/approval gates remain closed. Continue autonomously
until genuine user access, decision or authorization is required.

## Preceding checkpoint: two-H CH2 control rescues ASP49 but creates distant ALA151 failure

Saved-only attribution was committed/pushed as `8225b95`. The subsequent
[prospective methylene control](../research/ASP49-methylene-control.md) and
[full evidence](../research/ASP49-methylene-control-evidence.json) execute all
18 cells: all six 1TIM source contexts × unchanged/local LEU48 CH2/all-carbon CH2.
No retries, blockers, timeouts or budget nonexecutions. All six unchanged nm
endpoints exactly reproduce their prior methyl-source-arm endpoints. All nonowned
initial coordinates and fixed trajectory atoms remain exact; mobile-H forces
equal negative native gradients at all 900 callbacks (max error 0.0), 936 nodes.

Only newly mobile carbon H are owned. Local changes just LEU48 HB2/HB3 initially;
all-carbon changes 620 H at 310 CH2 parents. Lexicographic heavy neighbors define
bisector/normal, H directions a*b +/-sqrt(1-a²)*n with a=-1/(3*cos(theta/2)); first
lexicographic H takes +normal. Original H radii and all heavy geometry stay fixed;
each H-heavy cosine is -1/3, but H-H is IMPLIED by the actual heavy angle, not
independently idealized. Undefined/collinear/infeasible frames refuse, no fallback,
phase/label/energy/native search. Same original Reference/tolerance 1.0/max 50.

Any focal ASP49 inward callbacks: unchanged 3/6, local 0/6, all-carbon 0/6. Final
focal inward: 2/6, 0/6, 0/6. Final GLOBAL all-acute HA count: 2/1/0. New local-arm
failure: seed 2/graph-tetrahedral-zero source -> ALA151 HA, N/C/CB-to-HA
63.558/87.155/71.232°, radius 1.073012 Å. Original LEU48 CA/CB are 27.429/26.700 Å
from ALA151 CA; newly placed H are 26.450/25.783 Å from ALA151 HA, beyond native
10 Å cutoff without pair exceptions. No direct initial local-pair force there;
indirect environment vs global optimizer history is NOT isolated by this result.
Preserve distant regression, not just the rescued focal HA.

Large final H strain persists: maxima across arms ~72.281/70.937/71.256° and
H-bond deviations up to 0.228/0.195/0.167 Å. Seed 1/unchanged methyl source H-angle
maximum worsens 52.365→70.160/71.256 under local/all-carbon. Neither narrow positive
HA descriptor nor all-carbon's zero final all-acute count qualifies full H geometry.
No heavy atoms move in H stage, so downstream CA behavior remains separate.

Plan SHA256 `36131c49174fb83ef637b908ba4d8bd2399a11dbf3ac28ffbca85d0673b4f74f`;
private `ASP49-methylene-control-{plan,run,audit}-01*`. Twenty-two focused tests
and Ruff pass. Full fast suite through initial control/audit: 2,999 passed, three
optional skips/33 slow deselected; later explicit implied-H-H test separately
passes. Independent publications byte-identical. Production/defaults untouched.

**Resume:** separately freeze whole native-phase continuation for ALL 18 outputs,
requiring all six unchanged starts AND both endpoints to match the corresponding
previous 1TIM/methyl-source-arm whole trajectories exactly. Reuse native H-label
repair/object-index bridge and original Systems/anchors/Reference/no-solvent/
10 kJ/mol/nm/max 1,000 per phase/tenfold phase-2 restraint reduction. Report ASP49,
ALA151, all global volumes and H strain plus original-heavy movement. No emergency
safety-net execution in this pre-recovery research comparison; no candidate/PDB.
Retain every adverse state, no favorable source/seed/phase choice. Original unforced
THR2 cause, broad qualification and empirical-table production/license approval
remain unresolved; continue until genuine user input is required.

## Preceding checkpoint: delayed ASP49 inward paths implicate nonmethyl H contacts

The diverse downstream checkpoint was committed/pushed as `bfd0b84`. The
[saved-only ASP49 attribution](../research/graph-methyl-ASP49-attribution.md) and
[independent evidence](../research/graph-methyl-ASP49-attribution-evidence.json)
retain ALL three 1TIM arms × seeds1/2: six contexts, 312 nodes, 1,172,184 exact pair
observations. No new minimization/integration/H/seed, force/exclusion edits or
native error score. Analytic 3D sums close to original native NonbondedForce with
max component error 1.048e−9 kJ/mol/nm; summed force groups also match original
callback gradients (max 5.821e−11). All raw pairs stay private/checksummed.

Initial seed-1 HA nonbonded outward forces are +16.424/+17.077/+18.744 in
unchanged/zero/selected, with same bond +2538.318 and angle +119.911 contributions.
There is NO large initial inward impulse like earlier 1SHG. Seed-1 first negative
outward cosines occur only at callbacks11/12 (unchanged/selected), zero never.
Unchanged callback10: LEU48 HB2 at1.240141 Å gives −145750.713 outward projection,
almost all LJ (−145715.849), versus all methyl H together only−29.343. The native
graph makes HB2 a methylene H (CB bonded CA/CG/HB2/HB3), untouched initially by
the methyl-only arms. Selected callback11: same HB2 at1.146865 Å gives−416822.384.
After crossing, dominant pairs change; retain every node, not a single-pair story.

Selected seed2 has an additional transient inward path at callbacks26–30, final
outward cosine+0.984010. ILE86 HD13 methyl-H dominates before that crossing, then
ASP49 CG; positive total outward force at crossing cautions against interpreting
a quasi-Newton/line-search path as one instantaneous scalar-force direction.
All six matched trajectories retained. This is exact attribution, NOT controlled
pair ablation or proof of original THR2 cause/production repair.

Attribution plan SHA256
`b7273ab667e61c394fb9717516d41ea7d3ca6b5f6adfc13c95e3f594e858c1e2`;
private `graph-methyl-ASP49-attribution-{plan,run,audit}-01*`. Six focused tests and
Ruff pass; independent publications byte-identical. Latest full fast suite through
this attribution checkpoint: 2,978 passed, three optional skips/33 slow deselected.
Production/defaults/thresholds untouched.

**Resume:** separately freeze a prospective initial-condition control on the TWO
generated LEU48 CB H, with all three existing source arms and both seeds. Compare
unchanged/local canonical CH2/all-carbon canonical CH2 under original Reference/
1.0/50 physics. Deterministic observed-heavy-plane construction and name-ordered
slot assignment, preserve each original C–H radius and every nonowned coordinate;
refuse undefined/infeasible frames, no energy/native-based slot selection. Keep
all controls and adverse outcomes; no iteration/seed tuning or production path.
Original unforced 1CRN THR2 Cβ causality and empirical-table production provenance/
license/approval remain unresolved. Continue autonomously until genuine user
access, decision or authorization is needed.

## Preceding checkpoint: source-diverse whole phases retain ASP49 flattening

The graph H-stage checkpoint was committed/pushed as `8ceab31`. The
[downstream note](../research/graph-methyl-whole.md),
[complete continuous centre-level evidence](../research/graph-methyl-whole-centres-evidence.json)
and preserved [compact evidence](../research/graph-methyl-whole-evidence.json)
retain all 48 continuation cells: 30 completed, 18 source refusals, zero
retries/blockers/timeouts/nonexecutions. Ten independent saved-State native
baselines precede continuations; all ten unchanged starts AND both endpoints
reproduce exactly in nm, every atom. All phase Systems/anchors agree with baseline.

Native H-label repair once, followed by actual object-index mapping to full
labels; all ten original starts preflight exactly. Original Reference/no-solvent/
tolerance 10 kJ/mol/nm/max 1,000 iterations per phase; only original positional
weights reduced tenfold in phase 2. No emergency reflection, observation restore,
new H/seed, force edits or public candidate/PDB. These are two native phases BEFORE
production's unchanged emergency safety net, not complete production outputs.

Important continuous result: 1TIM/A/ASP49/seed1 starts at heavy volume 2.575429 Å³
in all arms. Unchanged/zero/selected phase-1 volumes are 1.065215/2.577428/1.060820;
phase-2 volumes are **0.500677/2.609578/0.498283 Å³**. Only selected crosses the old
0.5 Å³ shape descriptor (rounded 0.497121), but unchanged is similarly flattened
and borderline above it (rounded 0.500967). Do NOT exaggerate their binary split
or retune the threshold. All raw/rounded CA/THR-beta volumes at all 90 states are
published so both continuous outcomes remain visible. No final negative CA or
THR-beta in any arm. Final all-acute HA count is zero, but unchanged/selected
ASP49 HA angles remain ~81/94/86°; that narrow H descriptor misses heavy flattening.

Original 1TIM heavy observations move: 1,870 atoms, RMS
0.288697/0.285958/0.289334 Å, max 1.892546/1.880936/1.932975 Å. This is NOT
observation-immutable repair. Global-only half-Å³ descriptors are not the earlier
masked-target gate or complete scientific acceptance. Earlier fixed-phase 1SHG
TYR15/negative-CA regressions still preclude qualifying zero phase.

Plan SHA256 `cb35b44b30346854d3b4a216271120d0fddec81959085f74fa5abb0ba1a680eb`;
private `graph-methyl-whole-{plan,baseline,run,audit}-01*` under research temp root.
Five focused native-baseline/transport/audit tests and Ruff pass; independent
expanded evidence publications agree byte-for-byte. Full fast suite through native
baseline integration: 2,969 passed, three optional skips/33 slow deselected, with
later audit tests separately exercised. Production/defaults/thresholds untouched.

**Resume:** separately freeze a saved-only exact force/pair attribution on 1TIM
ASP49 using ALL three arms and both existing seeds, no additional optimization,
phase/seed tuning or native score. Determine actual contributors, not presumed
methyl causality. Preserve matched negative controls and every source refusal.
This is a second protein-context H-placement/flattening observation, NOT resolution
of original unforced 1CRN THR2 Cβ causality. Production/public repair and original
table provenance/license/approval gates remain closed; continue autonomously until
genuine user access, decision or authorization is required.

## Preceding checkpoint: graph-methyl selected-phase benefit does not generalize

Admission/preparation was committed/pushed as `a555f17`. The completed
[source-diverse H control](../research/methyl-diversity-control.md) and
[full evidence](../research/graph-methyl-control-evidence.json) retain all 48 cells:
30 completed, 18 source refusals, zero retries/blockers/timeouts/nonexecutions.
All 10 unchanged native nm endpoints exactly reproduce. Nonowned initial and all
heavy trajectory coordinates remain exact. Every mobile-H force matches negative
native gradients at all 1,500 callbacks (error 0.0); all 1,560 nodes and 15,168 original
initial-energy observations independently audited (maximum grid error 1.770e−8).

[Native capture evidence](../research/methyl-diversity-capture-evidence.json):
ten contexts, 632 methyl group instances (five protein views × seeds 1/2), complete
object-index/bond lineage through label restoration. Per-seed groups:
1CRN 27 / 1CTF 50 / 1TIM 154 / 1EMV 42 / 1FR2 43. First V1 collector failed ALL 10 attempts
BEFORE H minimization because its ForceField wrapper hid ArgTracker defaults;
retain `methyl-diversity-capture-01/` and immutable V1 script/plan. Explicit V2
signature-preserving amendment succeeds without bypassing argument checks or
changing scientific parameters. Both versions and every failure stay bound.

New adverse result: 1TIM/A/ASP49/seed1 is final all-acute HA in unchanged AND
energy-selected arms; zero-phase prevents it. Any/final all-acute contexts:
1/10 unchanged,0/10 zero,1/10 selected. N/C/CB-to-HA selected angles
59.505/88.192/67.971° versus zero108.094/107.788/109.720°. Initial energy falls
738938.992→449038.657kJ/mol in selected but does not prevent failure. Thus the
earlier 1SHG selected-phase benefit does NOT generalize. Zero also remains unsafe
to qualify: old negative-Cα/TYR15 findings persist, and current 1TIM H strain is
severe (maximum final H-angle deviations across arms71.533/72.281/62.324°).

Graph-H plan SHA256 `a3a23d88c5666e3f790223447a432d1e36b04256420decd265cf603fc32df349`;
private `graph-methyl-control-{plan,run,audit}-01*` under research temp root.
Forty-two focused tests and Ruff pass; independent publications byte-identical.
Initial full fast suite: 2,963 passed, three optional skips/33 slow deselected, with
later focused additions separately tested. Production/public repair stays off.

**Resume:** separately frozen all 30 whole native-phase continuations are underway,
preceded by ten independent saved-State baselines. All 10 H-label repair→object-map
original starts preflight exactly in nm. Use `graph-methyl-whole-plan-01.json`,
`graph-methyl-whole-baseline-01/` and `graph-methyl-whole-run-01/`. Retain source
refusals and all adverse states; unchanged must reproduce both native endpoints
exactly. No solvent, new seeds/forces or production safeguard changes. Two phases
stop BEFORE emergency repair, no PDB/candidate publication. Original THR2 causality,
broader scientific qualification and empirical-table production/license approval
remain unresolved; continue until genuine user access/decision/authorization needed.

## Preceding checkpoint: fixed graph-methyl cohort admission and preparation

The source-diverse [graph-methyl study](../research/methyl-diversity-control.md)
has an immutable [eight-input enrollment](../research/methyl-diversity-enrollment.json),
committed/pushed as `d43d4f5` before admission/implementation/outcomes. Enrollment
SHA256: `716060b015b1e553ef127318e573f1b60ff7cf601d3e9fd117502e16a9a896b5`.
These are convenience fixtures previously used in other research, not blinded or
temporal holdout inputs. No replacement/cleanup to gain admission.

[Independent admission/preparation evidence](../research/methyl-diversity-admission-evidence.json):
1CRN/1CTF/1TIM/1EMV/1FR2 admitted and prepared once each. 1UBQ/2VLN/2VLQ refused
with the shared combined altloc/occupancy/numbering reason. All originals and
first-model/first-ATOM-chain projections retained; no rejected view prepared.
Five preparations pass Reference PDBFixer/H seed1 controls, random-state/hook
restoration, unchanged borrowed input and canonical DatRecord decoding. No
execution blocker/timeout/budget nonexecution. Runtime PDBFixer1.12.0/OpenMM8.5.2,
31 bound resources. This is NOT new H-control/scientific acceptance.

Private source plan/projections/prepared artifacts are under
`tmp/pdbfixer-research/methyl-diversity-*`. Nineteen focused enrollment/admission/
preparation/audit tests and Ruff pass; independent publications byte-identical.
Production/defaults/thresholds and all earlier adverse outcomes are unchanged.

**Resume:** implement/test graph methyl geometry with actual native topology,
newly mobile-H ownership and deterministic second-shell heavy frame. Capture
exactly H seeds1/2 and first full-System starts, no whole minimization. Preserve
actual index/object/bond lineage and both H/full labels through variant restoration.
Freeze implementation and captures before evaluating all three enrolled arms.
Selected policy canonicalizes all groups to zero then makes one lexicographic
parent pass over 24 initial-energy phases, no revisit/native/final scoring.
Retain every refusal/blocker; require exact unchanged endpoints and adverse
geometry reporting. Original THR2 causality and production qualification/license
gates remain unresolved; continue autonomously until real user input is needed.

## Preceding checkpoint: downstream gains are limited and fixed-phase regressions remain

The H-stage checkpoint was committed/pushed as `3e9362d`. The
[whole-system continuation](../research/methyl-whole-control.md) and
[complete evidence](../research/methyl-whole-control-evidence.json) execute all
72 cells / 216 states with zero retries/blockers/timeouts/nonexecutions. All 24
unchanged starts AND both original endpoints reproduce exactly, heavy atoms and H.
The original post-H repair bridge is reproduced exactly in nm; all source forces,
anchors, reference volumes and the 10-fold phase-2 weight reduction remain native.

Original Reference / no-solvent / tolerance 10 kJ/mol/nm / 1,000 iterations per
phase. NO new H, seed/profile search, observation restoration, emergency reflection
or public candidate/PDB. These are deliberately TWO NATIVE PHASES BEFORE the
production post-phase-2 reflection safety net; never equate that with a weakened
production zero-D-Cα guarantee. Plan SHA256
`1ff54b50dd4f879dedccb94620a626a4c53645ecf7cfcd4e582af3d2dd0b24a9`;
private `tmp/pdbfixer-research/methyl-whole-{plan,run}-01*`.

Unchanged/fixed/selected raw AND rounded target+global-stereo heavy passes:
10/24, 11/24, 12/24. By source arms: raw 2/8 all three; distal-only 3/8,3/8,4/8;
full-generated 5/8,6/8,6/8. The known GLN16 volume improves
0.482378→2.544202/2.545041 Å³, with final HA angles near tetrahedral and CA–HA
radius ~1.087 Å. Selected also rescues distal-only THR24 seed 2 (target volume
0.009561→2.480606). Still SIX target THR failures + SIX VAL23 naming failures in
selected arm, one protein only. No native RMSD/accuracy or broad chemistry claim.

Explicit adverse fixed-phase findings in ALREADY-REJECTED cells: distal-only
THR24 seed 2 adds observed TYR15 near-planarity (+0.425081 Å³); distal-only seed 1
changes target THR24 to NEGATIVE CA volume −0.014602 Å³. Fixed arm has one final
negative-Cα case; unchanged/selected have zero; no final negative THR-Cβ in any
arm. Both regressions stay reported despite aggregate pass gains. No invalid PDB
is published or candidate accepted. Production reflection safeguard is unchanged.

Actual observed-heavy movement is reported, not restored: focal GLN16-rescued case
has 469 original heavy observations, RMS 0.219522/0.219133/0.218328 Å and max
1.077288/1.098176/1.084132 Å (unchanged/fixed/selected). Heavy gates cover target
geometry and global Cα/THR-Cβ shape/stereo, not exhaustive complete-complex/H checks;
all H-bond/H-angle deviations are retained descriptively.

Ten new focused tests and Ruff pass. Full fast suite through continuation:
**2,906 passed**, three optional Torch/Biotite skips, 33 slow deselected. Independent
publications regenerate byte-identically. User requests continued autonomous
research until a genuine access/decision/authorization requirement.

**Resume:** separately enroll a source-diverse, graph-defined generated-methyl
control with fixed input, ownership, canonical-frame and initial-energy rules.
Retain refusals/adverse contexts; do not select proteins/angles by final success,
enable a public repair or weaken the production chirality safety net. Original
unforced 1CRN THR2 Cβ causality, broader qualification and empirical-table
production/license approval remain separate and unresolved.

## Preceding checkpoint: tetrahedral methyl H-stage control avoids the inward endpoint

The negative rigid-control checkpoint was committed/pushed as `21d68b1`.
[Tetrahedral H-stage evidence](../research/tetra-methyl-control-evidence.json) and
[note](../research/tetra-methyl-control.md) execute all 72 declared cells:
24 original starts × unchanged / fixed-phase canonical / initial-energy-selected
canonical methyl. All complete. All 24 unchanged endpoints exactly reproduce.
Every nonmethyl initial coordinate and heavy coordinate at all 3,600 callbacks
is exact; mobile-H forces equal negative native gradients with error 0.0.

Original C–H radii stay fixed; only THR24 HG21/HG22/HG23 get canonical 109.471°
angular geometry in a fixed CB/CG2/CA frame. Native original Reference / tolerance
1.0 / max 50 iterations and all forces/masses remain unchanged. Canonical phase
zero is independent of a 24-phase initial-energy-selected arm (0–345°, 15° grid,
same frozen tie rule). No new seed/H/force/exclusion or whole-system minimization
in this H-stage experiment. Plan SHA256
`97d88724b5f984f87512d500357ff2b8c4096c62590c5d4b3198340e0cfd2102`;
private `tmp/pdbfixer-research/tetra-methyl-{plan,run}-01*`.

Any-inward callbacks: 2/24 unchanged, 1/24 fixed phase, 0/24 energy-selected.
Any/final all-acute: 1/24 unchanged, 0/24 in either canonical arm. Retain the
fixed-phase transient crossing. Focal selected phase is 30°; initial HA nonbonded
outward force reduces −56105.474→−13734.044 kJ/mol/nm. This narrow positive result
does NOT qualify complete H geometry: focal final CA–HA is 0.971753 Å (about
0.118 Å compressed vs FF reference), with N/CB/C-to-HA angles
103.584/89.102/129.787°. Initial callback bond compression also remains severe.

Sixteen new focused tests and Ruff pass; full fast suite through H stage:
**2,896 passed**, three optional Torch/Biotite skips and 33 slow deselected.
Independent publications regenerate byte-identically. Production is unchanged.

Autonomous continuation has separately preflighted the original post-H repair
bridge: all 24 original whole-system starts reproduce EXACTLY in nm, not merely
rounded Å. A whole-system control for all 72 H outputs has its own frozen plan
and raw/audit artifacts under `tmp/pdbfixer-research/methyl-whole-*-01*`; do not
infer its acceptance from this H-stage count. Transport and whole-system code
have separate focused tests. Old original THR2 cause, remaining target failures,
broader chemistry/platform qualification and licensed-table approval stay open.

## Preceding checkpoint: energy-selected rigid methyl rotation is negative

Following `c41dfcc`, the [rigid methyl-start control](../research/methyl-start-control.md)
and [complete evidence](../research/methyl-start-control-evidence.json) execute
ALL 48 cells (24 original starts × unchanged/initial-energy-selected methyl
rotation). All complete; all 24 unchanged endpoints reproduce exactly. All
nonmethyl initial coordinates and existing-heavy coordinates at every callback
are exact. Mobile-H physical forces equal negative native gradients at all 2,400
callbacks, maximum component error 0.0 kJ/mol/nm; 2,496 nodes retained.

Only generated THR24 HG21/HG22/HG23 rotate rigidly about CG2→CB, preserving their
original C–H radii and internal geometry. Fixed grid: 0–345° in 15° steps, minimum
original initial System energy, 1e−8 kJ/mol tie tolerance, shortest signed rotation.
No native/final-geometry scoring, H regeneration, force/exclusion edits or seed
search. Original Reference / tolerance 1.0 / max 50 iterations remains fixed;
NO whole-system minimization or observation rewrite. Plan SHA256
`d52e8c56ca650e04a212737c9f9110169f15063c52de3d0a86d7a62643560461`;
private `tmp/pdbfixer-research/methyl-start-{plan,run}-01*`.

Both focal full-generated THR24 seeds choose −30°. Seed 1 initial energy decreases
145674.319→145561.857 kJ/mol and nonbonded outward HA force weakens
−56105.474→−46482.938, but still ends all-acute/inward. Final all-acute count stays
1/24; any-inward callback stays 2/24; any-all-acute callback INCREASES 1/24→2/24
(selected seed 2 transiently all-acute, though finally outward). Keep this negative
result and the regression, not just energy improvement.

Rigid rotation preserves the native initial methyl-H cluster: focal H–CG2–H
angles are 6.13/18.55/18.88° for seed 1 and 4.22/2.97/7.10° for seed 2. This is
not a trial of canonical tetrahedral initialization. A separate next control is
being prepared: unchanged / fixed-phase tetrahedral / energy-selected tetrahedral,
same original C–H radii, original forces and all matched starts. Freeze/test its
frame, 24 phases, energy/tie rule and budgets BEFORE any new execution.

Seventeen focused tests and Ruff pass; independent publications byte-identical.
Full fast suite through this checkpoint: **2,880 passed**, three optional
Torch/Biotite skips, 33 slow deselected. Production/public repair and empirical-
table license/approval gates remain unchanged. Original THR2 causality remains
separate and unresolved. User explicitly requests autonomous continuation until
a genuine decision/access/authorization requirement, not a stop at each checkpoint.

## Preceding checkpoint: exact pairs identify experimental THR24 methyl repulsion

The native-path checkpoint was committed/pushed as `29ca966`. The subsequent
[pair attribution](../research/ha-nonbonded-pairs.md) and
[complete evidence](../research/ha-nonbonded-pairs-evidence.json) check ALL 24 cells /
1,248 nodes / 1,190,592 HA-other-particle observations. Analytic full-vector sums
match the already measured UNMODIFIED native NonbondedForce with maximum error
8.731e−10 kJ/mol/nm. Independent streaming audit checks every raw pair inventory,
identity, vector/projection/class sum and fixed three-most-negative inspection
rule. No new minimization, coordinate perturbation, pair deletion or native score.

Plan: `tmp/pdbfixer-research/ha-nonbonded-pairs-plan-01.json`, SHA256
`7bb42ad94178f7d4817bfeaa5d1c6d1c025871aceb2d2618d2d0b492a9294735`.
Full private pair archives: `ha-nonbonded-pairs-run-01/` and `-02/` beside it;
no full raw force/coordinate archives or FF parameter tables go into Git.

Known full-generated THR24 seed 1: initial GLN16 HA–THR24 HG23 distance 1.356899 Å,
outward projection −43245.578 (LJ −43207.680, Coulomb −37.897); HA–CG2 contributes
−6494.745 and HA–HG21 −3687.309 kJ/mol/nm. These are ordinary inter-residue pairs.
The target methyl belongs to the EXPERIMENTAL rebuilt heavy candidate, not the
unmodified deposited source. Local heavy geometry alone did not ensure safe H
environment. At callback 0, the dominant inward contact is now its own GLN16 O
at 1.408459 Å, an ORIGINAL exception pair contributing −51601.067; do not describe
one initial contact as the only force source throughout the optimizer path.

Seed 2 initially has even stronger methyl-H contacts (HG23/HG21/HG22) but returns
outward, so pair attribution is NOT a universal outcome predictor or controlled
pair ablation. Original THR2 causality, remaining target THR/VAL failures and
broader scientific/platform qualification stay separate. Fourteen new focused
tests and Ruff pass; combined angle/H-stage/pair controls pass 73 tests. Both
complete private pair runs and independent publications are byte-identical.
Full fast suite: **2,863 passed**, three optional Torch/Biotite skips and 33 slow
deselected. Agent-map and generated CLI/GUI checks pass. Production is unchanged.

The next step declared at this checkpoint was a new environment-aware initial-H/candidate-geometry mechanism control
would need separately declared settings and all matched cases retained. Do not
silently alter exclusions, ablate pairs, tune seeds/iterations to the known result,
or promote a neighbor-HA rewrite. Preserve the empirical-table production/license
gate. Pair forces identify recorded contributors, not unique optimizer causality.

## Preceding checkpoint: large early native H-minimizer displacement

Following checkpoint `3c8c2ff`, the
[native callback path](../research/ha-minimizer-path.md) and
[complete evidence](../research/ha-minimizer-path-evidence.json) retain ALL 24
cells / 1,200 callbacks / 1,248 initial+callback+final force-readback nodes.
Every original start/endpoint reproduces exactly; all existing heavy coordinates
are exact at every callback, and reported physical HA forces match negative native
objective gradients with maximum difference 0.0 kJ/mol/nm. All cells execute;
zero retries/blockers/timeouts/nonexecutions. Original Reference / 1.0 tolerance /
50-iteration settings stay fixed; callbacks only copy data and return False.
No constraints, PBC, virtual sites, integration or whole-system minimization.

Plan: `tmp/pdbfixer-research/ha-minimizer-path-plan-01.json`, SHA256
`4d159a5be0823c18577e0cb821c9a534361d6a41c28cb572617369137c0cfb2e`;
private path coordinates/gradients in `ha-minimizer-path-run-01/` beside it.
Observe the native Context directly; do not reconstruct a supposedly exact nm
start by round-tripping saved Å decimals. Old source/evidence stays immutable.

For known THR24/full-generated/seed-1 HA, outward cosine first becomes negative
at callback ordinal/native index 0 (+0.962637→−0.912389), while CA–HA radius grows
1.000000→2.506358 Å, then to 4.430951 Å at callback 1. First all-acute callback
is 17; final cosine −0.994893. Initial inward-axis force is −53175.879 kJ/mol/nm,
dominated by nonbonded −56105.474, opposed by bond +2464.937 and angle +464.623.
Initial sphere-tangent force is OUTWARD (+5524.842), including nonbonded +3808.931;
do not call the displacement ordinary bond-preserving inward rotation. Callbacks
are optimizer observations, not physical dynamics, every line-search evaluation,
integrated work or proof of convergence. All cells use their full 50 callbacks.

Only two full-generated THR24 starts show negative-outward callbacks; seed 2
returns outward and never has an all-acute callback. Heavy Cα remains fixed here;
its later flattening belongs to the separately reproduced paired intervention.
No particular nonbonded pair has yet been identified. Eleven new focused tests
and Ruff pass; combined angle/H-stage controls pass 59 tests. Repeat publication
is byte-identical. The new full fast suite passes **2,849 tests**, with three
optional Torch/Biotite skips and 33 slow deselected. Agent-map and generated
CLI/GUI checks pass. Production is unchanged.

The next step declared at this checkpoint was to separately freeze exact HA-touched nonbonded pair decomposition and
check full pair sums against the unmodified force before naming a clash partner.
Do not choose a neighbor solely by proximity, suppress failures, increase the
iteration budget or integrate an outward-H repair. Original THR2 cause, remaining
target THR/VAL failures, platform/diversity qualification and licensed-table
production approval remain open and distinct.

## Preceding checkpoint: inward HA originates during native H-only minimization

The paired-HA checkpoint was committed/pushed as `ae0a08c`. The next
[original H-placement trace](../research/ha-placement-trace.md) and
[complete evidence](../research/ha-placement-trace-evidence.json) replay all 24
original cells against the preserved historical source. All complete; all H
endpoints and first whole-system starts/System parameters reproduce exactly.
Plan SHA256 `92e057f80deb827070721c2832f4de0c96aca52ee214dedd0ef6227342be8919`;
private artifacts: `tmp/pdbfixer-research/ha-placement-{plan,run}-01*`.

Observe only the existing OpenMM 8.5.2 H-only minimization (Reference, tolerance
1.0, maximum 50 iterations), then stop before whole-system phase 1. For the known
THR24/full-generated/seed-1 case, initial N/CB/C-to-HA angles are
117.734/91.407/114.548° (all obtuse), ending 64.207/72.748/78.228° (all acute).
All 472 existing heavy atoms are EXACTLY immobile. The inward HA therefore
appears during native H-only minimization; the heavy Cα does not flatten at this
stage. Its later full-system flattening is the separately reproduced paired
trajectory, not an event during this heavy-fixed stage.

Initial descriptors: 12 all-obtuse / 12 mixed / zero all-acute; final: 18/5/1.
All 24 H endpoints and whole-system starts match saved originals at 0.0 Å,
including hydrogens. No alternative H placement/settings/seed or force change,
no full-system minimization, source rewrite, public repair or production change.
Eight new focused tests and Ruff pass; independent publication regenerates
byte-identically. Previous full fast suite is 2,830 passed; this increment has
focused acceptance only, not an inferred new full-suite pass.

The next step declared at this checkpoint was to separately freeze an original-H-System intermediate-step or force
readback before identifying a microscopic force/line-search mechanism. Do not
increase iterations or promote the one-H construction as a fix. Original 1CRN
THR2 causality, remaining target THR Cα failures, diversity/platform qualification
and empirical-table production/license gate remain unresolved and separate.

## Preceding checkpoint: paired one-HA intervention prevents the GLN16 failure

The angle-attribution checkpoint was committed/pushed as `d972aff`. The subsequent
[paired HA-start study](../research/ha-start-control.md) and
[complete evidence](../research/ha-start-control-evidence.json) execute all 48
declared cells (24 original starts × unchanged/outward HA), with zero execution
blockers/timeouts/nonexecutions. Plan SHA256:
`ec5500f041906f721a5ad10aecebe6c230304f75905687060cfb66ad8346d1ae`.
Private plan/run: `tmp/pdbfixer-research/ha-start-plan-01.json` and
`tmp/pdbfixer-research/ha-start-run-01/`; no private force exports in Git.

All 24 unchanged controls reproduce BOTH saved native-phase endpoints exactly,
all heavy atoms and hydrogens (maximum 0.0 Å difference). Independent audit
confirms only generated GLN16 HA changes at the paired starts, CA–HA distance is
preserved, and every saved System matches its original phase's parameters,
anchors and reference volumes. No H regeneration, force-strength search,
reflection, observation restoration, solvent or production code change.

For full-generated THR24 seed 1, final GLN16 Cα volume changes from 0.482378 to
2.546770 Å³, with HA angles restored to 109.488/107.942/108.227° (N/CB/C).
This supports a causal H-start effect under the reproduced physics for this
specific GLN16 failure, NOT original 1CRN THR2 causality or a production repair.
The original 469 heavy observations move during minimization: final RMS
0.219522/0.218545 Å and maximum 1.077288/1.097194 Å (unchanged/outward).
Never claim observation-immutable final coordinates.

Whole-cohort raw/rounded passes change only 10/24→11/24, specifically full-heavy
5/8→6/8. Raw stays 2/8 and distal-only stays 3/8. Seven remaining THR Cα failures
and six VAL23 naming failures remain retained. No native RMSD or earlier accuracy
result was recomputed. All controls are one protein, not independent diversity.

Driver/auditor: 18 focused tests passed; combined angle/HA set 40 passed; targeted
Ruff passed; independent audit/publication regenerate byte-identically. The new
full fast suite passed **2,830 tests**, three optional Torch/Biotite skips and
33 slow deselected. Prior angle fast suite: 2,812 passed, three optional skips,
33 slow deselected. Driver enrollment initially failed on an absent historical
trace `settings` key before writing a plan or running any cell; corrected before
freezing, with historical CLI/OpenMM defaults independently read back.

The next step declared at this checkpoint was to inspect the original H-placement stage's local geometry provenance
with separately frozen settings/observations, not a public HA-fix integration.
Keep the upstream H-placement origin, remaining target THR flattening and original
unforced THR2 event separate. Original empirical-table production/license gate
and all historical failure/source bindings remain unchanged.

## Preceding checkpoint: HA-containing angles dominate the GLN16 flattening component

On explicit user request, the complete earlier research checkpoint was committed
and pushed as `421b2a3` on `research/pdbfixer-reconstruction`. Private raw runs,
environments, upstream rotamer tables and credentials were excluded. This is a
research-source/evidence branch, NOT production qualification, package release,
public completion capability or legal-clearance claim. Main remains untouched.

Read [angle-term attribution](../research/stereo-angle-terms.md),
[complete angle evidence](../research/stereo-angle-terms-evidence.json) and
[initial-context identity checks](../research/stereo-angle-context-evidence.json).
All 24 cells / 96 states / 240 focal centres and 6,096 touched-angle observations are
retained; zero blockers. Analytic term sums match the old OpenMM angle-group
projections to 3.606e−12 kJ/mol/nm. Independent finite differences check all 6,096
terms at step 1e−7 nm; maximum absolute difference 7.460e−6 kJ/mol/nm. No Context,
integration, minimization, H placement or native scoring on the cohort here.

For full-generated THR24 seed 1 observed GLN16, final H-containing projection is
−3793.179 and heavy-only +1365.516 kJ/mol/nm. The three dominant negative terms
are N–CA–HA, HA–CA–CB and HA–CA–C: initial 64.207/72.748/78.228°, final 85.944/
85.697/89.067°, references 109.5°. Thus the heavy-only angle group is restoring
at this endpoint; the large negative component comes from H-containing angles.
Do not conflate a force projection with constrained optimizer motion/stationarity.

All 24 initial local GLN16 heavy-coordinate inventories have one hash, and all
touched-angle parameter inventories have one hash. HA angle triads are all obtuse
in 18 cases, mixed/boundary in five, all acute in one (the known full-generated THR24 seed-1
case). These are geometric descriptors, NOT new gates/CIP calls/population rates.
The initial raw/full harmonic-angle difference is entirely in H-containing terms
at identical local heavy coordinates and angle parameters. No causal replay or
original unforced 1CRN THR2 explanation is claimed. Old failures/accuracy/source-H
refusals and original frozen-source/publication hashes remain unchanged.

Next bounded test: matched saved-System replays changing ONLY generated GLN16
HA's initial position to an outward geometric construction, preserving its bond
length, every other starting coordinate and every force/anchor/reference. Keep
all 24 starts plus unchanged controls; no native-error selection, new FF strengths,
H regeneration, branch swaps, source-heavy rewrite or production correction.
Use the same original two-phase 10 kJ/mol/nm / 1,000-iteration Reference settings and
measure all original-observation motion and raw/rounded whole-context failures.
The replay and its source/force invariants must be declared before execution.

## Research return after diagnostic-only transfer (2026-11-10)

The user explicitly requested transferring the new diagnostics BEFORE continuing
this research branch. Production warning-only diagnostics and their runtime fixes
were completed separately, pushed, and merged into main at `3fc378e`. Original
commits `039e1ba` / `9e71372` were cherry-picked here as `471844a` / `8eeeada`;
there was no full main merge, rebase, stash or bulk research commit/push.
The transferred files match `9e71372` byte-for-byte. All 324 pre-existing modified/
untracked research files were checked unchanged immediately after transfer.
The active diagnose/runtime suite passed **164 tests in 20.86 s**.

Old frozen scientific plans bind the pre-diagnostics production source. Do NOT
rewrite their expected hashes to the new source or hide mismatches with skips.
An explicit detached historical worktree at `4a328a6`, with the 324 existing
research files copied and checked, is retained at
`tmp/pdbfixer-research/pre-stereo-baseline-01/`. It shares the original evidence
directory through `tmp/pdbfixer-research`; use it as a READ-ONLY historical-verifier
cwd, not a place to launch experiments or rewrite prior results. Entry manifest
and transfer checks are retained in `tmp/pdbfixer-research/research-return-20261110/`.
The prospective-geometry and start-geometry final verifiers passed there with
unchanged, byte-identical publications. New research binds the active new source
and explicitly identifies the historical source root for prior-plan verification.

The preceding bounded continuation was `scripts/audit_stereo_force_projection.py`:
static Reference force evaluations of ALL 24 saved cohort cells, both sides of
both phases. Focal centres are every target Cα, THR target Cβ only, and the known
observed GLN16 Cα context centre. Enrollment is RETROSPECTIVE to known geometry
outcomes, frozen before force evaluation; not a new blinded/native benchmark.
No integration, minimization, H placement, force/profile/seed search, native
accuracy, coordinate repair, production threshold/default changes or empirical
library tables. The script checks energy/force additivity, exact restoration of
group-only metadata changes, and the canonical Cα wall's analytic derivative.
New plan: `tmp/pdbfixer-research/stereo-force-plan-01.json`; run/readback prefix:
`tmp/pdbfixer-research/stereo-force-readback-01/`. Seven analytic/control tests pass.
Interpret force projections as local Cartesian tendencies, not a constrained
minimizer trajectory or proof of the original unexplained THR2 event.

Readback completed: see [static force analysis](../research/stereo-force-projection.md)
and [complete evidence](../research/stereo-force-projection-evidence.json).
ALL24 cells/96 states/240 focal centre-states evaluated, zero blockers, all24
identical-coordinate phase switches checked. Repeated complete static evaluation
and independent publication readback are byte-identical; not new scientific seeds.
All30 old final verifiers pass from the preserved historical source cwd.

Specific evidence: raw THR24 starts at vref0.025934Å³ and ends at0.008107Å³,
q0.312608, angle0.152°, so the qmin0.2 wall is inactive despite nearplanarity.
Observed-context GLN16 in full-generatedTHR24seed1 ends at q0.190783, volume
0.482378Å³/angle8.964° (unrounded). Its wall correctly restores (+2487.696
kJ/mol/nm projection) against a large flattening harmonic-angle component
(−2427.663); net+67.586 is NOT a constrained/stopping/trajectory certificate.
At the identical phase switch, positional-restraint projection weakens
+1350.708→+135.071, other group projections stay fixed, and net changes
+16.654→−1198.983. This explains local force competition, not unique THR2 cause,
H-only attribution, public repair qualification or a new default force.

The first active focused test run exposed the expected old static-inventory hash
failure (838pass/1fail; retained log stereo-force-focused-tests-01.log). Do NOT
overwrite `pdbfixer-usage-inventory.json`: its original probe and hashes stay
frozen. `amend_pdbfixer_usage_inventory.py` creates the explicit current-source
companion `pdbfixer-usage-inventory-stereo-diagnostics.json`, requiring identical
consumer paths/imports/calls/line locations/counts and only the changed diagnose
file hash. Current inventory tests check that companion; old scientific verifiers
remain anchored to their original context. Seven new math/group controls and seven
independent readback controls pass, plus seven inventory tests after amendment.

Checkpoint checks: **840 active focused research tests passed29.03s**, transferred
diagnose/runtime suite164passed20.86s; ALL30 historical final verifiers pass from
the original baseline, and the new force final verifier passes. TargetedRuff,
17-file librarymypy, CLI/GUI generated checks, agentmaps, whitespace,1,277 local
Markdown targets and all seven documentation-plan files pass. Notfullrepo/slow/
platform acceptance or repair qualification. `validate-stereo-force-final.py`
also proves the original324-file archive intact; only four active research notes
and the explicitly versioned inventory test changed among those original files.
Both force evaluations, the force-audit publication and the source-inventory
amendment regenerate byte-identically. At this checkpoint only the requested
diagnostic transfer had been committed. The user subsequently explicitly
authorized committing and pushing the complete research checkpoint; retain
private raw runs/environments and all licensed upstream tables outside Git.
That publication is research-source/evidence only, not a production integration,
package release, public completion capability or license-clearance claim.

Next bounded mechanism question: disaggregate harmonic-angle terms and H/context
mediation at retained starts/endpoints, with ALL declared cells retained and no
native-error/profile tuning. Do not replace this by a production wall-strength
change or claim that the unforced original THR2 Cβ crossing was reconstructed.

## Earlier checkpoint: new-input qualification fails broad repair claims (2026-11-10)

Read [prospective qualification](../research/prospective-geometry.md),
[geometry ledger](../research/prospective-geometry-evidence.json),
[all-output accuracy](../research/prospective-geometry-native-evidence.json) and
[saved-force disaggregation](../research/start-geometry-mediators-evidence.json).
Standing research-only scope maintained, production/default/dependency/library/
public capability unchanged, no commit/push. No PyRosetta tables used here.

Pre-download1FNA/1SHG/1VII enrollment/unchangedkernel policy,10 masks. Six heavy-
scope refusals retained:1FNA4 (150 depositedH),1VII2 (301H; two countshortfalls).
Original H NOT stripped/normalized to gain admission. These are valid deposited
structures outside prototype scope, not invalidPDBs. Only1SHGfour masks execute;
do NOT call this three independent successfully tested proteins. All4 preparations
complete with recorded PDBFixer/H Reference controls, seed1/defaultlegacy physics.
All24 refinementcells/48 nativephases complete, raw/beta_only/full_generated,
oldseeds1,2, all original observations retained before source rewrite, no newgrid.

Final commongatepasses per8executedcells: raw2, beta_only3, full_generated5.
THRpasses raw0/4, beta1/4, full3/4. This INVALIDATES a general12/12 guarantee from
the old two-protein controls. No acceptance-gate relaxation or source substitution.
VAL23bothseeds/allarms fail val_branch_name_orientation: equivalent methyl names,
NOT an inverted beta stereocentre. RawTHR24/32 andsomebeta arms have flatCA/angle
failures, CB stillR. FullTHR24seed1 targetpasses but OBSERVED GLN16 CA flattens:
2.528406→0.529063→0.482378 Å³ across nativephases, positive NOT D-Cα. Globalgate
correctlyrefuses. No claim whole-trajectory stereo/completeforcelevel cause.

ALL geometry outcomes committed BEFORE nativeanalysis. All24 completed outputs,
including refusedgeometry, reported; all3 armpairs×all/THR/VAL×allcompleted/common-
passing=18 rows fixed in hashed auditor before report. No native-based choices.
Same4 THRcell exactRMSD median raw0.509956, beta0.519561, full0.927699 Å. Full
improves detectedgeometry but NOT coordinateaccuracy. RawlowerRMSD doesn't make
flatCA acceptable. Common-passing raw/fullTHR intersectionEMPTY. All4 nativecontrols
pass; no originalnativegeometry-defect exclusion/filter. Original spontaneous
THR2 beta inversion still not replayed/explained or replaced by artificial stress.

Savedforceanalysis RETROSPECTIVE, nooptimizer/nativecoords: all8 oldTHRbeta_only
contrasts have exact entire serialized force lists, H/positions differ. All8
oldfullgenerated THR contrasts change only newCB anchor and targetCAguardvref;
weights/otherforces identical. EightVAL contrasts no-ops. Narrows mediationcaveat,
does not isolate H vs heavypositions or uniqueoriginalbugcause.

Code: `prospective_geometry_cohort.py`, `audit_prospective_geometry.py`,
`audit_geometry_mediators.py`, tests. Private prefix `tmp/pdbfixer-research/
prospective-geometry-`: enrollment01,inputs01,masks01,prepared01,refinement-plan01,
refinement01,ledger01,native01; forceanalysis `start-geometry-mediators-01.json`.
Newsourcebytes/provenance retained INPROJECT, not registered as trackedfixtures.
Focused825 passed27.85s, targetedRuff/17-file mypy pass. Three independent saved-
only publication readbacks are byte-identical; not fullrepo/platform suite.
Final batch complete:30 current/historical verifiers,55 historical fixture hashes,
1,265 local Markdown targets and fifteen documentation plans pass. Generated
CLI/GUI, agentmaps and whitespace checks pass. Retained log:
`tmp/pdbfixer-research/prospective-geometry-final-checks-01.log`.

### Historical next step, now completed separately: production diagnostics only

The proposal below is retained as historical context; the warning-only work was
subsequently authorized, implemented separately and transferred as described above.
Coordinate repair and empirical-library production approval remain separate gates.

Do NOT integrate the coordinate repair as qualified. A useful next scope is an
independent detection-only extension of `dvbfixer diagnose`: distinguish canonical
THR-Cβ stereochemistry and near-planar Cα geometry from D-Cα. Preserve coordinates,
identities and current prepare/minimize behavior. No automatic correction, branch
swap, new prepbackend, changed defaultforce/profile or public completion route.

This changes production behavior and is NOT authorized by research continuation
alone. Stop and ask the user before implementation. The current report uses all
`chirality` findings as "D-isomer error"; do NOT put Cβ/planarity findings under that
label without separating their semantics/JSON/report handling. Proposed work must
calibrate/test diagnostic criteria and positive/negative controls, not silently
promote the research0.5Å³ gate into an unvalidated productionfatal threshold.
Use production puregeometry/identity/variant helpers, never import researchmodules
or PyRosetta. Warning-only diagnostics and coordinate repair are separate approvals.
New trackedfixture registration would need explicit55-fixture historical-registry
amendment; never silently rewrite oldscientific expectedresults/counts. No push/
release/commit implied. Original-table/license approval remains a DIFFERENT gate
for the empirical PyRosetta route, not a blocker to puregeometry diagnostics.

## Earlier active protocol: prospective new-source qualification (2026-11-10)

User asks to continue autonomously until a real access/scope/production decision.
Saved-system disaggregation completed without optimization/native references:
`scripts/audit_geometry_mediators.py`, private
`tmp/pdbfixer-research/start-geometry-mediators-01.json`.24 arm contrasts: all8
THR beta_only contrasts retain the ENTIRE serialized force list exactly; H/start
coordinates differ. All8 full_generated THR contrasts change only the generated
CB positional anchor and target CA guard reference volume; no restraint weight
change.8 VAL contrasts are identical controls. This narrows the earlier potential-
mediation caveat without claiming original THR2 causality. Publication pending.

BEFORE download froze1FNA (fibronectinIII),1SHG (spectrinSH3),1VII (villin/H and
ensemble retention control), firstMODEL/firstATOMchain/firsttwoTHR+firsttwoVAL,
whole-sidechain masks. Keep originalH/altloc/occupancy/touching metadata/endpoints;
unchanged shared admission, no cleanup/substitution to make a case pass. Up to12
masks/72 refinementcells, raw/beta_only/full_generated andseeds1,2 unchanged.
Reference for PDBFixer/prepH (seed1); refinement explicitly rebuildsH with same
Reference andseed1,2 plus Reference native Context. This is research control,
not production default. One prepare per admittedmask, all failures retained.

Enrollment: `tmp/pdbfixer-research/prospective-geometry-enrollment-01/enrollment.json`.
Official downloads are running to `prospective-geometry-inputs-01/`; preserve
URL/time/originalbytes/hash, do not replace structures. No new-source native
accuracy read or target choices based on outcomes. All completed outputs, including
geometry failures, will be reported only AFTER geometry ledger commit; separately
predeclared common-passing matched arm pairs, no preferred-arm/native selection.
`scripts/prospective_geometry_cohort.py` is a new transport around frozen prior
geometry/worker functions. Scoped PDBFixer resource/Reference preflight passed now;
no environment changes or PyRosetta/table dependency in this continuation.

## Latest completed checkpoint: fixed-chi source geometry yields12/12 final passes in both rebuild arms (2026-11-07)

Read [source geometry factorial](../research/start-geometry.md) and
[independent evidence](../research/start-geometry-evidence.json).
Frozen ALL six original prepared/dat pairs×seeds1,2×raw/beta_only/full_generated,
no native coords/errors, rotamer selection, extra seeds or force/profile grid.
THR beta_only rebuilds OG1/CG2 at OLD CB; full_generated rebuilds CB from observed
N/CA/C then branches. Both preserve original measuredchi1 before PDB rounding.
VAL repair arms explicit exact source no-ops, no beta stereocentre assigned.

Only dat-owned NEW targetheavy atoms may move in source rewrite; all other
records/observedheavy coordinates/metadata/H records preserved, dat byte-identical.
Old H on private rewritten copies deliberately stripped/rebuilt BEFORE physical
minimization. ALL36 cells use --rebuild-h --no-solvent --atom-naming standard
--platform Reference plus explicit Reference INSIDE addHydrogens (CLI platform
alone does not select its internal Context). Fixed Python seed aroundH only,
state/hooks restored; same defaultphysicalsettings, no newbeta force or reflection.

36/36 complete,72 nativephases,0 emergencyCA repair messages. Final commongeometry
raw10/12 (THR6/8), beta_only12/12 (THR8/8), full_generated12/12 (THR8/8). Raw two
failures are alreadywrongTHR6 CB S; both repair arms yieldTHR6 R/finalpass under
bothseeds, retaining sourcechi/observations. Original unforcedTHR2 R→S AGAIN NOT
REPRODUCED in raw Reference/H-control, not uniquely caused/prevented by this test.

Eight VAL no-op contrasts finalheavy EXACT. Max sourcechi roundingchange0.039884°,
unroundedretained1e−8°; roundingmeasured not corrected away. Beta_only initialCA
stillpoor/nearplanar andsomeclashes; full_generated THR1 retains initialOG1–A37O
clash. All initialfailures retained andALLcells executed; finalsuccess is NOT an
inputadmission/generatorclashfree/intermediate-state certificate.
Actualenrollment is SIX masks/twoproteins, fourTHR+twoVAL, withtwoHseeds, not12
independent structures. ILE/PTMs/otherchemistry/partialobservedsidechains unqualified.

Sameprocedure is notidentical H or fixedpotential: changing heavy source geometry
mediates Hpositions; full_generated alsochanges newCBrestraintanchor/CAguardvref.
Distal-only leavesCB geometricanchor fixed. No force-isolation/uniquecauseclaim.
Actualoriginalobservation motion measured against MASKED ORIGINAL source excluding
generatedtargetatoms, separately fromallpreparedheavy motion; no restoration.
Actual-frame raw/rounded measuredchi geometry+globalCA/THRbeta+RDKit audit, H and
everycontextangle notqualified. Source record preservation≠final fixedobservations.

Scripts `research_start_geometry.py`, `audit_start_geometry.py`, tests. Private
`tmp/pdbfixer-research/start-geometry-plan-01.json`, `start-geometry-run-01/`,
originalaudit01/readback01 retain inputchanges/ownership/chi/H/phase/System evidence.
Fresh saved-only publication byte-identical; focused820 passed28.41s, targeted
Ruff/17-file mypy pass, notfullrepo/platform acceptance. No production/default/
dependency/library/public capability change, commit/push; licensegate unchanged.

**Continuation:** positive bounded source-preparation result, NOT publicbackend
qualification or exactoriginalTHR2 causal reconstruction. Do not erase rawfailures,
initialclashes/nearplanarity, change chi/nativeprofile or move observed sourceatoms.
Next meaningful work should separately qualify input-geometry admission/repair and
disaggregate downstreamH/geometry-derivedguards, or enroll new input chemistry;
not another native-error-driven search. Current finalchecks complete:28 final
verifiers,55 fixture hashes,1,256 local Markdown targets and fourteen completed
documentation plans pass. Generated CLI/GUI, agentmaps, whitespace, targetedRuff
and17-file mypy pass. Retained log:
`tmp/pdbfixer-research/start-geometry-final-checks-01.log`.

## Earlier checkpoint: private beta wall preserves sign but not full target geometry (2026-11-07)

Read [private beta-wall stress](../research/beta-wall.md) and
[independent evidence](../research/beta-wall-evidence.json). Production unchanged.
Before execution froze ALL12 previously saved complete starts (no new H/seed search)
×four arms: unchanged, target-beta wall, explicit artificial inversion drive,
drive+wall. Reference research-only, same native10 kJ/mol/nm/max1000×two phases,
only original positional restraints divided10 inphase2. No branchswap/reflection/
observation restore, native coordinates/RMSD, force-strength grid or mode selector.

Fixed dimensionless beta q=volume/(2.5 Å³), finitewall .5*10000*max(0,.4-q)^2,
drive .5*1000*(q+1)^2. Target CA–CB–OG1–CG2 order centresCB. One predeclared
wall targetfloor1.0 Å³, existing commongeometrynearplanarity gate0.5 unchanged.
Reference primitiveenergy/36 Cartesian derivative comparisons pass. Finiteforce
is NOT a mathematical chirality/floor guarantee, not a production-ready recipe.

Perarm12 enrolled: six execute (THR27/THR1/THR2 both old labels); two initialTHR6
CB S refuse, fourVAL nonapplicable (no beta stereocentre).48 declaredcells total,
24 executedcells/48 nativecalls, eight initialrefusals/16 nonapplicable, no blockers.
THR27 labelstarts duplicate: three targetmasks/two proteins/five distinct starts,
NOT six independent structures. Exact fullinitial heavy AND H for every arm.

Final CB R / commonraw+roundedgeometrypass: unchanged6/6 and6/6; wall6/6 and6/6;
drive0/6 and0/6; drive+wall6/6 but ONLY2/6 geometrypass. UnderdrivealltargetCA S,
beta flipsS withoutwall; beta remainsR withwall. This demonstrates a CB failure
mode independent of CA protection UNDER AN ARTIFICIAL DRIVE, NOT original unforced
THR2 bug replay/complete cause. Original failed output and nonexact replay remain.

Drive+wallfinal beta0.707–0.712 Å³, below requestedsoftfloor1.0: finitewall doesn't
enforcefloor exactly. THR2bothlabels pass; THR27/THR1bothlabels fail CA–CB–CG2 angle.
Allfour remain refused_geometry, no strength/budget/angle-gate adaptation after
results. AlreadywrongsourceTHR6 never "repaired" by preservingS/referencevolume.

Independent actual-frame abstract heavy-only views use captured externalbonds,
raw and ALL-coordinate PDBrounding; existing measuredchi targetgeometry plus
globalCA/THRbeta volume checks and RDKitCIP. Hgeometry not qualified; actual
original-coordinate motion retained, no observation-preserving-output claim.
WholeoriginalSystem equality verified with only newresearchforces removed; only
originalrestraintstrength changes inphase2, allmass/constraints/FFterms preserved.

Scripts `research_beta_wall.py`, `audit_beta_wall.py`, tests. Private
`tmp/pdbfixer-research/beta-wall-plan-01.json`, `beta-wall-run-01/`, originalaudit01,
freshreadback01. Publication independentlysaved-only byte-identical. Focused815
passed27.66s; not fullrepository/platform acceptance. Production/default/dependency/
library/public capability unchanged, no commit/push; PyRosetta licensinggate intact.

**Continuation:** do not call finite-wall signpass a full geometry/repair guarantee.
No posthocstrengthgrid or nativeaccuracy selection. A meaningful next safety step
needs explicit source/refusal policy plus independent finalgeometry validation,
not CA-only detection or a hidden THRbranch swap. Alreadywrongbeta repair remains
separate from protecting an admitted R centre; ILE/PTMs/otherchemistry unqualified.
Current finalverification complete:27 current/historical finalverifiers,55 fixture
hashes,1,251 local Markdown targets and thirteen documentation plans pass.
Generated CLI/GUI, agentmaps, whitespace, targeted Ruff and17-file package mypy
pass. Retained log: `tmp/pdbfixer-research/beta-wall-final-checks-01.log`.

## Earlier checkpoint: complete starting states separate two variability sources (2026-11-06)

User explicitly asked to investigate the current THR-Cβ issue first, WITHOUT
production changes. Read [complete-start controls](../research/legacy-start-controls.md),
[capture evidence](../research/legacy-start-capture-evidence.json) and
[replay evidence](../research/legacy-start-replay-evidence.json).

Frozen before execution: ALL six original prepared inputs/dat×labels1,2. Seed
Python random ONLY around native Modeller.addHydrogens, restore prior state.
Unchanged native minimize/default physics/auto platform/--no-solvent/standard.
Capture full System/Integrator/State XML, checkpoints, atom-map/bonds and actual
platform properties before/after each original native call. No native references/
accuracy, new FF/restraint/tolerance profile, seed expansion or production edits.

12/12 workflows,24 nativecalls/48 boundary snapshots, all OpenCL Apple M3 Pro
SINGLE precision. Crucial branch distinction: both1ENH masks keep H (zero Hcalls);
their seed labels are numerical repeats, NOT H interventions. Four1CRN masks do
eight total Hcalls. Acrosslabels heavy startcoords EXACT; named H RMS1.068–1.084 Å,
max2.022–2.048 Å. Named metrics include equivalent H permutations, not all physical
displacements. THR2 unique hydroxylHG1 difference0.134544 Å. All captured THR2
boundaries remainCB R;THR6 remainsS, no near-planarTHR1 endpoint in this control.

Original saved auditor incorrectly demanded heavycoords exactly unchanged across
addHydrogens; failed on1.283e−6 Å numerical movement. Original source/failure log
retained. Explicit `legacy-start-audit-amendment-01.json` binds revised auditor:
report exact:false/measured motion, NOT observation-preservation acceptance by
widening a gate. No scientific rerun, thresholds/seeds/coordinates changed.

Before any control results, also declared THR2+THR27 two labels×two repeats×two
modes=16 phase1 native replays: captured OpenCL checkpoint vs Reference XML start,
same serialized system/settings. Separate replay manifest binds ALL saved state.
16/16 complete; INITIAL HEAVY AND H coords EXACT for all16. All4 OpenCL repeat
pairs have nonzero heavy endpointdifferences0.005183–0.030307 Å, despite identical
starts/checkpoints; all4 Reference pairs repeat exactly includingH. THR27 label1/2
share same starting coords, so four declaredpairs are NOT four independentstarts.
Cross-platform heavy endpointdiff0.001981–0.079692 Å. Reference is research control,
NOT new production default/preferred scientific profile/general determinism claim.

ALL16 target CB R. Original THR2 R→S failure STILL NOT REPRODUCED/CAUSALLY EXPLAINED.
Original failed output remains immutable, not negated by these controls. Force
inventory identifies N–CA–C–CB signedvolume wall centeredCA, no target CB-centered
periodic improper. Other angle/bond/proper torsion/nonbonded terms still exist;
inventory is not a force-isolation/barrier certificate. Explicit CA invariant
does not also implement a signed-volume CB invariant. StartingTHR2 beta+0.198390 Å³
remains a poorly shaped/unprotected source condition, not proof of unique cause.

Scripts `research_legacy_start_controls.py`, `audit_legacy_start_controls.py`,
tests. Private `tmp/pdbfixer-research/legacy-start-controls-` plan01,capture01,
replay-plan01,replay01; `legacy-start-capture-audit-02.json`,
`legacy-start-replay-audit-02.json`, archivedoriginalauditor/failure and amendment.
Both publications fresh saved-only readback byte-identical. Focused810 passed28s,
not full repository/platform suite. Source/default/dependency/library/public
capability unchanged, no commit/push; PyRosetta production-license gate unchanged.

**Continuation:** retain original failed output, nonexact legacy replay and bounded
complete-start controls separately. Two variability sources now have direct case
evidence; neither establishes why original THR2 crossed. Do not grow a seed/profile
search until a wanted failure appears, restore observations, swap distinct THR
branches, or silently force Reference in production. Any next intervention should
be separately frozen and independently audit CB and CA geometry, not just detect
D-Cα or compare native RMSD. Current final verification complete:26 final verifiers,
55 historical fixture hashes,1,246 local Markdown targets and twelve documentation
plans pass. Generated CLI/GUI, agent maps, whitespace, targeted Ruff and17-file
package mypy pass. Retained log:
`tmp/pdbfixer-research/legacy-start-controls-final-checks-01.log`.

## Earlier checkpoint: logging-only legacy phase replay is nonexact (2026-11-06)

Read [phase trace](../research/legacy-stereo-trace.md) and
[native-free saved evidence](../research/legacy-stereo-trace-evidence.json).
Before execution froze ALL six original prepared inputs/.dat, one replay each,
same unchanged production minimize --no-solvent/standard naming/default physics/
auto platform. Supported Python Simulation.minimizeEnergy observer only, forwards
args/return unchanged and records ALL atom coordinates before/after each native
call. No forces/coordinates/seeds/settings changed, no native accuracy read.

All6 completed, two calls each→12 nativecalls/24 Context snapshots. Hooks restored.
ZERO6 final PDBs exactly repeat original full-comparator outputs; max heavy endpoint
differences0.047–2.130 Å retained, NOT hidden or rerun until failure. Initial Context
heavy identities match prepared inputs; numeric input differences1.28–3.34e−6 Å.
No exact-trajectory/observer-equivalence/automatic-platform determinism claim.

THR2 original preparedbeta +0.198390 Å³; tracedphase1 +2.432270,phase2 +2.440301,
roundedfinal +2.437493, CB R throughout sampled boundaries. Original comparator
CB S output remains immutable/verified; THIS REPLAY DOES NOT LOCALIZE THE ORIGINAL
R→S TRANSITION. H reconstruction/platform/observer numerical effects are not
separated. No repeat-until-match/failure or posthoc platform change.

Distinct THR1 outcome: Cα +0.403414 prepared→+0.080152 afterphase1→+0.080089 after
phase2→+0.082859 rounded. RDKit assigns CB R but NO CA CIP postphase/final; signed
volume stays positive, NOT D-Cα. Consistent with relative-volume guard design, not
causality from a force-isolation control or absolute tetrahedral-quality guarantee.
No negative CA at sampled boundaries, not whole-trajectory inversion certificate.

Scripts `research_legacy_stereo_trace.py`, `audit_legacy_stereo_trace.py`, tests.
Private prefix `tmp/pdbfixer-research/legacy-stereo-trace-`: plan01,run01,audit01,
fresh readback01. All source/snapshot/output hashes bound, native references never
read. Saved-only publication repeats byte-identically. Focused803 passed26.26s,
not full repository/platform suite. No production/default/dependency/library/public
capability change, commit/push; original PyRosetta production license gate intact.

**Continuation:** preserve this failed exact replay and original six comparator
outcomes separately. A meaningful next mechanism control needs declared starting
complete-state/system/H evidence and platform controls; do not infer a cause or
pick new numerical settings from known native errors. No branch swap/observation
restore or silent production fix. Current final verification is complete:25 final
verifiers,55 historical fixture hashes,1,239 local Markdown targets and eleven
documentation plans pass. Generated CLI/GUI, agent maps, whitespace, targeted Ruff
and17-file package mypy pass. Retained log:
`tmp/pdbfixer-research/legacy-stereo-trace-final-checks-01.log`.

## Earlier checkpoint: full legacy comparator exposes THR-CB failures (2026-11-06)

Read [full legacy comparison](../research/legacy-comparator.md) and
[actual-frame evidence](../research/legacy-comparator-evidence.json).
Separately frozen ALL six shared-admitted masks, one unchanged production CLI
prepare legacy→minimize --no-solvent run each, explicit standard naming/.dat,
unchanged physics/PROPKA/Reduce/default auto platform. All6 completed, no native-
error/stock-success filter or observation restore. Production sources/defaults/
dependencies/public capabilities remain unchanged; no commit/push.

Prepare original-heavy motion0, common2/6; three small POSITIVE target CA volumes
below research near-planarity bound, NOT D-Cα. Minimize common4/6, all target CA S,
but1ENH THR6/1CRN THR2 CB S rather than R. THR6 alreadyS afterprepare; THR2 R after
prepare→S afterminimize. Separate SIDECHAIN stereo issue, no complete mechanism
or production correction established. All original heavy identities retained;
minimize moves all at PDB precision, RMS0.120–0.216 Å/max0.710–1.306 Å. Heavy-only
actual-frame remasking preserves REAL relaxed coordinates, not original ones;
common observed_preserved refers only to that diagnostic view. H not covered.

Original-frame target RMSD median all6: prepare1.915500 Å, minimize1.163708 Å;
committed synthetic baseline same6 is0.191971 Å. Small correlated retrospective
two-protein comparison, no general ranking/new qualification. Earlier bare stock
matched-two counterexample remains. No settings/arm chosen from native errors.

Initial savedv1 reused fixed-observation helper for auxiliarychi, mixing original
anchors with relaxed additions. Retained unchanged. Explicit hash-bound
`amend_legacy_chi_reporting.py`/`legacy-chi-amendment-plan-01.json` correct ONLY
actual-framechi; all RMSD/motion/geometry/coverage/choices unchanged, no resampling.
Do not cite old hybridchi; publicationv2 binds correction/superseded hashes.

Scripts `research_legacy_comparator.py`, `audit_legacy_comparator.py`, amendment
and tests. Private `tmp/pdbfixer-research/legacy-comparator-plan-01.json`,
`legacy-comparator-run-01/`, `legacy-comparator-audit-01/`, originalreport01,
`legacy-comparator-readback-01/` retain raw logs/outputs/dat. Fresh saved-only
audit ledger,v1 and corrected publication repeat byte-identically. Focused suite
800 passed in26.25s; not full repository/platform acceptance.

Legacy final checks passed:24 final verifiers,55 historical fixture hashes,1,234
Markdown targets, ten completed documentation plans; generated CLI/GUI, agent
maps, whitespace, targeted Ruff and17-file package mypy pass. Log:
`tmp/pdbfixer-research/legacy-comparator-final-checks-01.log`.

**Resume:** separately predeclare actual THR-CB
transition diagnostics using frozen inputs/settings. No observation restore,
chemically distinct THR branch swap, production threshold changes or further
native-error-driven numerical grid. Preserve both final failures and movement.
Research-only PyRosetta/source/license/explicit production approval gate unchanged;
legacy itself does not depend on that library.

## Earlier checkpoint: prospective new-input comparison remains negative for means/refinement (2026-11-06)

Read [new-input comparison](../research/prospective-chi-guard.md),
[native-free selection ledger](../research/prospective-chi-guard-evidence-selection.json),
[accuracy report](../research/prospective-chi-guard-evidence.json) and
[ALL matched pairs](../research/prospective-chi-guard-pairs.json). No production/
library/dependency/fixture-registry/public-capability change, commit or push.
Installed PyRosetta remains research-only with mandatory original-table/source/
license/explicit user approval BEFORE production use of that route.

Before any new source download, froze1ENH/1CRN/1QYS, first MODEL/first encountered
ATOM chain, first two VAL/first two THR whole-sidechain masks, up to12 masks;
retain all shortfalls/refusals, NO replacement structures or source sanitization.
`scripts/prospective_chi_guard_cohort.py` and four synthetic tests freeze complete-
pool selection: ANY execution blocker blocks a pool; ANY guard exhaustion blocks
the whole guarded pool, never select only its successful starts. Stable index
ties, no native error or cross-profile best-arm selection.

Six declared arms: stock seeded PDBFixer Reference (NOT full legacy), synthetic
own fixed-ref2015, NCAC/library-mean prior, mean fixed-ref2015, unchanged mean-C,
mean bounded endpoint guard. Same frozen kernels/settings/gates; no new numerical
grid/angles/weights. Native accuracy only AFTER all arm choices commit. New
deliberate sources, not randomized/blinded/training-independent;1QYS Top7 is
Rosetta-designed and explicitly NOT independent of Rosetta design/training.

Artifacts inside `tmp/pdbfixer-research/`:
`prospective-chi-guard-enrollment-01/enrollment.json` (PRE-DOWNLOAD enrollment),
`prospective-chi-guard-inputs-01/` (all three official RCSB HTTP200 downloads,
URL/time/bytes/SHA256), `prospective-chi-guard-masks-01/` (frozen projection/masks).
`prospective-chi-guard-preparation-01/` completed through unchanged NCAC kernels.
Transport adapters were finished AFTER download, preserving the pre-download
scientific kernels/policy; do not call this a fully preimplemented/blinded study.
Library query01 FAILED before any API calls: view adapter eagerly imported the
OpenMM-dependent enrollment module into isolated PyRosetta. Original source/log
retained. Separate `research_prospective_guard_query.py`/query-plan02 fixes only
dependency isolation, no new dependency/kernel/angle/weight change or resampling.
Five complete three-well queries/periodic controls passed against same installed
StpDwn_0-0-0 binary; raw per-case means/weights/sigmas remain private.

10 masks (1ENH0VAL/2THR shortfall retained);1CRN4 masks;1QYS4 masks all share the
same missing/mismatched LINK/SSBOND-endpoint refusal in frozen projection. Do NOT
call originalTop7 deposition invalid or sanitize it into admission. Six admitted
requests, five mean pools:1CRN THR1 lacks observed phi/psi context, no inventedangle.
Only TWO proteins actually execute.30 NCAC stage snapshots/45 exact score imports
independently audited, fixed-score context identical within/across arms. Own
common14/18; means15/15; C/guard final15/15, protected motion0, stereo findings0.
Guard16 innercalls for15 starts, all15 actual-gradient criterion met, no exhaustion.

ALL six-arm choices committed BEFORE native accuracy; all15 unordered matched
pairs separately frozen before report. Selection counts stock2/10, own6/10,
mean-prior/fixed/C/guard5/10 each. All6 admitted native controls common-pass.
On SAME5 mean-supported masks, own exact/sym median0.135648 Å vs mean-prior0.283035,
mean-fixed0.278574, C/guard0.281341. No library-mean/refinement selected-output
rescue. C/G select same indices and byte-identical rounded heavy outputs5/5;
guardextra call is nonselected. Conversely SAME2 stock-supported VAL masks:
stock0.167449 vs own0.305992 Å. Keep that counterexample; coverage is not template
accuracy inferiority. Bare stock comparator is NOT production legacy prep/min/ZBS.
No numerical mode is adopted or weights/tolerances tuned from these native errors.

Adapter scripts: `research_prospective_guard_views.py`, explicit isolation
`research_prospective_guard_query.py`, all-arm independent
`audit_prospective_guard_cohort.py`, pre-report
`summarize_prospective_guard_pairs.py`. Private query/means/scores/C/guard trees
use the `prospective-chi-guard-` prefix. Key frozen artifacts:
`prospective-chi-guard-audit-manifest-01.json`, `prospective-chi-guard-selection-01.json`,
`prospective-chi-guard-report-01.json`, `prospective-chi-guard-pairs-plan-01.json`.
Fresh saved-only select/report readback must not alter kernels/choices or re-query
tables. Current full focused suite795 passed in27.04 s; not full repo/platform.
Seven-file routine plan `prospective-chi-guard-docs-plan.md`.

Known existing `dvbfixer` interpreter passes the required PDBFixer/Reference
preflight NOW (PDBFixer1.12.0/OpenMM8.5.2); fresh enrollment pins its source and
required template hashes. This is scoped readiness, NOT proof all historical
damaged resources were repaired or of their cause. No installation was changed.
Broad external environment-directory listing was denied and not bypassed; the
already-used interpreter suffices, no wider environment discovery is necessary.
Original new source bytes are retained privately inside the project, not yet
registered as tracked fixtures; do not count them as production fixture acceptance.
Installed PyRosetta remains research-only, no tables/production capability/push.

Prospective publication verification is complete: all23 final verifiers,55 historical
fixture hashes,1,228 local Markdown targets and nine documentation plans passed.
All three saved-only publication readbacks are byte-identical; generated CLI/GUI,
agent maps and whitespace checks passed. Log:
`tmp/pdbfixer-research/prospective-chi-guard-final-checks-01.log`.

**Resume:** separately freeze the full legacy comparator. Further numerical
grids on known native errors are inappropriate. Scientifically meaningful next
work is a separately declared production legacy prepare→minimize comparator or
source-context/broader chemistry qualification, NOT a claim stock2/6 proves legacy
inferiority. Actual legacy relaxation may move observations: measure/report it,
never restore them into a supposedly observation-preserving result. Every dev
minimize/ZBS call MUST include --no-solvent. Known interpreter has PROPKA and
Reduce executable available with its bin directory prepended in CHILD PATH;
no environment install/mutation or broad directory listing is required.

## Earlier checkpoint: bounded actual-endpoint guard passes39/42 gradients,37/42 combined (2026-11-06)

Read [bounded endpoint-gradient guard](../research/chi-endpoint-guard.md) and
[independent evidence](../research/chi-endpoint-guard-evidence.json).
Scripts `research_chi_endpoint_guard.py`, `audit_chi_endpoint_guard.py`, matching
tests. No production/library/fixture/default/dependency/public capability change,
commit or push. Research-only PyRosetta and mandatory original-table/provenance/
license/explicit user production approval gate unchanged.

Frozen before execution: ALL42 original mode-C starts/all15 pools/seven inherited
failures, no native/error/gradient/geometry filtering. Same inner mode-C L-BFGS,
absolute1e−6/gmax1e−4 perDEGREE/max200 each call; actual returned-pose gradient
measured on independent pose/ScoreFunction clones. Stop if actualmax<=gmax,
otherwise reenter identical inner optimizer, maximum THREE TOTAL CYCLES. No new
starting angles, weights, threshold grid, adaptive budget or candidate ranking.

All42 first cycles reproduce original complete mode-C returned atoms/raw scores
EXACTLY. Four starts continue;49 inner calls total/seven extras, zero execution
blockers. Final actual-gradient criterion39/42; three budget exhaustions retained:
THR22/i2 final0.000124301;THR27/i1 0.000126400;THR76/i2 0.000494669 perdegree.
VAL56/i0 reaches2.814014e−6 aftercycle2. No posthoc budget increase or convergence
label on the three exhausted cases. More cycles are NOT silently authorized by
these results. Endpoint threshold passes themselves are not global/general
convergence/tolerance-independent-coordinate certificates.

Independent EVERY-cycle raw/rounded common audit passes47/49 phases; final40/42,
same rejected THR55/i0 and VAL56/i2. Combined final gradient+geometry37/42.
Protected motion0.0 Å, all49 CA S/THR CB R, zero global stereo findings. Geometry
is audited afterward, not used to choose continuation. Lower raw H-containing
score is not native accuracy; no native coordinates/RMSD or profile selection.
No observation restore, reflection, context repair, scientific output ranking/PDB
publication or production integration. Original mode-C results remain unchanged.

Private plan/run: `tmp/pdbfixer-research/chi-endpoint-guard-plan-01.json`,
`chi-endpoint-guard-run-01/`. All cycles retain complete atoms/scores/gradient/
decision. Focused suite782 passed in27.43 s; not full repository/platform acceptance.
Four-file routine plan `chi-endpoint-guard-docs-plan.md`.

**Resume:** keep the three exhaustions and separate numerical criterion from
geometry and native accuracy. Do not grow an adaptive restart/tolerance grid or
pick numerical settings using known native errors. Numerical mechanism evidence
now supports pre-step-gradient stopping and the bounded guard helps one case but
does not qualify a scientific rebuilder. Further scientifically meaningful work
should separately freeze independent evaluation or missing-context/chemistry/full-
legacy comparison, retaining observation ownership and every failure. Original-
table/license/user production gate remains mandatory, not a blocker to authorized
local research. Endpoint-guard readback/quality is complete: byte-identical saved-
only readback, all22 final verifiers,55 fixture hashes,1,218 local Markdown targets
and eight current documentation plans pass. Generated CLI/GUI, agent maps,
whitespace, targeted Ruff and17-file package mypy pass. Retained checks/quality:
`tmp/pdbfixer-research/chi-endpoint-guard-final-checks.log`,
`chi-endpoint-guard-docs-quality.md`, `validate-chi-endpoint-guard-final.py`.

## Earlier checkpoint: native checkpoints support pre-step-gradient stopping (2026-11-06)

Read [native χ checkpoints](../research/chi-checkpoints.md),
[quiet norm evidence](../research/chi-checkpoints-evidence.json), and
[verbose component evidence](../research/chi-checkpoint-components-evidence.json).
Scripts `research_chi_checkpoints.py`, `research_chi_checkpoint_verbose.py`,
`audit_chi_checkpoints.py`, `audit_chi_checkpoint_components.py`; matching tests.
No production/library/fixture/default/dependency/public-capability changes, commit
or push. Installed PyRosetta remains research-only, mandatory original-table/
provenance/license/explicit user production gate unchanged.

Supported native AtomTreeMinimizer.deriv_check_result, NOT Python C++ virtual
callbacks/empty poses. Frozen ALL42 original mode-C starts/all15 pools/seven
inherited failures; no endpoint-gradient/native/geometry filtering. Same numerical
settings/MoveMap/objective/weights. Quiet native checker uses its fixed signed
.001/.002/.004/.008/.016° probes; all42 original initial/before/after/result bytes
and scores repeat EXACTLY.1097 native derivative checkpoints; protected movement0.

Native individual DOF rows only populate above0.001 analytic-or-numeric magnitude;
invalid defaults are ABSENT, never zero gradients. L2 norm is not max component
for two chi. THR22 last norm0.0001033724 by itself was inconclusive, not a pass.

Separate frozen verbose/stdout amendment keeps ALL42 starts and verifies all16
exposed numerical option getters unchanged. Checked native flushes/synthetic
real-tracer BEGIN/END controls in every capture; adapters restore. Actual native
diagnostic flags in verbose-adapter.json override the quiet outer wrapper's
recorded requested flags. Again42/42 endpoints/scores EXACT and all numerical
checkpoint files identical to quiet run.1097 complete vectors/1957 component
values validated against native full-precision L2 norms; printed-value uncertainty
kept separate from scientific thresholds. No full function/iteration/return trace.

Last computed components per DEGREE at the two formerly unexplained endpoints:
2CI2 THR22/i2 [-9.55821122e−5,-3.93689471e−5],20 checkpoints, versus actual endpoint
max0.001844293;1CTF VAL56/i0 [5.66662531e−6],9 checkpoints, versus endpointmax0.007424845.
Both last vectors are belowgmax1e−4 but final gradients are not. The two explicit
cap cases THR27/i1 and THR76/i2 each have201 checkpoints and above-cutoff last
components, consistent with prior max200 warnings.

Pinned Minimizer.cc computes Gmax BEFORE line minimization, then checks that
previous gradient after the line step (retained source674–744). Source+exact native
checkpoints support PRE-STEP-GRADIENT stopping for the two no-warning cases, not
missing gradients, an internal Rosetta unit error or moved observations. Still
NOT a directly recorded return-branch ID or full function trace/certificate.

Private plans/runs: `tmp/pdbfixer-research/chi-checkpoints-plan-01.json`,
`chi-checkpoints-run-01/`, explicit `chi-checkpoints-plan-02.json`,
`chi-checkpoints-run-02/`. Native MinMover source/provenance pinned under
`minmover-source-1a56185/`. Focused suite774 passed in26.71 s; not full repository/
platform acceptance. Five-file routine plan `chi-checkpoints-docs-plan.md`.

**Resume:** separately freeze a bounded ACTUAL-endpoint-gradient guard experiment
on all42 starts before any numerical control changes. Do not label the last
pre-step gradient a final stationary certificate. Preserve explicit caps versus
pre-step stopping evidence, all earlier profiles/geometry failures; no native-RMSD
choice or adaptive optimizer grid. Parent common40/42/output bytes unchanged.
Broader chemistry/context/full-legacy/independent accuracy gates remain open.
Native-checkpoint readback/quality is complete: both publications byte-identical,
all21 final verifiers,55 fixture hashes,1,212 local Markdown targets and seven
current documentation plans pass. Generated CLI/GUI, agent maps, whitespace,
targeted Ruff and17-file package mypy pass. Retained checks/quality:
`tmp/pdbfixer-research/chi-checkpoints-final-checks.log`,
`chi-checkpoints-docs-quality.md`, `validate-chi-checkpoints-final.py`.

## Earlier checkpoint: unit-aware analytic derivatives leave two unexplained endpoints (2026-11-06)

Read [unit-aware derivatives](../research/chi-derivatives.md) and
[independent evidence](../research/chi-derivatives-evidence.json). New private
`research_chi_derivatives.py`, `audit_chi_derivatives.py`, matching tests. No
production/library/fixture/default/dependency/public-capability change, commit or
push. Research-only installed PyRosetta and mandatory original-table/provenance/
license/explicit user production approval gate unchanged.

Pinned installed-revision AtomTreeMinimizer/AtomTreeMultifunc/MinimizerMap and
atom_tree_minimize.cc source: minimizer PHI variables are DEGREE-scaled (180/pi),
their derivatives score/DEGREE. Earlier finite slopes are correctly per RADIAN;
do not directly compare those numbers with gmax1e−4. Equivalent per-radian cutoff
is0.00572958. Finite-descent floor1e−7 energy was never a convergence tolerance.

Frozen ALL42 mode-C endpoints/all15 pools/seven inherited failures, no gradient/
native/geometry filtering. Exact original full-pose/score replay; saved mobile-
only endpoint materialization remains exact, energy recomputation bound unchanged.
Clone-only same MoveMap/MinimizerMap/setup_for_minimizing; call AtomTreeMultifunc
dfunc directly, NO Minimizer.run. Verify all62 CHI/PHI target/vector identities and
180/pi scale. Fixed signed .01/.005° variable probes,248 atom traces, no new grid.

All42 audited,0 blockers;62 analytic gradients/124 central-difference comparisons.
Protected movement0.0 Å. External set_chi versus minimizer-vector numeric slopes
match0/124 mismatches after unit conversion. Analytic/numeric reporting bound
frozen1e−6+1e−3*max(abs(a),abs(n)) perdegree;1/124 mismatch retained:1R69 VAL24/i0
chi1 at .01°, error1.1004971e−6 vs bound1.0590693e−6. At .005° error2.7511073e−7
passes. Do not drop coarse-step failure, tune bound or declare all derivatives
validated; finite-step sensitivity alone does not prove its cause.

Only FOUR endpoint analytic maxima exceed gmax1e−4 perdegree:2CI2 THR22/i2
0.001844293;1R69 THR27/i1 0.000296561;1CTF VAL56/i0 0.007424845;1CTF THR76/i2
0.000951026. Middle THR27 and last THR76 have explicit max200 warnings.15 of19
finite-descent endpoints already lie at/below optimizer gradient cutoff once
units agree. Remaining TWO above-cutoff/no-warning endpoints THR22/i2 and VAL56/i0
still unexplained; endpoint derivative values are not a recorded return branch.

Private source/provenance: `tmp/pdbfixer-research/atomtree-source-1a56185/`,
`atomtree-derivative-source-1a56185/`; no tables downloaded. Plan/run:
`chi-derivatives-plan-01.json`, `chi-derivatives-run-01/`. Focused suite761 passed
in26.76 s, not full repository/platform acceptance. Four-file routine plan
`chi-derivatives-docs-plan.md`.

**Resume:** a separately frozen original-start L-BFGS call/gradient trace may
resolve the two remaining returns; preserve all42 starts, settings/objective and
require exact endpoints. Do not infer return cause from final gradient alone or
change cap/tolerances/weights using native error. No stationary/global optimum,
probe geometry/CIP or scientific accuracy qualification. Broader chemistry,
context support/full-legacy comparison remain open. Derivative readback/quality
is complete: byte-identical saved-only readback, all20 final verifiers,55 fixture
hashes,1,204 local Markdown targets and six current documentation plans pass.
Generated CLI/GUI, agent maps, whitespace, targeted Ruff and17-file package mypy
pass. Retained final evidence/quality:
`tmp/pdbfixer-research/chi-derivatives-final-checks.log`,
`chi-derivatives-docs-quality.md`, `validate-chi-derivatives-final.py`.

Subsequent call-trace API feasibility controls are NOT cohort executions:
`tmp/pdbfixer-research/chi-calltrace-synthetic-01.log` shows the bare Multifunc
binding refuses subclass construction; an empty-pose-backed AtomTreeMultifunc
callback control in `chi-calltrace-synthetic-02.log` terminated SIGSEGV before a
usable trace. Do not infer a real-structure defect or callback dispatch behavior
from that unsupported empty setup, or repeat it. No installed/production state
changed; all42 derivative evidence remains intact. See
`tmp/pdbfixer-research/chi-calltrace-api-feasibility.md`. A supported and validated
instrumentation route is still required before any claimed full call trace.

## Earlier checkpoint: logging-only replay identifies two max200 stops (2026-11-06)

Read [termination diagnosis](../research/chi-termination.md) and
[independent evidence](../research/chi-termination-evidence.json). New private
scripts `research_chi_termination.py`, explicit native-flush amendment
`research_chi_termination_flush.py`, independent `audit_chi_termination.py`;
synthetic diagnostic/flush tests. Research-only permission and mandatory original-
table/provenance/license/explicit user production approval gate remain unchanged.
No production/fixture/dependency/library/default changes or commit/push.

Frozen before replay: ALL42 original mode-C starts/all15 executable pools/seven
inherited nonexecutions; no finite-descent/native/geometry filtering. Same numerical
settings/objective/weights/ownership; only add tracer unmute for exact upstream
Minimizer channel. Require initial snapshot, before/after atom traces AND complete
result bytes/numerical scores identical to mode C. Both attempts reproduce42/42.

Original `chi-termination-run-01/` capture FAILED per-case attribution: native
buffering emitted two warnings only after fd captures ended. Its zero per-case
counters must NOT be interpreted as no warnings. Preserve outer log and evidence;
auditor refuses unamended attribution. Synthetic first attempt failed unsupported
TracerProxy << str binding; second real-tracer write/native-fflush control passed.

Separate frozen `chi-termination-plan-02.json`/`chi-termination-run-02/` adds checked
native-stdio fflush at fd boundaries plus real tracer synthetic BEGIN/END controls
in every capture. Missing controls/failure blocks; synthetic controls contain no
optimizer-warning text. All42 control pairs verified, both adapters restored,
endpoints/results EXACT. No numerical changes or broader runtime-logging patch.

Controlled native logs identify TWO explicit LBFGS MAX CYCLES200 NOT CONVERGED:
1R69 THR27 index1 and1CTF THR76 index2. Zero explicit line-search-abort messages.
Both are among19 finite-descent endpoints; this is saved-only full-inventory
association, not new filtering or physical causality. Other17 finite-descent
endpoints still unexplained. Forty cases without explicit messages are NOT a
normal-convergence certificate; unlogged return paths exist. No analytic derivative
or full iteration trace, increased cap, native RMSD or preferred-profile selection.

Focused suite756 passed in26.37 s, not full repository/platform acceptance.
Private artifacts under `tmp/pdbfixer-research/`: both frozen plans, both runs,
outer logs and synthetic diagnostics. Four-file routine plan
`chi-termination-docs-plan.md`. Parent common40/42/results remain unchanged.

**Resume:** preserve the two explicit cap events separately from17 unexplained
finite-descent endpoints. Next numerical-mechanism study must separately freeze
analytic-derivative or full iteration/return-path instrumentation on ALL original
starts, without native tuning or an adaptive optimizer grid. More iterations are
not automatically a remedy for the other17. Independent scientific accuracy,
missing context, broader chemistry and full-legacy comparison remain open.
Termination publication readback/quality is complete: byte-identical saved-only
readback, all19 final verifiers,55 fixture hashes,1,197 local Markdown targets and
five current documentation plans pass. Generated CLI/GUI, agent maps, whitespace,
targeted Ruff and17-file package mypy pass. Retained final evidence/quality:
`tmp/pdbfixer-research/chi-termination-final-checks.log`,
`chi-termination-docs-quality.md`, `validate-chi-termination-final.py`.

## Earlier checkpoint: finite local probes find residual descent at 19/42 endpoints (2026-11-06)

Read [finite χ stationarity diagnosis](../research/chi-stationarity.md) and
[independent evidence](../research/chi-stationarity-evidence.json). New private
scripts `research_chi_stationarity.py`, `audit_chi_stationarity.py`, matching tests.
Installed PyRosetta remains research-only; original-table/provenance/license/user
production approval gate unchanged. No production/fixture/dependency changes,
public capability, new optimizer/profile selector, commit or push.

Frozen BEFORE probes: ALL42 mode-C endpoints, all15 executable pools, including
both parent geometry rejections; seven inherited nonexecutions retained. Every
target chi including proton chi; independent clones at ±.01° and ±.005°. No
native inputs, RMSD, adaptive steps, filtering by motion/error or optimization.
Mode C is the last predeclared numerical control, not an accuracy-selected mode.

All42 initial complete snapshots/numerical scores replay EXACTLY. Materialize
saved endpoint by set_xyz ONLY mobile atoms; every endpoint atom matches saved
trace exactly, topology unchanged. Protected observations are never moved then
restored. Endpoint energy is recomputed within separately declared absolute1e−8
bound, actual max difference5.68e−14. Do NOT call this exact final-score replay or
another minimizer execution. Every signed probe retains complete private atom
trace; independent ownership readback checks all248 traces, protected motion0.0.

All42 endpoints audited, zero blockers, 62 target chi/248 signed probes. At the
predeclared −1e−7 score-change reporting floor,19/42 endpoints have finite descent;
14/42 have descent on a heavy-moving chi (overlapping counts). Largest absolute
central slope0.425411 score/rad; best signed delta−7.25595e−5. Both step sizes,
central slope/curvature and step discrepancy retained. The other23 are NOT
certified stationary. No analytic-gradient comparison or termination trace;
iteration limit/line search/derivative/parameterization cause remains unknown.

Probes are energy diagnostics, NOT individually geometry/chirality-qualified
proposals; no raw/rounded common gate run for them and none is selected/output.
Parent mode-C common40/42 remains unchanged. No native accuracy/profile choice,
weight/angle/tolerance changes or global/local convergence certificate.

Private plan/run: `tmp/pdbfixer-research/chi-stationarity-plan-01.json`,
`chi-stationarity-run-01/`. Focused suite751 passed in26.38 s; targeted Ruff and
17-file package mypy pass, not full repository/platform acceptance.

**Resume:** diagnose termination only via a separately frozen native-free
iteration/line-search or analytic-derivative control on the complete inventory.
Do not infer the cause from these finite slopes or choose new numerical settings
by known native accuracy. Preserve all previous modes/results and failures.
Independent accuracy, missing context, broader chemistry and full-legacy comparison
remain separate open gates. Stationarity publication readback/quality is complete:
byte-identical independent readback, all18 final verifiers, 55 fixture hashes,
1,190 local Markdown targets and four current documentation plans pass. Generated
CLI/GUI, agent maps and whitespace pass. Retained checks/quality:
`tmp/pdbfixer-research/chi-stationarity-final-checks.log`,
`chi-stationarity-docs-quality.md`, `validate-chi-stationarity-final.py`.

## Earlier checkpoint: predeclared stopping controls change outputs, no native profile choice (2026-11-06)

Read [stopping-policy controls](../research/chi-stopping-controls.md) and
[independent native-free evidence](../research/chi-stopping-controls-evidence.json).
Scripts `research_chi_stopping_controls.py`, `audit_chi_stopping_controls.py`;
matching tests. Same research-only PyRosetta authorization, original-table/license/
explicit user production gate. No production/library/fixture/default/dependency
changes or new candidate starts, no commit/push.

Pinned public Rosetta Minimizer.cc to the exact reported revision1a56185: fractional
energy criterion for original lbfgs_armijo_nonmonotone, absolute variant _atol;
L-BFGS also uses gradient cutoff. Installed default gradient cutoff verified1.0.
Source/provenance privately retained at `tmp/pdbfixer-research/minimizer-source-1a56185/`.
No original rotamer table downloaded, no licensed data moved into production.

BEFORE controls, froze ALL42 original starts/15 pools, seven inherited failures,
max200 fixed, same objective/weights and mobility. Ordered controls: A absolute
tol.01/gmax1; B absolute tol1e−6/gmax1; C absolute tol1e−6/gmax1e−4. Baseline→A
criterion only, A→B tolerance only, B→C gradient cutoff only; sequential design,
NOT full factorial. No native-error/movement/score-margin filtering, no best-mode
selector. All126 initial snapshots/numerical scores repeat exactly; runtime getter
checks verify every minimizer's planned configuration and research hooks restore.

All126 return, zero execution/mobility blockers; protected motion0.0 throughout,
CA S/THR CB R, no global stereo findings. Each mode common40/42, selects15 pools.
A selected indexchanges0, rounded outputchanges1 vs baseline; B indexchanges1,
roundedchanges12; C indexchanges1, roundedchanges15. The changed initial index in
B/C is 1PGB THR11 0→2, but NO native coordinates/accuracy files read or metrics
computed/reported by this suite. Do not infer its accuracy or choose a mode by it.
These are126 repeated starts, not independent scientific samples/cohorts.

Raw energy drops vs baseline up to19.370642 in A,29.269935 in B/C. Max generated-
heavy movement1.553417 Å in A,2.550827 in B/C, no observed movement. A rejects
THR22 i2/THR55 i0; B/C reject VAL56 i2/THR55 i0. All unchanged-gate rejections
retained. B→C still changes four selected rounded structures and max per-state
score0.052735 units although selected indices stable. No convergence/global-
minimum/tolerance-independent coordinate certification; lower score not accuracy.

Private suite `tmp/pdbfixer-research/chi-stopping-plans-01/`, three
`chi-stopping-<id>-run-01/` trees. Expanded focused tests742 passed in26.95 s,
not full repository/platform acceptance. Six-file routine plan `chi-stopping-docs-plan.md`.

**Resume:** retain all modes; do not choose a profile by native RMSD or grow an
adaptive tolerance/angle grid on these known cases. A separate derivative/iteration
or mathematically declared convergence control may clarify termination; none is
certified here. Native accuracy qualification, missing-context support and full
legacy comparison remain independent open gates. Production remains disabled
until user-supplied original tables/rights/explicit approval AND scientific gates.
Final native-free control checks are complete: all 17 final verifiers, 55 fixture
hashes and 1,184 local Markdown targets pass; all three documentation plans are
complete. Both publication readbacks are byte-identical. Focused tests742,
targeted Ruff, 17-file package mypy, generated CLI/GUI, agent maps and whitespace
pass. See `tmp/pdbfixer-research/chi-stopping-docs-quality.md`.

## Earlier checkpoint: bounded generated-only χ relaxation does not rescue choices (2026-11-06)

Read [generated-only χ relaxation](../research/generated-chi-relax.md),
[independent selection ledger](../research/generated-chi-relax-evidence-selection.json)
and [accuracy evidence](../research/generated-chi-relax-evidence.json).
New script `research_generated_chi_relax.py`, independent auditor
`audit_generated_chi_relax.py`; matching synthetic mobility/audit tests. Research-
only PyRosetta permission and mandatory original-table/license/user approval
before production remain unchanged. No source/fixture/installed-library/default/
dependency/public-capability changes, no commit/push.

Protocol frozen before minimization: ALL 42 feasible mean starts from 15 pools,
seven inherited nonexecutions retained in the original 22-mask denominator.
Every complete initial fixed-heavy snapshot/numerical score repeats exactly.
Enable only target χ, all BB/jumps/other χ disabled; only generated distal heavy
atoms and target sidechain H mobile. Generated CB, source observations, backbone
H, context H, extras and virtuals protected. No packing, Cartesian minimization,
reflection, context repair, source rewrite or coordinate restore. Same installed
library hash/weights. lbfgs_armijo_nonmonotone, tol0.01, max200, neighbor-list off;
NOT a convergence/global-minimum claim or native-tuned optimization protocol.

All 42 bounded refinements return, zero execution/mobility blockers. Actual
protected motion0.0 Å everywhere (declared numerical bound1e−9 Å). All raw CA S,
all THR CB R, no global stereo findings. Independent raw/rounded common gate
passes41, rejectsone nonselected 2CI2 THR55 index0 for OG1–I/52 O clash despite
lower score. Only three nonselected starts have substantial heavy motion:
1CTF THR76 i1 0.276725 Å; 1R69 VAL6 i2 0.152352; THR55 i0 0.778725.

All fifteen selected indices and rounded heavy PDB outputs remain EXACTLY the
parent means' choices. Old-eleven selected exact/sym median0.198491 Å; new-four
exact1.930166/sym1.463808. No selected-output accuracy rescue under this bounded
configuration; do not generalize to all possible refiners or claim convergence.

Score semantics are explicit: final RAW refined pose with target sidechain H,
NOT fresh rounded-heavy PDB reimport scoring. Independent raw/rounded geometry
must pass before ranking; any execution blocker makes a pool incomplete and
blocks strict selection. Complete-pool ledger commits BEFORE native reporting.
This is coupled generated-heavy/sidechain-H relaxation, not an isolated H-only
causality test. No raw library parameters/tables in tracked evidence.

Private protocol/run: `tmp/pdbfixer-research/generated-chi-relax-plan-01.json`,
`generated-chi-relax-run-01/`. Expanded focused suite738 passed in27.68 s, not full
repository/platform acceptance. Routine six-file plan `generated-chi-relax-docs-plan.md`.

**Resume:** if diagnosing the very small motion, separately predeclare convergence/
derivative controls, with no native-error tuning of tolerance/angles/weights.
The current fixed200/tol0.01 result is a bounded negative experiment, not proof
of a converged energy minimum. Alternatively pursue genuinely new input validation
or separately scoped context completion/full legacy comparator. Production remains
disabled pending original-table/license/user gate AND scientific qualification.
Final χ-relax checks complete: both publications repeat byte-identically, all
16 old/new final verifiers and 55 fixture hashes pass; 1,178 local Markdown
targets valid. Ruff, 17-file package mypy, generated CLI/GUI checks, agent maps
and whitespace pass. Final verifier/quality:
`tmp/pdbfixer-research/validate-generated-chi-relax-final.py`,
`generated-chi-relax-docs-quality.md`.

## Earlier checkpoint: cross-arm context matches; model-pair preferences localized (2026-11-06)

Read [context/pair diagnosis](../research/mean-context-pairs.md),
[all-pool comparability evidence](../research/mean-context-comparability-evidence.json)
and [independent pair attribution](../research/mean-pair-energies-evidence.json).
New scripts `audit_mean_comparability.py`, `research_mean_pair_energies.py`,
explicit JSON-only replay adapter `research_mean_pair_replay.py`, independent
`audit_mean_pair_energies.py`; tests use matching names. No production/library/
fixture/scientific-threshold/dependency change, no commit/push. PyRosetta remains
research-only with mandatory original-table/license/user approval before production.

All 22 original masks retained, fifteen executable pools/seven inherited failures.
Read 84 mean+relative scored snapshots (42 each). All 15 pools match exactly
across arms in source bytes/ownership/target, observed+extra heavy coordinates,
non-target H, virtual coordinates, residue/atom type inventory and graph metadata.
This closes the earlier cross-arm context caveat on these data, not chemical
context/protonation calibration. No new score or selector in that comparability
audit. Target generated atoms and target H correctly excluded from context.

Generated-only VAL symmetry is NOT a rescue: new-four prior exact/sym median
1.210568/0.879635 Å, score 1.930166/1.463808. VAL29 score exact1.770814→sym1.499963
still much worse than prior0.078348; VAL21 exact2.336914→sym1.427653 still poor.
Old-eleven medians0.198491 unchanged by symmetry. No renaming/selection change.

Pair plan enrolled ALL two prior/score disagreements, not native-error filtering.
Four original states repeat with exact numerical scores/terms and complete
compressed pose snapshot bytes before EnergyGraph readout; all atoms stay fixed.
Initial two runs failed metadata container equality for terminal OXT (live tuple
vs JSON list), NOT coordinate/numeric mismatch. Both failures and diagnostics
retained. Separate frozen v2 JSON transport normalizes only containers; no added
numeric, snapshot or coordinate tolerance. Research hook restores.

Independent graph sums close for fa_atr/fa_rep/fa_sol plus fa_dun onebody. Explicit
cached-energy bookkeeping bound 1e−5+1e−7*sum(abs(terms)); actual delta residuals
~1e−13. For these three pair terms, zero changed non-target edges; fa_dun delta
only target. This is recorded model attribution, not physical/native causality.
VAL29 attraction/solvation gains chiefly PHE30/ALA26/THR18 outweigh repulsion and
prior penalty. THR26 repulsion relief chiefly VAL24/PHE44/ARG43 overrides prior
usefully. Do not blanket-upweight prior based on VAL29 or claim missing fa_dun.

Private frozen/readback artifacts under `tmp/pdbfixer-research/`:
`mean-comparability-plan-01.json`, `mean-pair-energies-plan-01.json`, failed
`mean-pair-energies-run-01/` and `-02/`, explicit `mean-pair-energies-plan-02.json`,
completed `mean-pair-energies-run-03/`. Individual library parameters remain private.
Expanded focused suite: 731 passed in26.12 s, not full repository/platform suite.
Routine seven-file plan: `mean-context-pairs-docs-plan.md`.

**Resume:** next causal intervention may be a separately frozen generated-only
target-chi relaxation, preserving CB/backbone/all observed atoms/context, auditing
every raw/rounded output and global chirality and retaining failures. It is NOT
implemented or qualified here. No new weight/angle tuning on these known native
errors. Context completion, genuinely independent accuracy validation and full
legacy comparator remain open. Finish final saved-context/pair readback/quality
before starting new intervention. Those checks are now complete: both evidence
files repeat byte-identically, all 15 old/new final verifiers and 55 fixture hashes
pass, 1,171 local Markdown targets valid. Targeted Ruff, 17-file mypy, generated
CLI/GUI, agent maps and whitespace pass. Final verifier/quality:
`tmp/pdbfixer-research/validate-mean-context-pairs-final.py`,
`mean-context-pairs-docs-quality.md`.

## Earlier checkpoint: library-mean proposals preserve geometry, new-input scoring worsens

Read [library-mean proposal results](../research/dunbrack-means.md),
[sanitized independent selection ledger](../research/dunbrack-means-evidence-selection.json)
and [accuracy evidence](../research/dunbrack-means-evidence.json). New script
`research_dunbrack_means.py`, auditor `audit_dunbrack_means.py`, synthetic tests
`test_dunbrack_means.py`. Research-only installed PyRosetta authorization and
mandatory original-table/license/user approval before production remain intact.

Frozen BEFORE new generation: all 22 original masks, fifteen executable parent
NCAC pools/seven inherited nonexecutions. Three API mean-chi proposals per
executable case, one per well. Rotate both generated distal branches rigidly
around saved CA–CB, preserving CB/backbone/observations. No minimization, packing,
sigma sampling, native-tuned means/weights or gate change. This is retrospective
on the same already evaluated cohorts, not a new input-validation cohort.

Executed 45 mean proposals, 42 common-feasible and three retained clashes
(2CI2 VAL32 wells 1/3, THR55 well 2), 42 exact fixed-heavy imports. Both prior-only
and unchanged ref2015 selectors commit before native readback. Auditor verifies
rigid transport/actual mean chi, source bytes/fixed CB, independent common audits,
complete feasible pool, physical/target/external graph import and score terms.
Scorer-owned context/extras/virtuals are identical within every mean pool; cross-
arm equality versus parent relative pools is not separately established here.
All selected outputs satisfy L-Cα/THR-CB and unchanged common gate, not accuracy.

Old-eleven mean prior/scored medians both 0.198491 Å vs parent relative scored
0.226284. But 1R69 THR26 differs: prior 2.022584 vs score 0.240256 Å. New-input
four-executable median prior/score 1.210568/1.930166 vs parent relative score
1.222795 Å. On 1PGB VAL29, score rejects the prior's good 0.078348 Å mean proposal
in favor of 1.770814 Å; THR11 remains poor (2.089519 Å) under both. Negative
generalization retained; neither family prior nor total ref2015 is qualified.
No table rows/means/sigmas/individual weights are in tracked publications.

Private artifacts: `tmp/pdbfixer-research/dunbrack-means-plan-01.json`,
`dunbrack-means-generation-01/`, `dunbrack-means-score-run-01/`. Focused suite
718 passed in 26.49 s, not full repository/platform acceptance. No installed
library, source fixture, production dependency/consumer/public capability change.
No commits/pushes. Six-file routine docs plan: `dunbrack-means-docs-plan.md`.

**Resume:** inspect retained score-term/context contributions for the discordant
THR26/VAL29 cases without refitting weights or thresholds to their native errors.
Do not densify sigma/angle proposals to manufacture a positive claim. Context
completion, independent evaluation and full legacy comparator remain separate
open gates; production remains forbidden pending the explicit user/materials gate.
Final mean-study checks complete: both sanitized publications replay byte-
identically; all 14 final verifiers, 55 fixture hashes and 1,163 local Markdown
targets pass. Targeted Ruff, generated CLI/GUI, agent maps and whitespace pass;
unchanged source package retains its passing 17-file mypy check.
Final verifier/quality: `tmp/pdbfixer-research/validate-dunbrack-means-final.py`,
`dunbrack-means-docs-quality.md`.

Finite saved-score diagnosis is now also complete: new
`scripts/analyze_dunbrack_mean_terms.py` analyzes ALL two prior/score disagreements,
no native-data enrollment filter, new scoring or weight modification. VAL29 total
selected-minus-prior −0.404808 Rosetta units: fa_atr −1.597526 and fa_sol −0.867703
outweigh unfavorable fa_dun +1.012951 and fa_rep +0.790042. THR26 total −9.906042,
dominated by fa_rep −9.448179 despite fa_dun +1.023393. Terms reconstruct totals;
this explains model ranking, not physical causality or a weight-refitting recipe.
Private `dunbrack-means-term-deltas-01.json`; two synthetic arithmetic tests pass.
Next work must retain these conflicting context/prior cases rather than blindly
upweighting the prior or claiming missing fa_dun explains the score failure.
Latest expanded focused suite including score-term arithmetic: 720 passed in
26.88 s. Independent mean/prior verifiers, agent maps, whitespace and note/handoff
links rechecked after the diagnosis update.

## Earlier checkpoint: installed-library access works, prior-only ranking does not rescue accuracy

Read [installed-library prior study](../research/dunbrack-prior.md),
[sanitized selection ledger](../research/dunbrack-prior-evidence-selection.json)
and [accuracy evidence](../research/dunbrack-prior-evidence.json).
New scripts `research_dunbrack_prior.py`, explicit amended adapter
`research_dunbrack_prior_weights.py`, auditor `audit_dunbrack_prior.py`;
tests `test_dunbrack_prior.py`, `test_dunbrack_prior_weights.py`,
`test_audit_dunbrack_prior.py`. No production/installed-library/fixture change.

Both already evaluated cohorts retained: 22 masks, 15 source-bound API queries
(eleven old/four new-input NCAC pools), seven inherited refusals/context blockers.
Queries use only masked-source phi/psi and canonical chi definitions, no PDB pose
import, coordinate generation, packing, H optimization or ref2015 calculation.
API type is single-chi Dunbrack library; active binary is
`rotamer/shapovalov/StpDwn_0-0-0/Dunbrack10.lib.bin`, SHA/size recorded, NOT
automatically the upstream recommended 5%-stepdown release. No bundled table or
per-case sample/weight/mean/sigma is included in tracked evidence/public artifacts.
Bounded private API responses remain research-only.

First frozen plan/run failed its normalized-probability requirement before any
complete selection ledger. API probability() raw values have totals
0.960705–1.041386 here and a value can exceed one. Original script/plan/failure
remain immutable. Separately frozen v2 adapter ranks finite nonnegative raw
weights with positive sum; no normalized probability, Gaussian offset or posterior
claim. Diagnostic normalized display is not used in selection. Unchanged common
geometry gates precede ranking, all three family wells retained, stable index tie.
No angle/source/geometry/score-weight tuning. API mechanism cause is unassigned.

All 15 phi+360/psi−360 controls pass; second independent process repeats entire
plan/environment/raw query/weight-semantics files byte-identically. Independent
auditor commits sanitized selections before native reporting. A cross-Python
diagnostic sum differed at 1e−16 and the auditor now checks those sums/displays
within 1e−14; original failed audit retained, no geometry/output-byte tolerance
relaxed. Focused suite 712 passes in 26.91 s (not full repository acceptance).

Prior-only results are negative as a rescue: on the four new-input pools every
choice equals the parent fixed-heavy choice, including bad THR11; prior/parent
median 1.222795 Å vs zero-angle 0.438949. Old-eleven prior/parent median both
0.226284, but 1R69 THR26 changes and WORSENS 0.337043→2.083064 Å. Aggregate
median equality hides it. All selected geometry passes, not accuracy qualification.

Private artifacts: `tmp/pdbfixer-research/dunbrack-prior-plan-01.json`, failed
`dunbrack-prior-run-01/`, amended `dunbrack-prior-plan-02.json`,
`dunbrack-prior-run-02/` and exact repeat `dunbrack-prior-run-03/`. Scientific
history remains separate; no commits/pushes or original-table redistribution.

**Resume:** separately freeze actual library-mean proposal generation, preserving
CB and all observations and comparing unchanged parents, if continuing this arm.
It has NOT executed here. Family-only ranking cannot be called a calibrated
context selector, and no public capability or production approval follows.
Keep original-table/license request to user mandatory before production, as above;
local PyRosetta research is already authorized. Final checks complete: both
sanitized publications replay byte-identically, all 13 old/new final verifiers,
55 fixture hashes and 1,156 local targets (including AGENTS.md) pass. Targeted
Ruff, 17-file mypy, generated CLI/GUI checks, agent maps and whitespace pass.
Final verifier/quality: `tmp/pdbfixer-research/validate-dunbrack-prior-final.py`,
`dunbrack-prior-docs-quality.md`.

## Earlier checkpoint: prospective inputs do not confirm score-selection benefit

Read [negative new-input result](../research/ncac-new-system.md),
[independent selection ledger](../research/ncac-new-system-evidence-selection.json)
and [accuracy evidence](../research/ncac-new-system-evidence.json). New scripts
`prospective_ncac_cohort.py`, `audit_prospective_ncac.py`; tests
`test_prospective_ncac_cohort.py`, `test_audit_prospective_ncac.py`.
All production consumers/dependencies/installed libraries and thresholds unchanged.
No commit/push. Standing instruction: continue without routine approval questions.

Enrolled 1PGB/1BPI/1CSP BEFORE download, first model/first ATOM chain and first two
VAL/THR identities. Whole masks only, same NCAC + 0/+120/−120° algorithm, unchanged
own baseline and ref2015 fixed-heavy scorer. No system replacement/repair/altloc
selection, native angle tuning or weight fitting. Twelve-mask ceiling becomes
eleven masks: 1BPI has one VAL. Original downloads now tracked source-only under
`tests/fixtures/protein_context/ncac_followup/`; manifest extends 52→55 inputs,
all old entries/digests preserved. Fixture README records source provenance.

Three 1BPI shared refusals (altloc/occupancy); four 1CSP NCAC protocol blockers
(missing non-target context sidechains trigger unexpected overlay). Stock also
refuses unrequested context additions on 1CSP. Retain all raw failures; no failure
counted as executed geometry success. Four 1PGB NCAC traces: 20 full snapshots,
four helper returns/eight fits; all twelve relative proposals pass geometry.
Stock passes two of four 1PGB common gates. Sources/CB remain fixed in admitted
selectable outputs. All partial failure files bound by independent retained hashes.

Thirty independently inspected scorer imports: 24 scored, six genuine import
blockers on 1CSP baseline. Rosetta adds sixteen GLU context heavy atoms at
A/3,21,36,66 and moves observed PRO58/CD. Contract rejects these; never restore
coordinates post hoc or exempt scorer-owned additions. All four 1CSP baseline
selections block; only four 1PGB cases have matched selected outputs.

Prospective result is NEGATIVE for the scoring gain: NCAC zero/selected medians
0.438949/1.222795 Å vs unchanged own-selected 1.172077 Å. Lone changed NCAC choice
THR11 +120° worsens 0.551633→2.119324 Å. VAL21 remains poor (2.360604 Å).
No accuracy superiority or three-system qualification. Native common passes six
of eight admitted cases: THR11 clash and VAL29 CB–CG2 length fail unchanged gates.
Denominators always 11 enrolled / 8 admitted / 4 selectable; no pruning of 1BPI,
1CSP, count shortfall, scorer movement or native-control caveats.

Private artifacts `tmp/pdbfixer-research/ncac-new-system-{enrollment-01.json,
download-01/,masks-01/,run-01/,score-run-01/}`. Geometry/selections audited and
published before native reference readback. Focused suite 696 passed in 27.66 s;
not full repository/platform acceptance. All changes remain uncommitted.

**Next boundary:** do not tune old angles/weights to manufacture a rescue.
Broader route needs separate context-completion admission and calibrated priors;
neither is established. Authentic Dunbrack raw-library/download/license lineage
still unavailable: official HTTPS page fetches and direct curl connection time
out; search snippets do not authorize extraction of bundled PyRosetta tables.
For local empirical-library research, use the installed PyRosetta route under
the current authorization above. Original upstream tables and documented rights
remain a mandatory request to the user before any production transition.
No claim that the library requires payment/restricted registration. Production
and public completion remain unchanged/disabled. Final readback/quality checks
are complete: all three new publication pairs repeat byte-identically; 12 final
verifiers, 55 fixture hashes and 1,126 local Markdown targets pass. Targeted Ruff,
17-file mypy, generated CLI/GUI checks, agent maps and whitespace pass. Final
verifier/quality: `tmp/pdbfixer-research/validate-ncac-new-system-final.py`,
`ncac-new-system-docs-quality.md`. All three documentation plans are complete.

## Earlier checkpoint: relative-rotamer selection improves four NCAC THR masks

Read [relative-rotamer results](../research/ncac-relative-rotamers.md),
[independent selection ledger](../research/ncac-relative-rotamers-evidence-selection.json)
and [accuracy evidence](../research/ncac-relative-rotamers-evidence.json).
Protocol `tmp/pdbfixer-research/ncac-rotamers-protocol-01.json` was frozen BEFORE
candidate generation. Scripts `research_ncac_rotamers.py`, `audit_ncac_rotamers.py`;
tests `test_ncac_rotamers.py`. Production, original fixtures and installed libraries
unchanged; no commit/push. Standing instruction remains autonomous continuation
until user access/resources or a substantial scope decision is genuinely needed.

Exactly 0/+120/−120° relative rigid rotation of BOTH generated distal branches
around saved CA–CB, preserving CB/backbone/observations. No minimization, empirical
library, native-tuned angles/weights or repacking. All eleven saved masks retained:
33 proposals, three clash refusals, 30 common-feasible scorer inputs. Runner plan
has no native reference/control and scores only the complete feasible pools with
unchanged ref2015 score-only imports, no H optimization. All eleven selections
commit before the independent auditor reads any native reference PDB.

Independent audit verifies source records, axial/radial geometry, signed relative
angles, rigid pair distances, common raw/rounded audits and all 30 physical/graph
import snapshots. Extras/virtuals/scorer context H are identical across alternatives
for all cases. Known atoms/CB stay fixed, selected CA S and THR CB R. Parent zero-
angle accuracy repeats exactly. No observed motion hidden by coordinate restore.

Four changed choices, all THR, four native accuracy improvements, zero worsened
on this retrospective set: 2CI2 THR22 +120 (1.945968→0.176830 Å), 1R69 THR26 −120
(1.848568→0.337043), THR27 +120 (1.970131→0.220035), 1CTF THR76 +120
(2.061262→0.331298). All-eleven zero/selected/previous-own median rounded exact
RMSD 0.322945/0.226284/0.293937 Å; seen-four 1.113282/0.198433/0.157531;
new-seven 0.322945/0.322945/0.383067. Seen-four remains worse than previous own.
Three rejected rotations (VAL32 both alternatives, THR55 −120) never scored.
Native-control failures stay unchanged. No independent/broad superiority claim.

Artifacts: `ncac-rotamers-score-plan-01/`, `ncac-rotamers-score-run-01/` under
`tmp/pdbfixer-research/`. Focused suite 689 passed in 26.46 s. Routine six-file
docs plan `ncac-rotamers-docs-plan.md`. Direct official Dunbrack HTTPS connection
also timed out (curl exit 28); no authentic asset downloaded, no PyRosetta table
extraction or license claim inferred. This is a network/provenance issue, not
proof the resource requires a paid/restricted license.

**Resume:** next accuracy gate is a separately enrolled, outcome-independent
VAL/THR cohort under this frozen algorithm, not angle densification on old errors.
Preserve four-seen/seven-new-variant history and whole-mask-only restriction.
Keep empirical-library provenance, calibrated chemistry, ILE/multi-χ and full
legacy comparisons distinct; public replacement remains disabled. Finish any
pending final readback/quality checks before enrolling that new cohort.

## Earlier checkpoint: all eleven NCAC whole masks pass geometry, accuracy mixed

Read [eleven-mask results](../research/ncac-whole-mask.md),
[geometry ledger](../research/ncac-whole-mask-evidence-geometry.json) and
[accuracy evidence](../research/ncac-whole-mask-evidence.json). Frozen plan/run:
`tmp/pdbfixer-research/ncac-expansion-plan-01.json`, `ncac-expansion-run-01/`.
Scripts: `research_ncac_expansion.py`, `audit_ncac_expansion.py`; tests:
`test_ncac_expansion.py`. Production, installed libraries, thresholds, original
fixtures and previous evidence remain unchanged. All work remains uncommitted.

All eleven original whole-sidechain VAL/THR masks were frozen before new runs,
with four already seen NCAC cases and seven new only to this variant. Twenty-two
stock/NCAC normal-relaxation arms produced 34 helper states, 170 stage snapshots,
67 target minimizer iterates, 22 checked NCAC fits and zero dynamics calls.
All 23 stock seed states reproduce the original untraced baseline; all four seen
NCAC complete traces/payloads repeat exactly. NCAC returns on seed 1 throughout.
Observed O remains fixed, merely excluded from fitting; all source observations
are preserved. Every NCAC construction/return has CA S and all THR endpoints CB R.

Independent auditor publishes geometry before opening native reference PDBs;
accuracy evidence binds its digest and stores a relative companion basename.
Fresh destination pairs are required. Existing knowledge of original cohort
native outcomes remains explicit; this is not blinded/new-family enrollment.
Common passes stock/NCAC/own-scored/native: 4/11/11/9. Native VAL32 clash and
THR76 angle failures remain visible with unchanged cutoffs. Own-scored selection
replays the previous fixed-heavy score choice, without native-based reselection.

Matched NCAC/own-scored rounded exact RMSD medians Å: all eleven
0.322945/0.293937; seen four 1.113282/0.157531; new seven 0.322945/0.383067.
Stock/NCAC mutual-pass count four, medians 0.215798/0.218371; all-raw stock
1.252770 is diagnostic, not an accepted-output accuracy comparison. No mutual
stock/NCAC passes among seen four: matched medians null, never zero.
Geometry rescue does not settle conformation; several NCAC THRs remain ~2 Å
from native. No broad superiority, full-legacy qualification or public remedy.

Expanded focused suite: 678 passed in 27.11 s. Six-file publication tracking:
`tmp/pdbfixer-research/ncac-expansion-docs-plan.md`. Original four-case and older
evidence stay separate. No public artifact, GUI integration or production change.

**Resume:** freeze a separate retrospective relative-rotamer test on these eleven
saved NCAC proposals: preserve CB and all observations, rotate generated distal
branches together by 0/+120/−120° around CA–CB, audit unchanged gates, score the
complete feasible pool with fixed-heavy Rosetta, commit selection before native
reporting. No native tuning of angles, thresholds or score weights; retain every
rejection and empty pool. This is not yet enrolled/executed. Verify official
Dunbrack download/license lineage separately; official page fetches timed out,
so search snippets alone have not established downloadable asset provenance.

## Earlier checkpoint: N–CA–C proper-fit control passes four mechanism cases

Read [used-anchor controls and bounded rescue result](../research/template-anchor-controls.md)
and [independent used-anchor/stage/common-gate evidence](../research/template-anchor-controls-evidence.json).
New scripts: `research_template_anchors.py`, `audit_template_anchors.py`;
tests: `test_template_anchors.py`. Old tracer/auditor files, installed scientific
libraries, source observations, thresholds and production code remain unchanged.

Before new generation, froze ALL four previous case roles and a 4x2x2 factorial:
proper fit with all supplied N/CA/C/O vs exactly N/CA/C, normal vs no relaxation.
The two confirmed-D 2CI2 whole-sidechain cases, stock-unassigned 1R69 THR27 and
identity-selected 2CI2 THR22 control all remain; no successful-rescue filtering.
O remains an observed fixed atom in the model/force field, merely excluded from
NCAC fit orientation/centering. Exactly four supplied anchors are source/resource-
bound; unexpected sets/order refuse. This is NOT support for retained distal atoms,
arbitrary partial residues, loops, ligands or PTMs.

The unchanged legacy tracer records FOUR supplied anchors and rotation class
`proper`. Separate `anchor-policy.json` records every actual USED name/index,
coordinate set, centroid and transform. Do not read the old `anchors=4` field as
used-anchor count. A process-local research fit hook supplies the subset and is
restored; parent source files and installed PDBFixer are not edited.

Executed 16 arms, 40 helper states, 200 full stage snapshots, 43 target reporter
iterates and 80 independently checked used fits. ALL eight four-anchor controls
reproduce their complete parent proper traces/payloads exactly (32 seed states).
Zero escape-dynamics calls. Each NCAC arm returns on seed 1, with all observations
preserved. Auditor checks transport against original sources/resources, proper
orthogonality/determinant, exact used centroids, optimum residual singular-value
bound, complete-template volume law and actual fixed-backbone hybrid construction.
Every returned seed receives the unchanged measured-chi common raw/rounded/global/
target geometry gate. No PDBFixer instance, generation or optimization in audit.

NCAC normal initial/final CA S: 4/4; final common-gate passes: 4/4 versus four-
anchor normal 0/4. THR22 was L but not a geometry-passing control: all-anchor
CA–CB length and N–CA–CB angle fail. Helper return is not independent gate pass.
NCAC initial→final signed volumes Å³: VAL28 +2.423129→+2.423129; THR55
+2.459139→+2.459139; THR27 +2.450350→+2.306998; THR22 +2.365314→+2.365314.
All THR NCAC endpoints are CA S/CB R. NCAC no-relaxation arms also pass and
preserve construction exactly, but skipping relaxation is NOT a user remedy.
The isolated anchor-set change supports this mechanism on these cases, not proof
that observed O is wrong, universal L safety or native conformation accuracy.

Artifacts: `tmp/pdbfixer-research/anchor-plan-01.json`, `anchor-run-01/`,
`anchor-audit-01.json`. Diagnostic arrays/transports/logs are private; no final
PDB/dat, public bundle or GUI artifact produced. No native accuracy, scoring
change or new independent-family qualification. 673 focused tests pass in
25.95 s and targeted Ruff passes.
Published used-anchor/common-gate audit repeats byte-identically. Package mypy
(17 files), generated CLI/GUI checks, agent maps and whitespace pass. All previous
evidence verifiers and 52 source-fixture checksums remain valid; all 1,107 local
Markdown targets across docs plus fixture README exist. Final verifier/quality:
`tmp/pdbfixer-research/validate-anchor-final.py`, `anchor-docs-quality.md`.

**Resume here:** separately freeze ALL eleven whole-sidechain masks from the
existing three-system cohort, not only four rescues. Distinguish four already
tested cases from seven not yet run under NCAC. Assess native accuracy only after
proposals, retain every raw stereo/geometry failure and native-control deviation;
do not relabel this a new blinded/family-independent cohort. General partial-
completion/public routes remain disabled; no installed-library or production patch.

## Earlier checkpoint: actual construction already D; proper-fit hybrid also D

Read [factorial stage tracing and hybrid counterexample](../research/template-initialization-causality.md)
and [independent stage/template evidence](../research/template-initialization-causality-evidence.json).
New scripts: `research_template_causality.py`, `audit_template_causality.py`;
tests: `test_template_causality.py`. Retrospective case mechanisms, not an accuracy
cohort or production remedy. Installed libraries, production code, original
comparison evidence and source observations remain unchanged.

Before new traces, froze four cases: confirmed-D whole 2CI2 I/VAL28 and I/THR55;
small-volume/CA-CIP-unassigned whole 1R69 A/THR27 (NOT confirmed D); and the first
whole case in original enrollment order, 2CI2 I/THR22, not success-selected.
Full 4x2x2 factorial: stock vs determinant-corrected proper overlay, normal vs
no relaxation. No-relaxation bypasses BOTH minimization and escape steps, is an
ablation only and is not a geometry remedy. Same observed-anchor set, same helper
seeds 1..5, no extra construction preview, no native coordinate use.

Executed 16 arms / 60 returned helper states, retained 300 full stage snapshots
and 75 target-coordinate minimizer reporter iterates. ALL four stock-normal arms
reproduce all sixteen prior untraced seed states exactly (identities/elements/raw
coordinates); observational hooks do not change those outputs. Hooks restore per
arm and installed source/resource hashes stay unchanged. Independently audited
all boundary observations and target-backbone reporter points. Zero dynamics calls;
returned states repeat exactly within every arm on pinned Reference. These seeds
are not independent conformational samples or evidence about other platforms.

Primary cases are already CA R at construction: stock VAL28 volume -0.38863 Å³
before minimize, -0.62739 after; THR55 -0.62022 before, -0.55516 after. With
relaxation disabled their constructed coordinates remain unchanged and wrong-side.
Thus optimization is NOT necessary to create these two D placements. CB is missing
in the source: no observed sidechain atom was reflected/moved.

Proper correction does NOT rescue them. Full transformed canonical templates stay
L (VAL volume +2.407856, THR +2.457112 Å³, CA S), but mixing their new sidechain
with the fixed observed backbone exactly reproduces D-side hybrids: proper VAL28
construction -1.555334, THR55 -1.275285, THR27 -1.083445, all CA R. This is checked
from original resource atoms and saved transforms, with independent determinant
and scalar-triple calculations. det(R)=+1 protects a fully transformed template,
not a hybrid whose original observed backbone was not transformed. No universal
anchor mechanism or successful determinant-only patch is established.

THR27 STOCK volume changes +0.06179→-0.08247 during minimization with CA CIP
unassigned at both boundaries: a small-volume sign change, NOT confirmed L→D.
THR22 control stays CA S/CB R, positive volume +3.09377, across all four arms.
Proper correction of THR27 instead starts CA R and is not a remedy. All raw
diagnostics remain private; no final PDB/dat, public bundle or GUI output is made.

Artifacts: `tmp/pdbfixer-research/causality-plan-01.json`, `causality-run-01/`,
`causality-audit-01.json`. Independent auditor uses saved stages/resources, no
PDBFixer instance or coordinate generation, minimization or dynamics. 664 focused
tests pass in 25.80 s and targeted Ruff passes. Endpoints/raw stages are not full
legacy outputs; downstream production chirality safety nets remain unchanged.
Published saved-stage audit reproduces byte-identically. Package mypy (17 files),
generated CLI/GUI checks, agent maps and whitespace pass. All previous evidence
verifiers and 52 fixture checksums remain valid; all 1,103 local Markdown targets
across docs plus fixture README exist. Final verifier/quality:
`tmp/pdbfixer-research/validate-causality-final.py`, `causality-docs-quality.md`.

**Resume here:** separately declare N–CA–C-only proper-fit control on ALL four
cases to test the fourth-anchor hypothesis, without selecting successful rescues.
Keep observations and independent hybrid/stereo/geometry audits fixed. Do not
ship an untested anchor change or determinant-only patch. Empirical priors,
calibrated context chemistry and broad/full-legacy qualification remain open.

## Earlier checkpoint: pre-download new-system accuracy/stereo controls

Read [enrollment, ranking gains and raw-template stereo failures](../research/rosetta-independent-controls.md),
[independent cohort evidence](../research/rosetta-independent-controls-evidence.json) and
[unchanged score readback](../research/rosetta-independent-controls-readback.json).
New scripts: `prospective_independent_cohort.py`, `audit_independent_cohort.py`;
tests: `test_prospective_independent_cohort.py`. No old generator/scorer/gate,
runtime dependency, public request/bundle or production consumer changed.

Before input download/inspection, enrolled 2CI2/CI2, 1R69/434 repressor amino-terminal
domain and 1CTF/L7-L12 carboxy-terminal domain. Rule: first coordinate MODEL, first
ATOM chain, up to first two VAL/two THR in file order, all four masks. Preserve
H/altloc/occupancy and relevant original links; no repair, replacement or outcome
filter. 2CI2 selects chain I, the other two chain A. 1CTF has only ONE observed THR:
retain its four-mask shortfall. Exactly 44 cases, not 48; all admit and both arms run.
Raw original RCSB inputs are under `tests/fixtures/protein_context/independent_followup/`
with README provenance and MANIFEST hashes. New to this study, not randomized,
blinded, temporally held out or independent of Rosetta training/templates.

Use original three-chi/anchored own proposals and the unchanged shared seeded
PDBFixer helper on Reference (PDBFixer 1.12.0/OpenMM 8.5.2). Preflight/resources usable;
no installation/repair. Same measured-chi common gate, fixed-heavy ref2015 import,
fresh unoptimized H, score weights/seed and complete feasible-pool selection.
All 176 poses score and pass independent exact physical/source/graph/context/score
readback; no context blocker. Source observations survive. There are 56 retained
template seed states. All reconstruction/environment/protocol blockers are zero.
An initial wrong `score` CLI invocation failed in argparse before science; retained
log and corrected `run` command, not a protocol amendment or a scientific failure.

Own scored/clearance coverage stays 41/44, template common-gate coverage 16/44;
native common-gate passes 37/44. Per-system own/template/native: 2CI2 15/4/13 of 16,
1R69 16/8/16 of 16, 1CTF 10/4/8 of 12. Three empty own pools retained: 2CI2 VAL32
CG1 and 1CTF THR76 OG1/CG2. Native failures are three correlated VAL32 CG1–C
1–4 contacts at 2.758456 Å (cutoff 2.8) and four THR76 CA–CB–CG2 angle failures
at 125.118487° (illustrative tetrahedral target ±12°). No native gate retuning;
filter-specific coverage is NOT accuracy superiority. On 16 mutually passing masks,
rounded medians own/scored 0.197650 Å vs template 0.190375 Å: no own superiority.

Six score changes improve native RMSD, representing THREE correlated residue events:
2CI2 VAL28 pair 2.435→0.124 / whole 2.093→0.168 Å; 1R69 THR26 pair 2.536→0.269
/ whole 2.056→0.383; 1CTF VAL66 pair 2.609→0.277 / whole 2.182→0.352. All 22 paired
matched medians 0.524→0.249 Å; 19 anchored medians unchanged 0.128 Å. Generated-only
VAL symmetry is reporting-only and preserves these gains; no truth selection.

Important chirality evidence after seeded retry exhaustion, NOT full legacy output:
whole 2CI2 VAL28 raw Cα volume -0.627386 Å³ + rounded CA R, and whole 2CI2 THR55
volume -0.555162 + CA R/CB S are two confirmed D-Cα proposals. Whole 1R69 THR27
volume -0.082474, CA unassigned/CB S is a small-volume collapse, NOT confirmed D.
All three retain five seeds/logs and fail the independent gate. All 41 selected own
outputs pass global/target stereo. No reflection or silent observation restoration;
production legacy downstream chirality defenses remain unchanged. Endpoints alone
do not establish whether initialization or optimization caused these outcomes.

Artifacts below `tmp/pdbfixer-research/`: `independent-enrollment-01.json`,
`independent-download-01.json`, `independent-masks-01/`,
`independent-comparison-01/`, `independent-score-plan-01/`,
`independent-score-run-01/` and saved readbacks. Auditor replays source/masks,
saved geometry/clearance, template observations/additions/snapshots and common
audits without rebuilding coordinates. 658 focused tests pass in 26.67 s and
targeted Ruff passes. Scoring readback retains historical benchmark flags; new
prospective provenance is separate and never rewrites its older scope labels.
Both publications reproduce byte-identically. Package mypy (17 files), generated
CLI/GUI checks, agent maps and whitespace pass; all old evidence verifiers remain
valid. All 49 previous fixtures plus three unchanged new sources verify (52 total);
1,098 local Markdown targets across docs plus fixture README exist. Final checks:
`tmp/pdbfixer-research/validate-independent-final.py`, `independent-docs-quality.md`.

**Resume here:** separately specify causal initial-overlay versus optimizer-stage
controls on the two confirmed-D whole-sidechain cases, with all observations fixed
and no installed-library or production edits. Keep the third unassigned case
distinct. More empirical priors, calibrated context geometry and full-legacy/new-
family qualification remain open. Do not enable public completion or approximate
coupled packing from these bounded results.

## Earlier checkpoint: exhaustive all-reference replay, feasible interactions but no selection gain

Read [all-reference/feasible-rectangle analysis and stopping result](../research/rosetta-reference-controls.md),
[complete analysis](../research/rosetta-reference-controls-evidence.json) and
[independent raw-cell readback](../research/rosetta-reference-controls-readback.json).
New stdlib-only scripts: `analyze_rosetta_references.py` and
`audit_rosetta_references.py`; tests: `test_rosetta_references.py`.
Retrospective descriptive analysis only: no sampling, Rosetta energy evaluation,
geometry change, native-accuracy selection, selected PDB/dat change or public API.
Every prior sampling/scoring implementation and evidence file remains unchanged.

A separate descriptive plan binds the already seen dense audit, input/raw/coarse
scores and implementations before this new analysis. Retain ALL 141 reference
choices over the same twelve matrices and ALL 381 distinct rectangles. Of those,
255 have four joint-feasible corners; every one has nonadditivity above the fixed
1e-9 reporting tolerance. Seven two-axis matrices each contain feasible rectangles;
five single-axis matrices remain uninformative, not evidence of zero interactions.
References/rectangles are correlated diagnostics, not independent biological cases.

Per-reference surrogate is E(i,b)+E(a,j)-E(a,b) on the same global feasible set.
Reference classes: 12 rejected anchors produce eight different choices; 25 feasible
anchors with rejected axis terms produce two; **104 references with ALL used terms
feasible produce ZERO**. The 104 include 86 genuine two-axis references and 18
single-axis references. A feasible anchor alone is not enough: inspect all terms.

Each of 129 feasible anchors also defines a supported subset with every four-corner
term feasible. Compare full and surrogate minima on THAT SAME subset: zero choice
differences/positive score gaps. The global winner remains in 121 supported pools;
eight exclude it and must NOT be reported as unchanged global coverage. Rejected
anchors have no supported subset.

THR43/THR51 feasible-only max rectangle deltas are 3.008091/0.401633 Rosetta units,
not the earlier rejected-corner first-reference 48.693/5.441. Other feasible maxima
are 0.0064–0.4273 units. Interaction effects are real under the frozen full-score
protocol, but all observed choice discrepancies rely on rejected terms. No general
independent/coupled-packer equivalence or practical packing advantage is established.
Earlier ZERO changed coarse-to-dense winners and all source-contact caveats stand.

Independent auditor uses ORIGINAL raw scored cells and sets of four eligible cells,
not the analyzer's audit-derived inputs/axis-pair enumeration. It replays every
reference, subset, minimum, rectangle/class and aggregate with no scientific-stack
imports. Both output copies reproduce byte-identically. Local records:
`tmp/pdbfixer-research/reference-plan-01.json`, `reference-analysis-01.json` and
`reference-audit-01.json`. 649 focused tests pass in 26.39 s; targeted Ruff,
17-file package mypy, generated CLI/GUI checks, agent maps and whitespace pass.
Both publications reproduce byte-identically. Dense/contact/joint final verifiers,
prior follow-up/comparison evidence and all 49 fixtures remain valid. All 1,091
local Markdown targets across docs plus fixture README exist. Final checks:
`tmp/pdbfixer-research/validate-reference-final.py`, `reference-docs-quality.md`.

**Resume here:** this bounded coupling arm has a stopping result. Do NOT keep
densifying these same four pairs to manufacture a positive claim or implement an
approximate graph/public routing from them. Further accuracy qualification needs
a separately specified independent-family/reference-free comparison with all
refusals retained. Context chemistry and empirical candidate priors are still
separate open questions; saved-matrix analysis is not new prospective qualification.

## Earlier checkpoint: separate nine-chi grid, no changed selected conformations

Read [dense-grid policy, exact repeats and conditional-reference caveat](../research/rosetta-dense-controls.md)
and [independent dense evidence](../research/rosetta-dense-controls-evidence.json).
New scripts: `research_rosetta_dense.py`, `rosetta_dense_diagnostics.py` and
`audit_rosetta_dense.py`; tests: `test_rosetta_dense.py`. Original generators,
validation defaults, joint runner/auditor, score policy, dependencies and public
consumers remain unchanged. This is a distinct private nine-chi protocol, not
old-grid generator admission or production materialization.

Before new sampling/scoring, enrolled ALL four frozen contact pairs and ALL
twelve mask cases. Paired grid ordered -60,60,180,-80,-40,40,80,160,-160 degrees
(original three plus +/-20-degree offsets); nine proposals/paired target, at
most 81 joint combinations/case. Anchored mode retains its one exact fixed-partner
closure. Original three proposals/audits and all observations must reproduce.
New grid requires explicit membership then existing measured-chi comparator
geometry with unchanged raw/rounded stereo, clash, bond and angle thresholds.
Production still rejects the additional angles. Standard DatRecord sidecars
explicitly label private dense ownership, never observed/scorer-created atoms.

Executed 141 complete source-feasible combinations; 12 raw joint rejections,
zero extra rounded rejections, 129 private eligible PDB/dat outputs. All 153 poses
(141 own +12 native) score with 306 exact target-view imports and independent
stereo/graph/source readback. Twelve private selected pairs match eligible cells.
All 38 old coarse/native PDBs, full imported snapshots and scores repeat exactly.

Seven matrices now have two varying axes and nonzero enrolled full-score
nonadditivity; five remain single-axis/uninformative. THR43/THR51 distal/whole
max deltas are 48.693/5.441 Rosetta units, but BOTH first-source-feasible reference
poses fail joint geometry. Other two-axis max deltas are 0.0055–0.4202 units.
This diagnoses richer interaction-sensitive tables, not calibrated atom-pair
energies or a generally qualified coupled packer.

Principal negative result: **ZERO coarse-to-dense score-selection changes**.
All twelve full-score winners already belong to the coarse grid; matched median
pooled exact-name RMSD stays 0.202370 Å. No added selected accuracy/coverage.
An enrolled conditional factorized surrogate E(i,j0)+E(i0,j)-E(i0,j0), restricted
to the same joint-feasible set, disagrees in two THR43/THR51 masks: full score
reductions 7.103/4.632 units and pooled RMSD 1.484→0.140 / 0.315→0.209 Å. BOTH use
the geometry-rejected reference: do NOT call this robust practical coupled-packing
benefit. Neither independent clearance nor independent conditional winners fail
joint geometry; no winner rescue. Truth never selects a candidate/reference.

Artifacts: `tmp/pdbfixer-research/dense-enrollment-01/`, `dense-plan-01/`,
`dense-score-01/`, audit revisions `dense-audit-{01,02}.json`; selected private
`selected-candidate.{pdb,dat}` files per case. The latest independent auditor
replays saved geometry/enumeration, not coordinate generation, and checks every
raw/rounded/source/graph/score/context/stereo/DatRecord binding and exact coarse
repeat. Native geometry 9/12 vs source-contact refusals stays a separate inherited
control, never native grid admission. Public diffusion/completion bundles remain
unmodified. 626 focused tests pass in 26.84 s; targeted Ruff, 17-file package mypy,
generated CLI/GUI checks, agent maps and whitespace pass. Complete saved-candidate
audit is byte-identical; old contact/joint/follow-up/comparison evidence and all
49 fixtures still validate. All 1,085 local Markdown targets across docs plus
fixture README exist. Final verifier/quality record:
`tmp/pdbfixer-research/validate-dense-final.py`, `dense-docs-quality.md`.

**Resume here:** specify a descriptive reference-robust analysis of saved matrices,
retaining ALL reference choices and whether fully feasible rectangles exist.
Do not resample, choose a favorable reference or relabel this retrospective
analysis as new independent prospective qualification. Coarse-to-dense selection
benefit is absent. Empirical sampling, context chemistry, independent-family/full
legacy qualification and production integration remain open. No approximate graph.

## Earlier checkpoint: input-backbone-selected close contacts

Read [close-contact enrollment, new–new clashes and limits](../research/rosetta-contact-controls.md)
and [portable independent evidence](../research/rosetta-contact-controls-evidence.json).
New scripts: `research_rosetta_contacts.py` (enrollment only),
`audit_rosetta_contacts.py` (input/mask and supplemental readback); tests:
`test_rosetta_contacts.py`. Original joint runner/auditor, source-package kernels,
scientific thresholds, score policy, production defaults and dependencies unchanged.

Before new reconstruction/scoring, select same-chain VAL/THR pairs by observed
CA distance <=8 Å, minimum three positions apart in residue file order, ascending
distance/file-order ties, up to four pairs per structure. Hidden sidechain
coordinates and reconstruction outcomes do not enter selection. All eligible
rankings/masks/code hashes are frozen. 1CRN has ZERO eligible pairs: retain that
result, never substitute/expand the cutoff. 1AKI has four pairs: THR89/VAL92,
VAL2/THR40, THR43/THR51 and THR51/THR69 (CA distances 5.169/6.355/6.549/7.972 Å).
Three original mask classes yield 12 admitted cases; no case replacement/filtering.

Exact runner produces 26 source-feasible combinations; three fail raw joint
geometry, leaving 23 canonical eligible outputs/DatRecord sidecars, no extra
rounded failures. All 38 complete poses (26 own +12 native) score and pass 76
exact target-view imports and stereo checks. Twelve private selected PDB/dat
pairs match eligible cells. No rejected combination acquires a selected output
or materialized sidecar; observation bytes and physical-heavy coordinates survive.

New result: THR43/THR51 source-only-feasible proposals collide after combination:
distal [0,0] CG2–CG2 at 2.455 Å, distal [2,0] OG1–CG2 at 2.550 Å, whole [0,0]
CG2–CG2 at 2.689 Å. Both target audits independently identify conflicts. This is
real joint-gate mechanism evidence, BUT the original independent clearance winner
already avoids every collision (zero independent-winner joint failures). All 12
matrices still have a one-option axis, so energetic coupling remains UNMEASURED.
Forty-two source-rejection finding occurrences across correlated case-targets are
all generated clashes, not 42 independent residues/failures.

Score changes three correlated raw-joint-clearance choices: THR89/VAL92 distal
1.782→0.238 Å, same pair whole 1.464→0.196 Å and VAL2/THR40 distal
1.738→0.344 Å pooled exact-name RMSD; 12-case medians 0.233→0.202 Å. This is
candidate ranking within old gates, not new coverage or a coupled-packing benefit.

Native geometry needs measured chi, not staggered-grid admission. Initial retained
`contact-supplement-01.json` sent metadata-free native candidates to the generator
audit: 0/12 is protocol incompatibility, NOT twelve geometric failures. Corrected
supplement uses existing comparator measured-chi geometry only, including the
other native target; no threshold change. Nine of twelve native geometry controls
pass. Three correlated THR89 masks retain the deposited OG1–ASN87 OD1 2.399187 Å
contact, clearance -0.640813 Å vs unchanged -0.6 cutoff. This is source geometry,
not reconstruction damage or proof of chemically impossible contact. Generator
materialization stays 0/12; native never becomes eligible. All 12 native poses
separately pass import/stereo checks. Keep those denominators distinct.

Artifacts: `tmp/pdbfixer-research/contact-enrollment-01/`, `contact-plan-01/`,
`contact-score-01/`, `contact-audit-01.json`, initial `contact-supplement-01.json`
and corrected `contact-supplement-02.json`. Portable supplemental evidence embeds
the original joint audit; final readbacks must use saved proposals, not resampling.
615 focused tests passed in 26.53 s; targeted Ruff, 17-file package mypy,
generated CLI/GUI checks, agent maps and whitespace pass. Complete joint and
supplemental saved-candidate audits reproduce byte-identically without generation.
Old joint/follow-up/comparison evidence and all 49 fixtures still validate; all
1,081 local Markdown targets across docs plus fixture README exist. Final checks:
`tmp/pdbfixer-research/validate-contact-final.py`, `contact-docs-quality.md`.

**Resume here:** separately predeclare a bounded richer-sampling study on ALL four
frozen pairs/masks, preserving empty/one-option pools, observations and geometry
gates. Diagnose coarse staggered sampling versus context-steric limitations before
any approximate graph. Do not filter on candidate success or tune a gate to native
failure. More context chemistry, independent-family/full legacy qualification and
production/public integration remain open. No new source fixture was needed.

## Earlier checkpoint: exact two-target controls, negative coupling evidence

Read [exact search, physical/sidecar checks and principal negative result](../research/rosetta-joint-controls.md)
and [independent evidence](../research/rosetta-joint-controls-evidence.json).
New research scripts: `research_rosetta_joint.py` and `audit_rosetta_joint.py`;
tests: `test_rosetta_joint.py`. No old source-package kernel, scientific threshold,
single-target scorer, runtime dependency or production consumer changed.

Twelve cases were enrolled before new joint reconstruction/scoring: first two VAL
and first two THR identities in file order of already-studied 1CRN/1AKI, split
into consecutive pairs; distal/distal, whole/whole, anchored-first/distal-second.
Pairs: 1CRN THR1/THR2 and VAL8/VAL15; 1AKI VAL2/VAL29 and THR40/THR43. Previously
known systems/constituents: this is NOT independent-family/prospective accuracy
qualification. Keep all cases, including known difficult anchored VAL8.

Unchanged mixed planning/generation enumerates the complete Cartesian product of
source-feasible pools, at most nine combinations/case. Exactly 34 combinations
were produced; all pass raw/rounded joint guards, canonical materialization and
independent target CIP. Twelve native controls remain ineligible. All 46 complete
poses score, with 92 fresh target-view imports: each target graph independently
checked in the same full physical-heavy/H/virtual/type/external-link pose, equal
finite scores required, no averaging. Fixed context is identical across each case.
Original records and exact heavy coordinates survive. All 34 DatRecord sidecars
contain exactly generated heavy atoms, not scorer H/OXT/virtuals. Eleven private
selected PDB/dat pairs match winning cells; the anchored VAL8 case has no source-
feasible first-target candidate and writes no selected artifact. No gate rescue.

Main negative result: ten candidate matrices have a single-option axis and one is
empty. Only 1AKI VAL2/VAL29 distal/distal is a genuine 3x3 matrix; observed CA
distance 13.128 Å, enrolled finite-difference max nonadditivity 1.7053e-13. That
score table is effectively additive on this grid, NOT evidence of generally absent
physical interactions. None of the independent per-target clearance choices fails
joint geometry, so new–new clash advantages are also not challenged scientifically.

Score changes two raw-joint-clearance choices: pooled rounded RMSD for 1CRN
THR1/THR2 distal/distal 1.739→0.183 Å and 1AKI VAL2/VAL29 1.707→0.114 Å; same
11-case medians 0.311→0.211 Å. The first winner already equals independent
per-target clearance; the second only repeats known VAL2 improvement, leaving
VAL29 unchanged. Do NOT call these a demonstrated coupled-packing improvement.

Local artifacts: `tmp/pdbfixer-research/joint-enrollment-01/`, `joint-plan-01/`,
`joint-score-01/` (per-case private `selected-candidate.pdb` / `.dat`) and initial
`joint-audit-01.json`. Final independent reader reuses saved candidates, not any
coordinate generator; it replays source/individual/joint/rounded guards, complete
enumeration, original clearance winners, both graph views, score terms/context,
global Cα/THR stereo, target CIP and canonical sidecars. Portable enrollment is a
protocol excerpt with original enrollment hash, not a rewritten frozen manifest.

Final verification: 611 focused tests in 24.97 s; targeted Ruff, seventeen-file
package mypy, generated CLI/GUI checks, agent maps and whitespace pass. Complete
saved-candidate/34-sidecar/92-view audit replays byte-identically. Old follow-up,
original/mixed evidence and all 49 fixtures remain valid. All 1,077 local Markdown
targets across docs plus fixture README exist. Final verifier/quality:
`tmp/pdbfixer-research/validate-joint-final.py` and
`tmp/pdbfixer-research/joint-docs-quality.md`.
**Resume here:** separately predeclare proximity-selected pairs/masks using INPUT
geometry, not reconstruction success, and retain collapsed/empty pools. This weak
coupling cohort does not justify implementing/qualifying interaction-graph
approximations. More feasible sampling, calibrated chemistry/protonation and all
production/public integration remain unfinished. No new source fixture was needed.

## Earlier checkpoint: predeclared scoring follow-up with explicit protocol amendment

Read [enrollment, corrected serialization and bounded results](../research/rosetta-scoring-followup.md),
[follow-up evidence](../research/rosetta-scoring-followup-evidence.json),
[independent score readback](../research/rosetta-scoring-followup-evidence-readback.json),
[reconstruction record](../research/rosetta-scoring-followup-reconstruction.json) and
[native audit](../research/rosetta-scoring-followup-native-audit.json).
New scripts: `prospective_branched_cohort.py` and `prospective_scoring_followup.py`;
new tests: `test_prospective_branched_cohort.py`. Old generators, scoring kernels,
selection/eligibility rules, scientific gates and production consumers are unchanged.

Before input download/inspection or new proposals/scores, enrolled 1CRN/crambin,
1AKI/lysozyme and 1BTA/barstar: first two VAL/two THR in first-model chain-A file
order, four masks each, no outcome filtering or replacement structures. Preserve
all selected ATOM bytes, occupancy/altloc/H and relevant link annotations. Raw
deposited sources are now tracked under
`tests/fixtures/protein_context/scoring_followup/`, with provenance and hashes.

Important protocol limitation: initial projection put CONECT before coordinates,
causing all 32 admitted template calls to fail in OpenMM `_current_model.connects`
before coordinate rebuilding. This was OUR serialization bug, not PDBFixer geometry
failure. A separate post-reconstruction enrollment amendment moves only CONECT
records after coordinates. All 48 identities, all non-CONECT bytes and CONECT
record multisets across 54 artifacts are identical; all 32 own results reproduce
exact additions/audits/selections. Original source/request hashes change as expected.
New-cohort scoring had not started; the chosen scoring rule/code did not change.
Do not call this a pristine unamended frozen protocol or randomized/blind benchmark.
Both enrollments, old projector bytes and failed/corrected runs remain retained.

Final outcomes: 48 declared, 16 source refusals (all barstar), 32 admitted. Barstar's
714 input H coordinates remain intact; no strip-to-pass. Own/template common passes 31/14,
zero template execution blockers after correction; native controls pass 32/32.
Crambin contributes own/template 15/9 and lysozyme 16/5. Own singleton 1CRN VAL8 CG1
still refuses its CG1–C test; template and native pass. No new gate or rescue.

All 128 inputs (64 own attempts, 32 template, 32 native) score with exact physical
heavy preservation, source records, checked external/target topology, and identical
within-case extra OXT/virtual/non-target-H context. Only eight cases have multiple
feasible own candidates (23 single-option, one empty, one two-option, seven
three-option pools). Same 31 selections; one changes: 1AKI VAL2 distal pair,
rounded exact-name RMSD 2.414→0.147 Å. Both-generated VAL swap sensitivity gives
1.771→0.147 Å, so not just methyl relabelling; this is reporting only.
Paired n=16 rounded median 0.247→0.242 Å; singleton n=15 remains 0.115 Å. This is
one useful new mask, not broad superiority, new coverage or family qualification.
The original template-mutual-pass raw paired subset still favors template
(n=6 own/template 0.248/0.187 Å); keep these denominators separate.

Local evidence under `tmp/pdbfixer-research/`: `prospective-enrollment-{01,02}.json`,
`prospective-inputs-01/`, `prospective-masks-{01,02}/`,
`prospective-reconstruction-{01,02}/`, `prospective-score-plan-01/` and
`prospective-score-run-01/`. Correction proofs:
`prospective-conect-amendment-audit.json` and
`prospective-own-replay-after-conect-fix.json`. The final score plan archives the
unchanged legacy builder manifest separately and adds truthful enrollment provenance;
generic readback's broad-benchmark qualification flag stays false.

Final verification: 591 focused tests in 24.61 s; targeted Ruff, seventeen-file
package mypy, generated CLI/GUI checks, agent maps and whitespace pass. All 49
fixture entries verify. The 128-input independent readback and 32 new native
audits replay byte-identically; the outer wrapper differs only in the output-path
reference. Amendment byte/own-result proofs recheck successfully. Prior scoring
and original/mixed evidence remain valid. All 1,072 local Markdown targets in docs
plus fixture README exist (code excluded). Quality/verifier:
`tmp/pdbfixer-research/prospective-docs-quality.md` and
`tmp/pdbfixer-research/validate-prospective-final.py`.
**Resume here:** separately enroll exact joint atom-fixed controls before any
interaction-graph approximation. More systems, feasible sampling, calibrated
sterics/protonation and production integration remain independent unfinished gates.

## Earlier checkpoint: atom-fixed scoring and descriptive topology diagnostics

Read [protocol, results and limits](../research/rosetta-fixed-scoring.md) and
[portable independent evidence](../research/rosetta-fixed-scoring-evidence.json).
New scripts: `rosetta_fixed_contract.py`, `research_rosetta_fixed_scoring.py` and
`summarize_rosetta_fixed_scoring.py`; tests: `test_rosetta_fixed_contract.py`.
This remains private research scripts only; no production consumer, package
dependency, coordinate generator, old threshold or old evidence changed.

Same 48 declared masks: 16 shared source refusals and 32 admitted cases. Score all
64 retained own attempts, 32 retained template proposals and 32 native controls
(128 inputs), including old rejected geometry as diagnostics. Each rounded-heavy
candidate is independently imported into existing PyRosetta 2026.06 with packing
and H optimization OFF; H is rebuilt for that candidate, not carried over.
Scoring imports no normal scientific stack and never calls the packer/minimizer.

The first executable inventory blocked all 128 BEFORE scoring because Rosetta
adds typed virtual PRO NV atoms in its heavy-index range. Preserve that run at
`tmp/pdbfixer-research/rosetta-fixed-run-03/`, including executed harness bytes.
Separate v2 plan inventories physical and typed-virtual atoms separately, with no
physical atom-name exemption or geometry waiver. Other planning/invocation errors
remain retained. Final plan/run: `rosetta-fixed-plan-03/` / `rosetta-fixed-run-04/`
under `tmp/pdbfixer-research/`; all 128 score successfully.

Exact input physical heavy identity/coordinates, source records, source-derived
external links and canonical target internal bonds are checked. Extra physical
atoms are only documented terminal OXT; virtual atoms are retained separately.
Physical heavy atoms stay exact through scoring. Independent readback also checks
score-term sums and identical extra-heavy/virtual/non-target-H coordinates across
each case. Context internal bond orders are NOT independently fully audited;
default Rosetta protonation/unoptimized H is NOT legacy PROPKA/Reduce equivalence.

Selection uses minimum total ref2015 score only among the unchanged source AND
common-gate-feasible own attempts; missing feasible scores block the whole choice.
Selection is committed before native control imports. Reference error is computed
afterward by the independent summarizer, never supplied to ranking. There are
22 one-candidate pools, two two-candidate, seven three-candidate and one empty pool:
only nine cases have a real choice. 31 selections remain; VAL39 CG1 still refuses.

Three selections change, all with lower rounded exact-name held-out error:
1EMV THR11 paired 2.483→0.167 Å, whole 2.000→0.183 Å; THR20 whole 1.977→0.036 Å.
These are correlated masks on two residues of one structure, with previous outcomes
known before the new policy; not prospective/general superiority. Matched own-only
paired n=16 median improves 0.221→0.176 Å; singleton n=15 remains 0.041 Å. Do not
compare these rounded/all-own-feasible metrics to previous raw or mutually-passing
own/template n=8/n=9 denominators. Historical bad 1A3N VAL1 is not fixed by this run.

Partial-graph diagnostics retain all 34 old rounded clash instances across 128
inputs: ten known three-bond/1–4, 24 without a short path in the incomplete graph.
The latter are unresolved, NOT proven nonbonded pairs. Three native refusals still
share VAL39 CG1–C. Other 1–4 contacts include template THR11 OG1–N at 1.989 Å, so
path classification cannot justify blanket exclusion or acceptance. No new gate.

Final verification: 575 focused tests in 25.90 s, targeted Ruff, seventeen-file
package mypy, agent maps and whitespace pass. Independent scoring/old 16-cell
summaries and all 32 old native audits replay byte-identically; old 48-case and
six original/eight mixed evidence records remain valid. All 1,056 local inline
documentation link targets exist. Quality record:
`tmp/pdbfixer-research/rosetta-fixed-docs-quality.md`.
**Resume here:** predeclare an independent-family scoring follow-up, then evaluate
joint atom-fixed combinations with exact small-case search. Richer feasible
sampling, calibrated chemistry-aware sterics, independent context chemistry and
production/public integration remain separate unfinished gates. Existing Rosetta
licenses do not authorize redistribution or exporting an empirical library.

## Earlier checkpoint: Rosetta packer source analysis and installed mechanism probes

User authorized source analysis and use of existing `pyrosetta` environment.
Read [mechanism, ownership, count-pair and implications](../research/rosetta-repack-mechanism.md)
and [portable source/probe evidence](../research/rosetta-repack-mechanism-evidence.json).
Installed PyRosetta 2026.06 identifies Rosetta `1a56185c2592611dec4c9c75ddc9468cd2227c1f`,
not inspected current upstream `eddc15b`. Matching pinned source manifest: 32 files;
raw upstream files stay local under `tmp/pdbfixer-research/rosetta-installed-source-01/`.
No installed environment changes, vendored source, exported empirical tables or
production dependency. Existing installation is not redistribution/license clearance.

Pipeline: import fills missing coordinates from residue internal coordinates;
optional `pack_missing_sidechains` selects residues by missing-atom mask and calls
FULL `pack_rotamers` despite stale comments/docstrings saying rotamer trials.
Tasks constrain residue identity/freedom; rotamer libraries build typed,
backbone-conditioned candidates with extra/proton chi; bump filtering precedes
interaction-graph energies and seeded simulated annealing. Final placement
replaces residues. IncludeCurrent is an option, not observed-atom immobility.

16 predeclared mechanism cells (four retrospective masks × auto-import packing
off/on × preserve-Cβ off/on), fresh single-threaded processes, seed 7: all execute.
Artifacts: `tmp/pdbfixer-research/rosetta-mechanism-{plan,run}-02/`. With automatic
packing off, all eight imports preserve observed heavy atoms but fill missing ones.
With it on, all four singleton imports move observed CB/partners. Explicit repack
moves observations in all eight singleton cells; all eight whole-sidechain cells
preserve observations. Backbone/residue identities unchanged in all 16; target-only
RDKit CIP passes 16/16, not independent global stereo qualification. No restoration
or splicing conceals these contract mismatches. Do not compare their raw RMSDs as
matched generated-only reconstruction accuracy.

Installed option is `packing:preserve_c_beta`, not older-doc `preserve_input_cb`.
It aligns the CA–CB direction without copying measured length/coordinates. Even
enabled, VAL39 observed CB moves 0.001234 Å and CG2 moves 0.435373 Å; THR11 CB moves
0.015855 Å and OG1 moves 0.118334 Å. Diagnostic sets are separately constructed,
not annealer traces. IncludeCurrent gives only one observation-preserving candidate
per singleton set when import packing is off; current includes imputed atoms.

Critical count-pair correction: `cp_half` is **0.2**, not 0.5. Crossover 3 counts
three-bond (1–4) pairs at 0.2; crossover 4 excludes them and counts four-bond pairs
at 0.2. Mode/connectivity/type/options choose behavior. Classic `fa_intra_rep`
(ref2015 weight 0.005) still counts native VAL39 CG1–C; this pair is not globally
ignored. Native raw pair repulsion 4.974661 yields a derived classic weighted
contribution 0.004975. Deliberate distance contractions raise pair repulsion but
are not physical structures or calibrated thresholds. See `rosetta-contact-02/`.

**Resume here:** preserve independent hard geometry/stereo/ownership checks;
first add topology-aware diagnostics, then separately specify any new steric/score
policy. Rosetta's useful lessons are richer feasible candidates, coherent scoring
and cached joint interactions, not wholesale adoption of residue-level repack or
universal 1–4 exclusions. Keep exact enumeration for small cohorts. An atom-fixed
partial-completion adapter and joint/accuracy qualification remain unimplemented.
Two initial CLI/binding protocol blockers are retained; no scientific failures
were removed. Prior frozen comparison/source hashes and production gates unchanged.

Verification: 559 focused tests in 25.15 s; targeted Ruff, seventeen-file package
mypy, agent maps and whitespace pass. Independent saved-snapshot summary replays
byte-identically; old 48-case/32-native and six original/eight mixed evidence
controls remain valid. All 1,050 local inline documentation link targets exist.
Detailed final checks: `tmp/pdbfixer-research/rosetta-docs-quality.md`.

## Earlier checkpoint: rule-selected 48-mask study and native 1–4 gate limitation

`scripts/unfiltered_branched_comparison.py` now separates plan and run. Rule fixed
before execution: first two VAL and first two THR identities in chain-A file order
of 1EMV/1TIM/1DQJ, all four paired/whole-sidechain/singleton masks. All 48 inputs,
original/projection bytes and code hashes are frozen under
`tmp/pdbfixer-research/unfiltered-plan-01/`; proposals never determine enrollment.

Read [selection, complete outcomes and gate caveat](../research/protein-completion-unfiltered.md),
[comparison evidence](../research/protein-completion-unfiltered-evidence.json) and
[post-execution native audit](../research/protein-completion-unfiltered-reference-audit.json).
Run artifacts are in `tmp/pdbfixer-research/unfiltered-comparison-01/`.
All original chain-A ATOM observations are retained, including altlocs/partial
occupancy; no source repair/selection/normalization was used to make admission pass.
These are explicit projections, not intact complexes. 1EMV is homologous to
previous 1FR2; do not count these files as independent families/blind sampling.

Outcomes: 16 shared source-admission refusals (all 1DQJ), 32 admitted masks. Own
common passes 31/32, seeded-template 18/32, zero template execution/protocol
blockers. Both pass 17, own-only 14, template-only one. Template runs even if own
generation fails. Shared passing-subset medians: paired own/template 0.184/0.130 Å
(n=8), singleton 0.055/0.137 Å (n=9). Source refusals remain declared cases, not
scientific generation failures or discarded outcomes. No new empirical data.

Critical result: own singleton 1TIM VAL39 CG1 is rejected while template passes.
Post-execution native-coordinate control passes only 29/32 admitted masks; all
three failures expose the same native CG1–backbone C topological 1–4 contact.
Native distance 2.762756 Å; own 2.767208 Å, error 0.018689 Å; template native error
0.198586 Å but common pass. Native CIP CA=S. The unchanged coarse C–C clash cutoff
is 2.8 Å, so this is a gate limitation, not proof that the closer proposal is
worse or three native stereo defects. **No threshold or old policy was changed.**
The new reference sanity audit is not calibrated external packing validation.

Reporting adds VAL methyl-swap RMSD only when BOTH methyls are generated; fixed
CB/partners and singleton names are never exchanged. THR branches cannot swap.
All actual symmetry scores equal exact-name scores; tests exercise swapped
controls. Truth is post-proposal reporting only, including failed own attempts;
there is no reference-guided reranking or automatic acceptance of rejected data.

Verification: **556 focused tests passed in 24.42 s**, seventeen-file package
mypy and targeted Ruff passed. Tests prohibit reconstruction during planning,
verify frozen inputs, baseline execution despite own failure, generated-only
symmetry and native zero-RMSD gate failure without false stereo claims. Only
research scripts/tests/docs changed; production/core policies remain untouched.

Final read-only validation also passed: 21 frozen source hashes/54 artifacts,
all 48 retained outcomes, byte-identical replay of all 32 native audits and
attempt scores, six old cohort/eight mixed evidence records, agent maps and
1,045 local inline link targets across docs. `git diff --check` is clean.
Final verifier and quality record are retained under
`tmp/pdbfixer-research/validate-unfiltered-final.py` and
`tmp/pdbfixer-research/unfiltered-docs-quality.md`; no reconstruction was rerun.

**Resume here:** declare/test a chemistry-aware 1–4 steric policy with native
and deliberately perturbed controls, retaining genuine new–new/nonbonded clashes.
Version any changed scientific policy; do not quietly loosen frozen gates or
rewrite old evidence. Broader independent families/masks, calibrated metrics,
ILE/multi-chi, licensed data and all production adapters remain open.

## Earlier checkpoint: matched seeded-template scientific comparison

Research scripts now compare 32 **retrospectively selected prior own-accepted**
single-target masks (16 paired, 16 singleton) against the shared seeded-PDBFixer
helper on Reference. Both arms use identical masked bytes/context; held-out truth
enters only after proposals/selection. PDBFixer 1.12.0/OpenMM 8.5.2 passed actual
soft.xml/template/platform preflight without environment repair. Every returned
seed state is retained; no reflection, H placement or downstream refinement.

Read [design, matched-subset results and caveats](../research/protein-completion-template-comparison.md),
[full evidence](../research/protein-completion-template-comparison-evidence.json)
and [corrected aggregate](../research/protein-completion-template-comparison-summary-02.json).
Raw artifacts: `tmp/pdbfixer-research/template-comparison-01/`. New scripts:
`branched_template_baseline.py`, `compare_branched_completion.py`,
`summarize_branched_comparison.py`. All source/resource hashes are retained.

Common physical/source/rounded gates passed own 32/32 vs template 15/32 (paired
6/16, singleton 9/16). These are own-success-filtered controls, not independent
own success-rate evidence. The research comparator accepts measured non-grid χ;
production and own sampling-policy validators do not opt in. The geometric gate
uses own illustrative reference angles, not an independent calibrated standard.
Separate RDKit CIP audits provide the independent chemical stereo check.

Raw exact-name medians including rejected proposals: paired own/template
0.179/0.992 Å, singleton 0.085/0.379 Å. On the SAME mutually passing subset:
paired 0.164/0.107 Å (template more accurate in 5/6 cases), singleton 0.069/0.133 Å
(own lower in 6/9). Do not hide the template advantage on accepted paired masks.
The own whole-sidechain 1A3N VAL1 remains 1.865 Å raw (1.864 rounded), versus
template 1.227 Å raw but rejected geometry. Exact-name scoring is not VAL symmetry-
minimized, and no reference-guided reranking or coordinate repair occurred.

Four template Cα gate failures include nearplanarity: only ONE is confirmed CA=R
by independent rounded CIP. Nine THR proposals have confirmed CB=S, whereas ten
fail the broader R-or-nondegenerate CB gate (one lacks a CB CIP label). Two VAL
name-orientation failures are not physical CB stereocentre inversions. The first
aggregate JSON used ambiguous inversion/count labels; preserve it but use the
corrected summary-02. Original raw proposal evidence was not changed.

Helper retries: 31 cases returned seed 1; 1FR2 THR11 whole sidechain exhausted all
five seeds with a retained D-Cα warning. All returned states/logs are present.
Temporary observed motion <=1.435e-14 Å. One extra OXT per case is recorded but
never inserted into masked proposals. This is NOT the full chirality-safe legacy
prepare/minimize workflow and establishes no new causal stage or production defect.

Verification: **551 focused tests passed in 21.69 s**, seventeen-file package mypy
and targeted Ruff passed. Production defaults, consumers, platform choice,
PDBFixer dependencies and downstream chirality guards remain unchanged. All work
is uncommitted; comparator PDBs are diagnostic proposals, not published acceptance.
All 64 saved-arm proposal audits/scores were independently recomputed without
resampling; all 36 seed snapshots and source/resource/artifact SHA checks passed.
The corrected aggregate reproduces exactly. Fourteen older cohort/mixed evidence
records still replay; agent maps, local links and `git diff --check` passed.

**Resume here:** predeclare broader structures/masks without conditioning on own
success, and add independent calibrated geometry/packing metrics. A full legacy/
refinement comparator is a distinct study; do not infer general superiority,
MD readiness or public replacement from these control results. Empirical data,
ILE/multi-chi, broader chemistry and all production consumer adapters remain open.

## Earlier checkpoint: replay-verified mixed evidence with all failure stages

`protein_completion/mixed_evidence.py` adds distinct schema-1 kind
`experimental-protein-completion-mixed-evidence`. Source, requested targets,
source-derived modes/partners/global THR checks, all individual+joint attempts and
rounded PDB/DatRecord are canonically replayed. Admission refusals retain original
bounded bytes/identities/reason with no fabricated successful plan; raw/rounded
refusals retain all applicable attempts and no accepted candidate files.

Read [mixed evidence contract/lifecycle](../research/protein-completion-mixed-publication.md)
and [portable record](../research/protein-completion-mixed-publication-evidence.json).
Rehashing edited modes/partners/stereo/audits/selection does not bypass replay.
Only known rounded geometry refusals are recorded as scientific materialization
failures; unexpected errors propagate. Hashes establish integrity, not authorship;
an external expected-source SHA can additionally bind consumption.

The writer reuses private artifact reads/durable writes/manifest-last commit,
preflights before staging and reads back staged payloads before committing the
index last. No overwrite/adoption of existing/crash-left/competing paths. Old
cohort/public diffusion/experimental atom-completion readers reject mixed evidence,
and the mixed reader rejects old formats. Public/qualification flags stay false.

Retained publication: `tmp/pdbfixer-research/mixed-publication-01/`, replaying all
eight original `mixed-cohort-02` controls: five accepted, one source-admission
refusal, one no-joint-solution and one rounded refusal. Original source/plans/
attempts and accepted outputs/sidecars were checksum-verified before commit;
no new masks/scientific seeds. The poorer 1A3N VAL1 whole-sidechain 1.864 Å case
remains unchanged, not improved by storing it.

Executed: **544 focused tests passed in 22.08 s**, including generic commit and
old/new evidence tests; seventeen-file mypy and targeted Ruff passed. Tests cover
rehashed semantic edits, mode/partner/stereo binding, canonical bool/int separation,
all refusal stages, source pinning, exact ownership, cross-reader rejection,
single-linked/flat inventory, staging failure, resource bounds and competing
reservations. Production defaults/callers/registries and safeguards are untouched.
All work remains uncommitted; no new public CLI/API/GUI or backend registration.
Eight independent committed replays, five prior PDB byte matches, source/artifact/
retained SHA checks, fixed-partner sidecar exclusion and all refusal stages passed.
All six older evidence records still replay. Agent maps, local Markdown links
and `git diff --check` passed. Historical source hashes remain unchanged.

**Resume here:** storage/replay is not scientific qualification. Broader held-out
mask/structure comparisons, explicit chemical graph/sidechain-stereo coverage,
ILE/multi-chi and genuine licensed empirical statistics remain open. Do not remove
PDBFixer or route production consumers on the strength of these control records;
keep old schema semantics and failed evidence intact.

## Earlier checkpoint: new mixed-method joint source contract

`protein_completion/mixed_cohort.py` now composes explicit standard ALA/SER/free-
CYS, paired VAL/THR and observed-partner singleton modes in one source-bound plan.
Policy: `mixed-canonical-paired-anchored-source-global-thr-beta-v1`. Targets/masks/
modes/partners/global source THR checks are rederived; up to four targets and 81
source-feasible combinations are independently audited for new–source/new–new
clashes. Readback is all-or-nothing with global Cα and complete R THR-Cβ output.

Read [mixed contract and retained controls](../research/protein-completion-mixed.md)
and [current portable record](../research/protein-completion-mixed-evidence-02.json).
The old cohort shares only the bounded search helper; its default plan/format and
valid replay semantics are unchanged. Wrong-kind old evidence calls now reject
explicitly before staging. Mixed results cannot be substituted even into an
otherwise valid all-standard old plan. No mixed evidence format or production
CLI/GUI/backend was added.

Retained execution: `tmp/pdbfixer-research/mixed-cohort-02/`, five accepted outputs
(four real projection cohorts plus a synthetic mixed-choice control), one source
occupancy refusal, one no-joint-solution and one rounding refusal. Every accepted
PDB/DatRecord has exact combined ownership; fixed partners remain observed.
Thirty independent RDKit rounded checks cover THR centres/VAL targets (repeated
centres, not independent samples). Native 1A3N VAL1 full-sidechain RMSD still
equals 1.864 Å: feasibility is not native recovery and truth was not selection.

Initial `mixed-cohort-01/` and its portable record are preserved. Run 02 repeats
the same declared deterministic cases after explicit old-writer type guarding,
not new masks/scientific seeds. Cases are explicit 1UBQ/1FR2/1A3N ATOM projections,
not intact-complex admission; some 1UBQ compositions were unit-tested beforehand.

Verification: **492 focused tests passed in 19.86 s**, sixteen-file mypy and
targeted Ruff passed. All 81 four-paired combinations and free-CYS composition
are tested. Tests retain clashes despite individual passes, selected-rounding
refusal, forged modes/partners/policies, source case/icode identities, resource
limits, edited outputs and old cohort/evidence isolation. Work is uncommitted;
production defaults, PDBFixer dependencies/consumers and chirality guards remain.
Current source/artifact SHA checks, all eight mixed outcomes and five combined
PDB/DatRecord ownership replays passed; run-02 accepted PDB hashes match run 01.
All 32 earlier paired/singleton source readbacks and six old evidence records
still validate. Agent maps, Markdown links and `git diff --check` passed.

**Resume here:** implement a separate mixed evidence schema/reader/publisher,
replaying source-derived modes, fixed partners and all individual/joint/rounded
audits. Preserve old schemas and rejected evidence. Do not infer empirical/native
qualification or registration readiness. Harder chemistry/masks, ILE/multi-chi,
licensed rotamers, scientific comparators and production adapters remain open.

## Earlier checkpoint: separate anchored source policy and rounded singleton output

`protein_completion/anchored_source.py` now adds a distinct source plan binding
full bytes/context/links, the exact observed distal partner, source THR-Cβ checks
and policy `singleton-val-thr-source-fixed-partner-global-thr-beta-v1`. Only one
distal VAL/THR atom can be generated; CB and partner stay observed. Independent
readback checks exact records/ownership/rounding/geometry, global Cα and complete
R THR-Cβ centres. No pre-generation target exemption survives output.

Read [anchored source contract and lifecycle](../research/protein-completion-anchored-source.md)
and [portable evidence](../research/protein-completion-anchored-source-evidence.json).
Shared source assembly has a private explicitly selected admission validator;
default/paired/cohort/evidence entry points keep their old validator and scope.
Paired and anchored plan classes/policy IDs cannot be silently interchanged.
No publication format, public CLI/GUI or production caller was added.

Retained run: `tmp/pdbfixer-research/anchored-source-01/`. All 18 cases were
declared first: 16 singleton masks on explicit 1UBQ A1–72, observed 1FR2 chain A
and 1A3N A1–59 passed raw/rounded/file/sidecar checks; two full-chain-A 1A3N controls
refused partial occupancy without normalization. Original inputs/plans/attempts/
refusals/PDBs/sidecars/audits are retained. RDKit independently checked rounded
THR/target centres 108 times (repeated fixed centres, not independent samples).
Fixed partners remain byte-identical and never acquire generated `.dat` ownership.

Rounded singleton exact-name error is 0.012–0.243 Å, scored only after selection.
This easier mask does not fix/supersede the whole-sidechain VAL1 1.864 Å case.
The source probe gained `--mode anchored`; default paired mode/semantics stay
unchanged. No empirical statistics or MD/native-quality qualification follows.

Verification: **473 tests passed in 17.65 s**, fifteen-file mypy and targeted
Ruff passed. Tests reject coherent partner/observation forgery, wrong policy,
observed D-Cα and wrong/undefined context THR-Cβ, shifted old/new atoms, rounded
O/C-tip clashes and exhausted serials; valid partner ANISOU, sparse TER and CONECT
survive. New RDKit controls independently check rounded singleton stereochemistry.
All 16 singleton PDB/ownership replays and current source/artifact SHA checks pass;
16 earlier paired readbacks and six older private evidence records still validate.
Agent maps, Markdown links and `git diff --check` pass. All work remains
uncommitted; production defaults/guards/callers are unchanged.

**Resume here:** any joint/evidence integration of paired and singleton source
policies needs a new versioned contract and cross-target/global-stereo/readback
audits, not edits to old frozen formats or a larger default whitelist. Broader
chemistry, harder masks, ILE/multi-chi, empirical data and production adapters
remain separate gates. Keep failed evidence and the poorer whole-sidechain case.

## Earlier checkpoint: observed-partner singleton kernel, still typed-only

`protein_completion/anchored_branch.py` adds an explicit one-distal-atom VAL/THR
entry. `geometry.two_angle_tip` algebraically intersects two cones around fixed
CB using retained CA and the other branch, choosing the correct named-volume
hemisphere. No observed atom is moved and no chi sampling/energy minimization is
used. Shared validation explicitly isolates this mask/continuous measured-chi
contract; default paired/kernel/source/cohort/evidence consumers still refuse it.

Read [anchored singleton scope and controls](../research/protein-completion-anchored.md)
and [portable evidence](../research/protein-completion-anchored-evidence.json).
Result revision is `canonical-val-thr-observed-branch-anchor-v1`; policy is
`observed-partner-two-angle-closure-v1`. There is no empirical-library parameter.
Independent audit rechecks both constraints, CA/THR-CB/named-VAL orientation,
ownership, measured chi and clashes. A solvable but clashing cone proposal is
retained as a rejection; impossible/collinear/nearplanar cones refuse without
moving observations or inventing a fallback.

Retained execution: `tmp/pdbfixer-research/anchored-branch-01/`. Sixteen declared
single-atom masks across the same explicit 1UBQ, 1FR2 and 1A3N projections passed
scoped geometry and independent RDKit CIP. Exact-name error is 0.011–0.243 Å;
truth was consulted only after proposal/selection. These easier masks retain CB
and one distal branch: they do not fix or supersede the earlier harder 1A3N VAL1
whole-sidechain 1.864 Å result. Full requests, attempts, masks, originals and
checks are retained, but no singleton PDB/.dat was generated or published.

Verification: **443 tests passed in 15.40 s**; fourteen-file mypy and targeted
Ruff passed. Tests keep synthetic algebraic successes separate from backbone
clash refusals, check the opposite solution despite corrected chi, independent
CIP when either branch is missing, metadata forgery and old API isolation.
All 16 singleton requests/attempts replay, source/artifact SHA checks and agent
maps/Markdown links/`git diff --check` pass. All 16 older paired PDB readbacks and
six earlier private evidence records still validate under the shared validator.
All changes remain uncommitted; production defaults/guards/callers are untouched.

**Resume here:** implement a distinct anchored-source policy with exact fixed-
partner/source/global-stereo and rounded-record binding before any new PDB writer
or evidence/cohort enrollment. Do not silently change paired-source v1 semantics.
More demanding masks, ILE/multi-chi, empirical statistics, comparative science
and all production adapters remain open.

## Earlier checkpoint: separate VAL/THR source path with global THR Cβ

`protein_completion/branched_source.py` now explicitly admits VAL/THR single-
target paired branches, binds source/input THR-Cβ checks/policy into its plan,
inserts only new atoms and independently audits rounded PDB and generated-only
DatRecord. It reuses source parsing/layout/link/graph/insertion/readback helpers;
default source/cohort/evidence target capability remains ALA/SER/free CYS.

Read [new source contract and limits](../research/protein-completion-branched-source.md)
and [portable record](../research/protein-completion-branched-source-evidence.json).
Policy `paired-val-thr-source-global-thr-beta-v1` requires every retained THR CB
to be defined and R (>0.5 Å³ named volume). Only the exact admitted incomplete
target can be pending before generation; no exceptions survive rounded output.
Global Cα checks remain unchanged. VAL naming is not a physical CB stereocentre;
ILE/other chemical centres are not globally qualified. No RDKit import in core.

Retained run: `tmp/pdbfixer-research/branched-source-01/`, all 18 cases declared
before execution. Sixteen paired masks on explicit 1UBQ A1–72, all observed 1FR2
chain-A ATOM records and 1A3N A1–59 passed raw/rounded/file/sidecar/global-THR checks.
Two full-chain-A 1A3N controls refused partial occupancy without rewriting it.
Original bytes, projections, requests, attempts, failures and outputs are retained.
RDKit independently checked rounded coordinates 108 times across accepted files
(repeated fixed centres, not independent scientific samples).

Held-out exact-name RMSD spans 0.027–1.864 Å. Do not conceal the poorer 1A3N VAL1
whole-sidechain result (1.864 Å; retained CB gives 0.281 Å) or call geometry success
native recovery. Truth entered only later reporting; no reranking, FF relax or
threshold adjustment was used. No intact-complex acceptance or empirical-library
qualification follows from these explicit projections.

Verification: **403 tests passed in 15.85 s**, thirteen-file mypy and targeted Ruff
passed. New tests cover wrong/nearplanar/undefined global THR CB despite correct
Cα, exact-chain/icode pending exemptions, policy/context tampering, ANISOU and
serial preservation/exhaustion, rounding-created clashes and all synthetic χ CIP
labels after readback. Original eight source controls reran under `source-adapter-04/`.
Their PDB hashes match the original exactly; all six earlier private evidence
records still replay. Current source/retained artifact SHA checks, all 16 independent
rounded/sidecar replays, agent maps, Markdown links and `git diff --check` passed.
No production caller/CLI/GUI/default backend changed; all work remains uncommitted.

**Resume here:** keep this explicit source mode separate. Any branched cohort or
evidence expansion needs a versioned joint/global-stereo/readback contract, not a
larger default residue set. More demanding masks, singleton anchored completion,
ILE/multi-chi chemistry, real licensed rotamers and comparative qualification stay
open. Preserve existing chirality guards and the visible poorer native-fit case.

## Earlier checkpoint: typed VAL/THR paired-branch kernel

The typed kernel now constructs paired VAL/THR branches with explicit named-Cβ
orientation and independent bond/angle/stereo/clash audits. Natural THR requires
CA=S, CB=R; VAL CB is not a physical stereocentre (two identical methyls), so its
negative named volume is an atom-label convention, not a D/R/S claim.

Read [branched implementation and measured controls](../research/protein-completion-branched.md)
and [portable evidence](../research/protein-completion-branched-evidence.json).
`branched.py` and a research-only RDKit CIP oracle are new. Both branches can be
absent, with or without CB; singleton distal masks/observed-distal closure refuse.
New result revision: `canonical-branched-sidechain-staggered-v1`. No authentic
rotamer data were used; supplying the SER/CYH authority explicitly refuses.

Crucial capability boundary: source/cohort/evidence targets still explicitly
admit only ALA/SER/free CYS. New typed support cannot auto-enroll those consumers
or their internal-graph admission. Production defaults/callers remain unchanged.
The typed THR audit covers its target, not all context/ILE/ligand stereocentres.

Executed: eight declared masks (VAL5/26, THR7/55, paired branches with/without CB)
on the explicit 1UBQ A1–72 projection under `tmp/pdbfixer-research/branched-controls-01/`.
All eight selected a geometry pass; 12/24 proposals passed geometry and all 24
passed independent CIP/named-volume controls. Selected exact-name RMSD ranges
0.027–0.472 Å; truth was read only for later reporting, never selection. These
are four-residue controls, not template superiority or broader qualification.
Original source, predeclared plan, full requests, attempts, refusals and scores
are retained. Source/cohort publication of these new residues remains disabled.

Independent RDKit checks (2025.09.6, Python 3.12.13) agree on eleven complete
reference residues and every synthetic χ1. Reflecting only THR's second branch
inverts CB R→S without changing CA=S or primary χ; both audits reject it, including
when it is fixed observed geometry on a no-op. Coordinate-transport atom maps are
cleared before independent CIP assignment. No RDKit import was added to the core.

Verification: **374 tests passed in 13.53 s**, twelve-file mypy and targeted Ruff
passed. Original ALA/SER controls reran under `own-kernel-controls-05.json`;
their coordinates/audits/selection match the original exactly. All six earlier
private evidence records still replay successfully. Current-source/raw-record SHA
checks, agent maps, Markdown links and `git diff --check` passed. Historical
evidence/hashes are retained. All changes remain uncommitted.

**Resume here:** qualify an explicit source-level THR Cβ invariant and rounded
paired-branch output before deliberately expanding source/cohort/evidence targets.
Add broader structural controls; do not infer global stereo, native packing or
public readiness from these cases. Anchored singleton completion, ILE/multi-chi,
licensed empirical data and all production consumer adapters remain open.

## Earlier checkpoint: private replay-verified evidence and manifest-last commit

`protein_completion/evidence.py` now snapshots and reads cohort source/targets/
full plans/all individual+joint attempts/audits, with candidate PDB and canonical
generated-only `.dat` only for accepted rounded results. Raw/rounded scientific
refusals retain evidence without accepted structures. The distinct experimental
kind keeps public/scientific qualification flags false; both existing public
diffusion and experimental atom-completion readers reject it.

Read [private evidence contract and lifecycle](../research/protein-completion-publication.md)
and [portable record](../research/protein-completion-publication-evidence.json).
The reader verifies bounded flat regular artifacts/digests, then independently
replays source planning, ALL attempts/selection, materialization and combined
readback. Rehashed tampering still refuses; canonical JSON comparison avoids
Python's bool/int equality loophole. `.dat` parsing remains owned by DatRecord.

Publication reuses the private-file reader and shared manifest-last commit helper.
Preflight precedes staging; existing/crash-left/competing paths are never adopted
or overwritten. Staging I/O failure preserves written diagnostics; generic commit
rollback may remove its own transferred entries, never caller inputs. Unknown
runtime/materialization failures propagate rather than becoming scientific passes
or fake geometry refusals. No public CLI/backend/GUI routing was added.

Retained publication: `tmp/pdbfixer-research/cohort-evidence-01/`. All six previous
`source-cohort-01` controls were checksum-verified and committed/read back: four
accepted and two rejected, preserving unchanged attempt and output content. No
new masks/scientific seeds; originals and older hashes remain historical.

Executed: **341 focused tests passed in 13.82 s**, including generic commit tests;
eleven-file package mypy and targeted Ruff passed. Public-reader compatibility
was actually run in the scientific interpreter (system Python lacks OpenMM),
without introducing a skip or repairing environments. Six independent evidence
readbacks, current/retained source-artifact SHA checks, agent maps, Markdown links
and `git diff --check` passed. All work is uncommitted.

**Resume here:** evidence integrity is not scientific qualification. Authentic
licensed rotamer data/independent empirical evaluation and broader chemistry stay
open. Before any production consumer adapter, extend capability and regression
coverage deliberately; preserve source ownership, rejected evidence and all
existing chirality guards. Do not relabel these control-only bundles as public
diffusion results, native packing or proof that PDBFixer can be removed.

## Earlier checkpoint: bounded multi-target control and joint readback

`protein_completion/cohort.py` adds immutable source-bound plans for 1–4 targets,
retained per-target/joint audits, deterministic enumeration of <=81 feasible
staggered combinations, and one-pass/all-or-nothing PDB/DatRecord materialization.
Every target is re-audited against every other target's new atoms before and
after rounding. Source bytes/ownership/links remain immutable. Public single-
target admission never inherits cohort missing-CB exemptions.

Read [cohort contract and limits](../research/protein-completion-cohort.md) and
[portable record](../research/protein-completion-cohort-evidence.json).
`source.py` now shares insertion and physical-readback helpers with the cohort;
its existing public-internal entry point retains the original single-target gate.
No default, production caller/registry, geometry threshold or dependency changed.

Retained evidence: `tmp/pdbfixer-research/source-cohort-01/`. Three explicitly
projected 1UBQ A1–72 multi-target cases (2 ALA, 3 SER, 4 mixed targets) passed joint
raw/rounded/file/sidecar readback. Synthetic controls retain a two-SER rejected
combination, a two-ALA no-joint-solution refusal despite individual passes, and
whole-output refusal when rounding breaks one target. All 81 four-SER combinations
are independently exercised in tests. These are not empirical/native-quality
benchmarks or full-source acceptance. Eight single-target outputs reran under
`source-adapter-03/`; historical evidence is retained rather than overwritten.

Executed: **288 focused tests passed in 12.49 s**; ten-file package mypy and
targeted Ruff passed. Agent maps, current source/artifact SHA checks, all six
declared cohort outcomes, eight-mask single-target PDB byte parity, Markdown
links and `git diff --check` passed. All changes remain uncommitted; no public API/CLI/backend
or filesystem publisher was added. The cohort is explicitly staggered-control-
only: no empirical-library parameter or authentic Dunbrack dataset is implied.

**Resume here:** design a separately audited evidence/publication contract before
any public adapter. Preserve all raw and joint failures, full source binding and
manifest-last acceptance. Broader chemistry/empirical packing and licensed real
rotamer data remain independent gates. Do not treat bounded feasibility search
as a physically scored packing model or enable a production replacement.

## Earlier checkpoint: strict declared-link and graph admission

The private source adapter now parses exact LINK/SSBOND endpoints and validates
supported ALA/SER/CYS internal CONECT/LINK against a shared canonical connectivity
graph. Explicit peptides must agree with observed same-chain file-order/TER-aware
C–N evidence. Context CYS SG disulfides require 1.8–2.3 Å and 1:1 partners across
SSBOND/LINK/CONECT. Target disulfides, nonidentity symmetry, absent endpoints,
bad geometry, unknown crosslinks and unsupported graph types still refuse.

Read [implemented link scope and controls](../research/protein-completion-link-adapter.md)
and [portable record](../research/protein-completion-link-adapter-evidence.json).
`links.py` is a new strict metadata parser, not a coordinate/CIF reader. The
supported target graph is shared with the existing candidate auditor; no bond
orders, FF terms or new geometry thresholds were introduced. Internal graph
metadata for other residue types now refuses, rather than being retained without
chemical validation as in the earlier slice.

Retained real controls: `tmp/pdbfixer-research/source-links-01/`, explicitly
projected 1DQJ ALA9 (CB withheld), CYS23 and CYS88. Baseline, deposited SSBOND and
SSBOND+canonical-target-CONECT completed with identical generated coordinates;
nonidentity symmetry and disulfide-target selections refused. Original bytes,
predeclared outcomes, inputs, requests, attempts, outputs and sidecars are retained.
This is not full-1DQJ acceptance, native-quality superiority or automatic
disulfide-detection qualification. All eight prior 1UBQ masks also reran under
`source-adapter-02/`, with separate evidence rather than overwriting old records.

Executed: **266 focused tests passed in 9.84 s**, nine-file package mypy and
targeted Ruff passed. Current source/retained-artifact hashes, eight-mask PDB
byte parity, agent maps, Markdown links and `git diff --check` passed.
Production defaults/callers/registries remain unchanged;
all work is uncommitted. No authentic rotamer data were loaded.

**Resume here:** extend only with declared graph authority and independent
tests. Multiple-target composition and evidence/publication remain separate
gates; authentic Dunbrack artifact/licensing and scientific evaluation remain
open. Do not enable a public replacement based on these admission controls.

## Earlier checkpoint: private source-bound PDB adapter

`protein_completion/source.py` now derives full observations/masks/peptide links
from source bytes, rebinds requests before materializing, inserts only new atoms,
and independently re-audits rounded PDB records and generated-only `DatRecord`.
Source and output protein Cα checks reject D/nearplanar/undefined centres outside
the permitted pre-generation missing target CB. No production consumer switched.

Read [source adapter scope and refusals](../research/protein-completion-source-adapter.md)
and [portable record](../research/protein-completion-source-adapter-evidence.json).
Eight real-fixture masks passed independent source/rounded-output/sidecar readback
under `tmp/pdbfixer-research/source-adapter-01/`. This is the explicitly retained
1UBQ ATOM A1–72 projection, not full-source acceptance: partial occupancy and
MASTER in the deposited source remain refused, without occupancy rewriting.
The source-bound controls use the unchanged staggered sampler; no authentic
Dunbrack rows were obtained. The publisher was retried and timed out again.

Current boundary: canonical heavy-atom protein/water only, one target. No
MODEL/altloc/measured H, PTMs/ligands/ions, LINK/SSBOND, target intra-residue CONECT
or nonlocal/nonpeptide external CONECT. Nearby CYS contexts conservatively refuse
CYS target completion. Other source records remain byte-identical, including
valid ANISOU and bare TER. Existing parser/formatter are reused; the adapter
itself imports the scientific stack but kernel imports remain separate.

Tests explicitly reject tampered context/links/records, global D/undefined Cα,
serial exhaustion, unsupported chemistry and rounding-created clashes. Raw χ
checks remain exact; only bounded PDB readback uses a separate 0.15-degree χ
tolerance. All original geometry/clash/chirality thresholds remain unchanged.
Executed: **231 tests passed in 9.08 s**, eight-file package mypy and targeted
Ruff passed. Agent maps, source/retained-artifact SHA checks, Markdown links and
`git diff --check` passed. All changes remain uncommitted; no push or registration.
No output filesystem publisher exists: the private service returns bytes;
only the research harness explicitly retains artifacts. Earlier source hashes
remain historical and are not overwritten to match the extended validator.

**Resume here:** authentic artifact/licensing and empirical qualification remain
open. Broaden the source adapter only with independent link/chemical-graph tests;
then design multiple-target composition and an evidence/publication contract.
Do not register a production backend based on these source/transport controls.

## Earlier checkpoint: single-chi adapter implemented, real-data gate open

The user's continuation now has an offline, checksum-checked traditional
single-chi table reader, periodic phi/psi probability/circular-moment lookup,
observed-link-based backbone angles, optional authority-driven proposals and
independent evidence rechecking. The staggered control remains separate. No
production PDBFixer consumer has been replaced.

Read [rotamer adapter and open data gate](../research/protein-completion-rotamer-adapter.md)
and [portable evidence](../research/protein-completion-rotamer-adapter-evidence.json).
**No genuine Dunbrack statistics were downloaded or sampled.** Publisher pages
time out on this host; the inspected old mirror ships restrictive academic terms.
The publisher's indexed page announces modern CC BY 4.0, but no exact data
artifact/license lineage has been approved. Do not infer that the old terms
revoke the modern license, or that modern landing-page text clears an arbitrary
mirror. No registration/license was accepted or data silently substituted.

Current implementation surfaces:

- Added `protein_completion/rotamers.py` and `backbone.py`.
- Extended private contracts/service/audit with optional `RotamerEvidence` and
  model digest; empirical result revision is `canonical-sidechain-rotamer-adapter-v2`.
- Added `tests/test_protein_rotamers.py` and independent Biopython sign checks in
  `tests/test_protein_rotamer_conventions.py`.
- Added `scripts/probe_rotamer_authority.py` (pages/terms only) and
  `scripts/probe_rotamer_transport.py` (explicit synthetic tables only).
- No CLI flag, registry entry, dependency, fixture, existing threshold, public
  diffusion schema or production adapter changed.

Executed: **185 focused tests passed in 8.18 s**; seven-file package mypy and
targeted Ruff passed. Agent maps, Markdown links, current-source/raw-record hash
checks, original-control parity and `git diff --check` also passed.
The synthetic transport ran four requests / 12 candidates
on two authored SER/CYH grids (1,296 cells each), all scoped passes with exact
observations and independent evidence readback. These are **format/transport
controls, not empirical Dunbrack scientific seeds**. Signs of real-source
phi/psi/chi agree with independent Biopython. Original eight-mask staggered
coordinates, audits, status and selection match the retained first-slice run
exactly; do not replace the earlier evidence or count the rerun as a new arm.

Raw evidence: `tmp/pdbfixer-research/rotamer-authority-01/` (all three publisher
page timeouts plus retrieved legacy terms/description), `rotamer-transport-01/`
and `rotamer-transport-02/` (synthetic, current source hashes in latter), and
`own-kernel-controls-04.json` / `own-kernel-summary-04.json` (control rerun).
Earlier first-slice hashes are historical and must not be overwritten to match
the extended sources. All new science remains private with no PDB publication.

**Resume here:** obtain an exact licensed publisher artifact or explicitly
reviewed local copy, ideally Dunbrack 2010 traditional 5% SER and free-CYS
categories; verify original headers, raw SHA, variant/subclass, terms and
attribution before extraction/bundling. Feed that pinned local authority to the
reader and crosscheck real entries independently. Then predeclare real-library
comparisons and broader residue/geometry work. The reader currently admits SER
and explicit CYH, requires complete external 36x36 grids, refuses missing
neighbours/cells instead of falling back, and stores sigma without subrotamer
sampling. No multi-chi, coupled packing or public capability is implied.

## Earlier checkpoint: first internal-coordinate implementation

After the usage audit, the user authorized starting the replacement work. The
first private own-code slice now exists at `src/dvbfixer/protein_completion/`:
ALA/SER/free-CYS sidechains, immutable exact-identity observations, explicit L-Cβ
construction, three staggered chi proposals, independent scoped geometry/stereo/
steric rejection, generated-only candidates and retained failed attempts. It
uses no PDBFixer/OpenMM/template coordinates. **No production consumer is switched.**

Read [first-slice implementation and remaining gates](../research/protein-completion-first-slice.md)
and its [portable evidence](../research/protein-completion-first-slice-evidence.json).
Eight predeclared masks on the labelled 1UBQ 1–72 projection have scoped geometry
passes (16/20 passing proposals); selected SER RMSD reaches 2.459 Å. This is
feasibility, not native quality, broad qualification or empirical rotamer support.
CYS has synthetic controls only. Missing backbone, OXT/cap creation, other
residues, distal closure, linked sidechains and unsupported context refuse.

Changed/new surfaces in this implementation continuation:

- Five private package files under `src/dvbfixer/protein_completion/`.
- `tests/test_protein_completion.py`, `scripts/benchmark_protein_completion.py`.
- `docs/research/protein-completion-first-slice.md` and evidence JSON.
- Research entry/contract/context in `docs/agent/`; current plan/research links.
- No CLI, public schema, dependency, version, existing fixture or production
  threshold change. No generated CLI/GUI specs need regeneration.

Verification: 131 focused tests passed in 7.72 s in the existing general science
runtime with project-local pytest; targeted Ruff and package mypy passed.
The new private kernel plus inventory suite separately passed 70 tests in 2.80 s
in the base interpreter without OpenMM. Local Markdown links, agent maps,
executed-source/raw-evidence digest assertions and `git diff --check` passed.
Base-interpreter broader collection/tests lacked OpenMM, then the same relevant
suite passed in the science interpreter; no skip/resource repair was introduced.
The new kernel also passes a separate `python -S` dependency-free import test.
Full suite not rerun and historical damaged prefixes remain untouched.

Raw immutable requests:
`tmp/pdbfixer-research/own-kernel-controls-03.json`. Source hashes in the portable
record match the executed modules. `own-kernel-controls-01.json` and `-02.json`
are retained from earlier milestones, not independent scientific cohorts. Final
admission includes the absent sidechain's reach in the metal exclusion envelope;
the eight source controls contain no metals and their outcomes remain unchanged.

**Resume here:** pin and evaluate a licensed empirical rotamer/geometry authority,
then extend residue/defect coverage under the T1–T14 matrix. Before any integration,
add full-source/digest/context/link audit and global chirality checks. Never equate
the typed kernel's selected-residue audit with an intact published-PDB audit.
Do not route all nine production calls to this restricted first slice. Original
causal experiments and inventory below remain historical context, not a claim
that the replacement is complete.

## Documentation and evidence plan

- [x] This plan/live handoff (execution order, decisions, continuation).
- [x] `docs/research/pdbfixer-reconstruction.md`: source map, causal experiments,
  feasibility and limits of a custom rebuilder.
- [x] `docs/research/pdbfixer-reconstruction-evidence.json`: portable measured
  results and source/resource identities, not an assertion of general validity.
- [x] Reproducible private research scripts and focused tests.
- [x] Validate local links, measurements and focused tests; record exact results.

Link tracking: [report](../research/pdbfixer-reconstruction.md) and
[evidence](../research/pdbfixer-reconstruction-evidence.json) exist; no pending
document links.

## Work packages and exit conditions

1. Pin an intact scientific runtime, PDBFixer/OpenMM source and resources.
2. Map detection, residue/atom construction, template fit, soft-force minimization,
   collision dynamics and final topology transfer. Separate whole gaps from
   within-residue completion.
3. Retain stage coordinates, signed volumes and first failure; controlled ablations
   must distinguish initialization from minimization/dynamics. Use Reference in
   research only. Do not invoke raw `addMissingAtoms` outside the existing seeded
   rebuild wrapper.
4. Map OpenMM hydrogen placement separately; reproduce or explicitly bound the
   known defect, inspect force terms and existing-H behavior, test interventions.
5. Assess an independent internal-coordinate/rotamer/constraint architecture with
   explicit scope, underdetermination, chemistry authority, ranking and rejection.
6. Deliver intelligible measured findings, not just a literature plan. Unresolved
   causal claims stay unresolved with an exact next experiment. Never relabel a
   mitigation as a cause or a geometrically valid conformation as native truth.

## Safety and continuation

Original licensing worktree is untouched. Keep inputs immutable and all raw failed
attempts. Research output lives in this worktree under `tmp/`; portable summaries
and scripts are tracked. No environment repair, backend replacement, public
capability change or upstream publication is authorized. Seeds do not guarantee
cross-platform coordinate identity. All original zero-D and identity guards stay.

## Executed checkpoint: meaningful causal results obtained

The report distinguishes three mechanisms: reflected SVD initialization; hybrid
fixed-anchor failure even after determinant correction; and L-to-D inversion in
the first optimizer stage without dynamics. The N/CA/C frame + wall sidechain
intervention has no final D in six declared cases. The wall-protected gap still
has nine limited geometry findings and is **not a successful loop repair**.
The independently derived L-Cβ primitive passes six withheld-atom controls.
The tested SER HG collapse is causally linked to one zero-LJ attractive Coulomb
pair in an already cramped fixed-heavy environment, not just a short iteration
budget. These are intelligible narrow findings, not general backend acceptance.

### File ownership

- Research scripts: `scripts/research_pdbfixer.py`,
  `scripts/research_chiral_completion.py`,
  `scripts/research_pdbfixer_hydrogen_pair.py`,
  `scripts/audit_pdbfixer_research.py`.
- Tests: `tests/test_pdbfixer_research.py`.
- Portable docs: this file, report and evidence JSON; plan index updated.
- Narrow `.gitignore` entry excludes only `tmp/pdbfixer-research/`.
- No `src/`, argparse, public protocol, production threshold or version changes.
  No new source fixture: existing fixture manifest needs no regeneration.

### Retained execution state

All paths below are relative to this worktree, not the original licensing tree.

- `tmp/pdbfixer-research/run-01/`: heavy instrumentation error retained;
  do not count its empty `final_d` lists as scientific passes.
- `run-02/`: corrected initial six heavy controls and 20 hydrogen placements.
- `run-03/`: six heavy controls with optimizer callback volumes; 20 hydrogen arms.
- `sidechains-01/`: initial six-site three-arm cohort; exposed hybrid-anchor issue.
- `sidechains-02/`: all 24 six-site arms, including backbone-only proper fit + wall.
- `h-pair-01/`: 15 pair ablations from saved `run-02` initial systems/coordinates.
- `executed-scripts/`: exact sampling-script bytes before import-only cleanup.
- `provenance.json`: pinned upstream/installed byte comparisons and Amber/H hashes.
- `test-tools/`: project-local pytest/Ruff only; no external environment was changed.

`run-03` and `sidechains-02` are the heavy outcome sources in the portable record;
earlier reruns are not independent scientific samples. Pair results retain their
original `run-02` ownership. All scientific requests use source projections, not
altered occupancies or hidden input repair. Observations remain fixed.

### Reproduction commands

Use an intact PDBFixer 1.12.0/OpenMM 8.5.2 interpreter as `$PY`; exact original
prefix is recorded below for this host only. Run from the branch's worktree.
Each output destination must be new. These scripts intentionally preserve bad
research outputs; **do not feed their `returned.pdb` into production as accepted
structures**. The seeded wrapper is called with `max_attempts=1` to expose one
attempt, so its existing exhausted-retry warning is expected in failed arms.

```bash
export PY=/private/var/folders/rj/5v2182513dx99_nlm140sfl00000gn/T/opencode/dvbfixer-apple-cc94-mps-v3/refinement/bin/python
PYTHONPATH=src "$PY" scripts/research_pdbfixer.py --output tmp/pdbfixer-research/repeat-main
PYTHONPATH=src "$PY" scripts/research_chiral_completion.py --output tmp/pdbfixer-research/repeat-sidechains
"$PY" scripts/research_pdbfixer_hydrogen_pair.py --run tmp/pdbfixer-research/repeat-main --output tmp/pdbfixer-research/repeat-pair
PYTHONPATH=src "$PY" scripts/audit_pdbfixer_research.py --run tmp/pdbfixer-research/repeat-main --sidechains tmp/pdbfixer-research/repeat-sidechains --pair tmp/pdbfixer-research/repeat-pair --output tmp/pdbfixer-research/repeat-audit.json
PYTHONPATH=src:tmp/pdbfixer-research/test-tools "$PY" -m pytest tests/test_pdbfixer_research.py tests/test_seeded_atom_rebuild.py tests/test_ffutils_chirality.py -q
```

Provenance was additionally obtained by comparing the installed-source snapshots
with the two commit-pinned raw GitHub URLs recorded in the evidence. This step
uses the network; scientific scripts do not require downloading templates or ML
weights for these canonical cases. Import sorting/unused-import cleanup after
execution changed source digests, not scientific algorithms; exact executed
script bytes remain under `executed-scripts/`.

### Verification and documentation quality

- Focused tests: **20 passed in 11.67 s** in the intact CPU runtime, with pytest
  loaded from project-local tools. Initial attempt in the general Standard runtime
  stopped because pytest was absent; no test was silently skipped.
- Ruff on four scripts and new tests: **passed** (only pre-existing obsolete
  `UP038` configuration warning).
- Agent documentation-map check, local Markdown link checks, portable-evidence
  count/source/fixed-coordinate assertions and `git diff --check`: **passed**.
- Independent readback audited all 30 heavy arms, including every failed D case;
  maximum observed movement was 0.0 nm in every arm.
- Documentation scope: three outcome/plan documents, four scripts and a focused
  test file complete. The report has algorithms, numerical controls, negative
  results, limitations, citations and a concrete next implementation boundary.
- Full scientific suite not run: no production code changed, and the intact
  isolated runtime is intentionally narrower than Standard. No general test or
  GPU/platform readiness claim is made.

## Next-session handoff: do not repeat completed research as a new plan

### User-directed priority update: usage and replacement boundary

Before implementing a complete own builder, inventory actual source call sites,
state dependencies and molecule classes. Deliver a separate usage/requirements
note and a reproducible source inventory, distinguishing stock upstream template
coverage from DVBFixer pass-through/other-engine support. Define minimal
coordinate-kernel replacement versus complete PDBFixer removal, and map required
tests to current coverage without claiming missing tests already pass.
No production replacement is authorized by this documentation/research slice.

- [x] Source-level direct consumer inventory and manual indirect-route review.
- [x] Offline stock-template/chemistry scope probe, with no CCD downloads.
- [x] Own-code functional requirements, molecule routing and concrete test matrix.
- [x] Pure inventory regressions, focused verification and updated handoff.

Outcome: [usage and replacement scope](../research/pdbfixer-usage-and-replacement-scope.md)
and [source/template inventory](../research/pdbfixer-usage-inventory.json).
`scripts/inventory_pdbfixer_usage.py` records 11 constructors and nine seeded
rebuild sites, not just the five historical preparation sites. Its optional
offline probe found 30 stock templates (20 amino acids, eight DNA/RNA nucleotides,
ACE/NME) and detected the withheld heavy atom in all 30 discovery controls.
No reconstruction coordinates were generated in this slice. Keep lists include
I/DI even though their stock coordinate templates are absent; unknown HETATM
silence is not evidence of molecular completeness.

The note defines F1–F9 for the minimal heavy-coordinate kernel, the additional
scope required to remove PDBFixer entirely, and T1–T14 acceptance families with
existing regression assets and explicit missing coverage. Existing nonprotein
Class A reconstruction, OpenMM H placement and whole-gap modeling stay separate.

Verification: `PYTHONPATH=src python -m pytest tests/test_pdbfixer_usage_inventory.py -q`
passed **6 tests in 2.51 s**; targeted Ruff and agent documentation-map check
passed. These are inventory regressions, **not six scientific builder successes**.
The earlier 20 scientific tests are retained historical results: the old temporary
CPU interpreter now fails during encoding/bootstrap initialization and was not
repaired. Only the inventory/discovery probe ran in the existing general
`/Users/lilpoundcake/micromamba/envs/dvbfixer/bin/python` (PDBFixer 1.12.0,
same source SHA); that does not establish its optimizer/resource readiness.

Reproduction of this inventory slice (use a new output path):

```bash
python scripts/inventory_pdbfixer_usage.py --output tmp/pdbfixer-research/new-static-inventory.json
/Users/lilpoundcake/micromamba/envs/dvbfixer/bin/python scripts/inventory_pdbfixer_usage.py --probe-templates --output tmp/pdbfixer-research/new-template-inventory.json
PYTHONPATH=src python -m pytest tests/test_pdbfixer_usage_inventory.py -q
```

**Next action:** choose the explicit v1 chemistry/defect strata and freeze the
T1–T14 subset for them before implementing the sidechain service. Protein-only
completion cannot replace all nine callers unconditionally while nucleic acids,
terminal defects, mutations, temporary FF completions and retained linked
components remain unresolved. No public/default change is approved.

1. Read the report/evidence first. Preserve both causal paths and the negative
   result that proper fitting alone can worsen partial THR/ILE initialization.
2. If developing the private builder, start with canonical sidechains and exact
   observed-atom masks; use explicit stereochemical internal coordinates. Add
   sidechain stereocentres, rings, authoritative residue-dependent geometry and
   constrained rotamer sampling before claiming a usable rebuilder.
3. Predeclare a source-disjoint cohort with 1–2 missing atoms, whole sidechains,
   branched residues, aromatic rings, proline, disulfides, retained partner/ligand
   context, insertion codes and case-sensitive chains. Score held-out coordinates
   only after generation; keep symmetry-aware scoring distinct from identity.
4. Compare all proposals under unchanged independent geometry/stereo/fixed-atom
   gates, with failure/runtime counts and candidate ranking. Compare whole gaps
   to full legacy MODELLER separately. No threshold tuning on the test outcomes.
5. Reproduce Reference mechanisms on the actual production platform before
   proposing a production fix. A proper-fit-only patch is explicitly insufficient.
6. Hydrogen continuation: control the source OG/OXT geometry and separately
   investigate sibling-H geometry. Pair-force ablations are explanatory tools,
   never silent production parameter changes.
7. Any public backend/CLI enablement, broad chemistry expansion or release is a
   separate approval/qualification gate. Do not remove retry, chirality barriers,
   final `assert_all_l`, hydrogen repair or independent geometry rejection.
