# Mac reconstruction comparison continuation: execution evidence

Date: 2026-11-11. Private research, not public backend qualification.
Protocol: [controlled continuation](../plans/mac-reconstruction-comparison-continuation.md).
Historical studies: [feasibility log](mac-reconstruction-feasibility-2026-11-11.md).

Objective clarification (2026-11-12): the user requires usable completion of real
missing regions, not recovery of hidden crystallographic conformations. RMSDs
below remain descriptive, never acceptance gates. Small synthetic panels are
workflow controls, not representative production qualification. The
[natural-gap follow-up](../plans/mac-reconstruction-real-gaps.md) prioritizes real
large incomplete structures, 1HE8 and a precise audit of why bounded joint
generated-plus-anchor minimization did not finish its repair.

## Full seeded MODELLER proposal baseline

Thirty proposals were declared before execution, five historical masks times
seeds 7/42/101/102/103/104. This is a retrospective baseline on already inspected
cases, not an independent benchmark. The private LoopModel adapter uses one
initial model and one fast loop output, not the default production two loops.
It reuses production chirality restraints and gap-only loop MD, but retains the
native proposal **before** any chirality reflection. Initial MODELLER optimization
can move observations. This is not the previous complete-proposal CG-only control.

| Projection / hidden residues | Source-restored raw geometry | Final custody/geometry/L | Median final heavy RMSD | Median final backbone RMSD |
|---|---:|---:|---:|---:|
| 1UBQ A1–72 / A10–13 | 0/6 | 6/6 | 1.497 Å | 0.678 Å |
| 1A3N A1–59 / A30–37 | 0/6 | 5/6 | 6.720 Å | 4.528 Å |
| 1FR2 A3–85 / A24–31 | 0/6 | 6/6 | 5.261 Å | 3.883 Å |
| 1A3N A1–59 / A30–45 | 0/6 | 3/6 | 12.907 Å | 11.940 Å |
| 1FR2 A3–85 / A24–39 | 1/6 | 5/6 | 11.663 Å | 10.567 Å |

RMSD medians describe **all six outcomes including geometry rejections**; rejected
structures are not accepted as scientifically usable accuracy results. No
reference-based candidate selection was performed. There are six distinct final
PDBs per case. All 30 final outputs passed independent global L readback and
byte-identical observed ATOM-record checks; five still fail geometry. Twenty-nine
source-restored proposals required the same existing common generated-only CPU
repair (1000 iterations, optimizer seed7, zero restarts, no preconditioning or
anchor release); the one clean proposal was not rebuilt.

The maximum native observed-atom displacement after observation-only fitting is
5.26/10.59/10.06/14.76/6.89 Å respectively. Thus source restoration is a substantive
intervention, not a negligible formatting step. Its induced junction damage is
retained in the native/restored audit reports. These outcomes do not justify a
general claim that MODELLER is inferior: they concern this explicit one-loop
protocol, projected inputs and restoration policy. Model families and structures
are not independent validation units here. The historical 16-residue ML pilots
have only seed7; comparing them with six-seed MODELLER medians is unmatched.

### Retained technical failure and recovery

Attempt01 finished `model.make()` for 1UBQ seed7 and then failed serializing the
opaque MODELLER `Values` object in `pdfterms`. The guard stopped before the other
29 jobs. Original outputs and partial metadata remain unchanged.

Attempt02 explicitly recovered the unchanged saved native loop PDB (SHA256
`9fe214f1c3f3a794e86d87c7f00025104e9b639ffee43a253166d73f82a2e7b3`),
validated its identity/inventory and performed independent downstream assessment.
Seed7 was **not resampled**. Its opaque loop status/score metadata was not
reconstructed or inferred. Other 29 jobs were first executions. Serialization
now records selected scalar metadata, and pre-serializes before exclusive writes
so unsupported values cannot create a new partial JSON artifact.

Evidence:

- `tmp/mac-reconstruction-comparison-2026-11-11/modeller-attempt01/` — original failure.
- `tmp/mac-reconstruction-comparison-2026-11-11/modeller-attempt02/plan.json` — frozen settings, input/code hashes and recovery amendment.
- Same directory `progress.json`, `completed.json` and `execution/*/` — all 30 native/raw/audit/repair outcomes and execution metrics.
- `tmp/mac-reconstruction-comparison-2026-11-11/modeller-independent-summary/` — independent final readback, per-case ranges and stage-aware summary.

## Prospective expansion: preparation is not execution

New masks were frozen before sampling:

1. 1UBQ A1–72, hidden A1–4: terminal mask passed source/reference checks.
2. 1UBQ A1–72, hidden A10–13 and A46–49: two-gap atomic mask passed checks.
   Spatial interaction is not assumed merely because two gaps exist.
3. 1A3N A1–59 plus B1–146, hidden A30–37: refused on partial occupancy at B2.
   No occupancy changes, source repair, automatic replacement or sampling occurred.

Evidence: `tmp/mac-reconstruction-comparison-2026-11-11/expansion-preparation-attempt01/`.
An additional separately declared barnase/barstar admission control, 1BRS A3–110
plus D3–89 with hidden A83–90, also refused partial occupancy (D14). Its frozen
mask/source and refusal are under
`tmp/mac-reconstruction-comparison-2026-11-11/1brs-interface-admission-attempt01/`.
It does not replace the failed 1A3N cell and was not sampled. Do not repeatedly
move masks or trim partner observations until a nominal interface cell passes.
An admitted mask does not establish backend consumer transport, public admission
or scientific accuracy. Existing single-gap smoke builders remain unchanged;
the separate private expansion encoder uses exact whole-gap ownership, observed
templates and zero-coordinate RFD3 dummy atoms. Templates contain no hidden native
coordinates. Both pair members must come from the same sampler invocation.

Six new seed7 pilots (two masks × three backends) are declared under
`tmp/mac-reconstruction-comparison-2026-11-11/expansion-seed7-pilots-attempt01/`.
Settings remain 200 requested timesteps per backend, explicit RFD3 MPS seeding,
serial pressure-v2 guards and no generic unsupported-op CPU fallback or remote
MSA. RFD3 requested timesteps are not automatically measured diffusion calls.
Native, restored and necessary-only common-repair stages remain distinct.
All six pilots subsequently completed. Independent readback verified exact
observed ATOM records, full requested inventory and global L-chirality; **6/6**
passed final geometry/custody. Only the Protenix terminal proposal was raw-clean;
the other **5/6** required common generated-only CPU repair.

| Mask, seed7 | Backend | Final heavy RMSD | Final backbone RMSD | Repair movement RMS |
|---|---|---:|---:|---:|
| N-terminal A1–4 | RFD3 | 1.225 Å | 0.749 Å | 0.888 Å |
| N-terminal A1–4 | Boltz | 1.245 Å | 0.790 Å | 0.639 Å |
| N-terminal A1–4 | Protenix | 1.243 Å | 0.259 Å | 0.000 Å |
| Pair A10–13 / A46–49 | RFD3 | 1.340 Å | 0.406 Å | 0.306 Å |
| Pair A10–13 / A46–49 | Boltz | 1.100 Å | 0.410 Å | 0.662 Å |
| Pair A10–13 / A46–49 | Protenix | 1.020 Å | 0.393 Å | 0.600 Å |

Boltz native observed maxima are 4.450/4.373 Å before restoration; the other two
models preserve observations to numerical precision in these pilots. These are
single-seed case observations, not a model ranking or terminal-gap qualification.
Independent evidence: `expansion-seed7-independent-readback/` under the new root.
The separate attempt02 readback also records the existing source-only flank graph:
the two gaps have a 2.840 Å minimum contact between their observed eight-residue
flank contexts and form one component at the unchanged 6 Å rule. This establishes
source-context proximity, not generated-loop entanglement or a benefit of joint
sampling. No generated/reference coordinates were used to choose that graph.

The unchanged private MODELLER protocol also completed both masks times six seeds
(12 new proposals). Source-restored raw passes are 5/6 for the terminal mask and
0/6 for the pair; final geometry/custody/L passes are **6/6 each**. All-outcome
heavy RMSD medians are 6.674/3.395 Å (backbone 5.524/1.510 Å). Matched **seed7**
MODELLER heavy errors are 5.923/2.376 Å; comparing six-seed medians with one ML
seed would still be unmatched. Evidence: `expansion-modeller-attempt01/` and
`expansion-modeller-independent-summary-attempt02/`. The first expansion summary
contained copied generic wording about 30 outcomes/seed7 recovery; the second
corrects metadata to 12 and no recovery. It does not change any scientific output.

The remaining five development seeds for the same two masks and three models
are now separately declared (30 new serial calls), retaining seed7 without
resampling. No inputs, masks, contexts, thresholds or per-backend sampler settings
are tuned from pilot RMSDs. Common repair remains necessary-only CPU, all generated
atoms together, no anchor motion. Evidence:
`expansion-development-extension-attempt01/`. All thirty new calls subsequently
completed, followed by necessary-only common CPU repair and independent readback.

### Completed six-seed terminal/pair distributions

The following includes original seed7 **without resampling** and the five new
seeds. Every cell contains exactly six outcomes. The 36 ML and 12 MODELLER final
structures all pass independent custody/geometry/global-L readback. ML raw-clean
count is 21/36; the other 15/36 required repair. Six distinct final outputs exist
per cell; RFD3 also has six distinct first-noisy arrays per mask, independently
checked against its consumer receipts. Native/restored outcomes remain separate.

| Mask | Backend | Raw restored geometry | Final passes | Median heavy RMSD | Heavy range | Median backbone RMSD |
|---|---|---:|---:|---:|---:|---:|
| N-terminal A1–4 | RFD3 | 5/6 | 6/6 | 1.164 Å | 1.024–1.265 Å | 0.303 Å |
| N-terminal A1–4 | Boltz | 2/6 | 6/6 | 1.291 Å | 1.005–1.356 Å | 0.585 Å |
| N-terminal A1–4 | Protenix | 6/6 | 6/6 | 0.882 Å | 0.869–1.243 Å | 0.236 Å |
| N-terminal A1–4 | MODELLER | 5/6 | 6/6 | 6.674 Å | 2.253–10.483 Å | 5.524 Å |
| Pair A10–13 / A46–49 | RFD3 | 4/6 | 6/6 | 1.202 Å | 1.036–1.422 Å | 0.326 Å |
| Pair A10–13 / A46–49 | Boltz | 0/6 | 6/6 | 1.067 Å | 0.964–1.140 Å | 0.373 Å |
| Pair A10–13 / A46–49 | Protenix | 4/6 | 6/6 | 1.206 Å | 1.020–1.709 Å | 0.365 Å |
| Pair A10–13 / A46–49 | MODELLER | 0/6 | 6/6 | 3.395 Å | 1.872–4.035 Å | 1.510 Å |

“Raw restored” is not native source-preserving acceptance: Boltz native observed
maxima reach 7.026/7.057 Å, versus numerical precision for RFD3/Protenix on these
masks. MODELLER remains the explicit private one-fast-loop comparator, not a
qualification or replacement of the shipped default pipeline.

Protenix is descriptively strongest on this terminal mask (lowest median and no
repair). Boltz has the lowest final heavy error on this pair but every proposal
needs restoration/repair; RFD3 has the lowest median pair backbone error. These
are different endpoints and do not identify a universal winner. Both tasks are
from one cropped protein; six seeds do not provide six independent biological
cases. Native topology, optimizer trajectories, training independence and broader
terminal support are not established. No heldout-based best-of-six selection or
threshold tuning was performed.

Independent extension evidence: `expansion-development-independent-readback/`.
Combined all-seed/per-candidate summary and checksums:
`expansion-six-seed-summary/{README.md,summary.json}` under the new evidence root.

## Source-only conditioning A/C/D

Thirty separately declared RFD3 calls yielded 24 assembled outcomes. All 24 passed
independent final inventory, exact observation, geometry and global-L readback.
Each arm uses the six development seeds. Baselines reuse the original samples;
two necessary CPU repair derivatives platform-match the full-heavy baseline without
replacing its historical Reference-platform outcomes.

| Arm | Raw restored passes | Final passes | Median heavy RMSD | Median backbone RMSD |
|---|---:|---:|---:|---:|
| Eight-residue flanks | 6/6 | 6/6 | 1.038 Å | 0.239 Å |
| Spatial context, 6 Å rule | 1/6 | 6/6 | 1.037 Å | 0.360 Å |
| Backbone-only observations | 3/6 | 6/6 | 1.031 Å | 0.312 Å |
| Independent pair proposals | 3/6 | 6/6 | 1.334 Å | 0.349 Å |
| Full-heavy single-gap baseline, CPU-matched | 4/6 | 6/6 | 1.019 Å | 0.273 Å |
| Joint pair baseline | 4/6 | 6/6 | 1.202 Å | 0.326 Å |

The spatial arm changes entity layout and chain cuts as well as observation count.
Independent pair calls omit the complementary unknown gap, changing task size;
their assembly preserves individual donor lineage and jointly refines both regions.
These are not isolated causal demonstrations of context or joint-sampling superiority.
Hidden reference coordinates were used only after proposal/refinement for scoring.
Evidence: `conditioning-acd-attempt01/`, `conditioning-acd-independent-readback/`
and `conditioning-full-heavy-common-cpu/` under the continuation root.

## Necessity-gated anchor calibration

The grid executed 36 CPU attempts on three remaining source-preserving junction
failures, with separate movement, metadata and observation-custody audits.
1A3N16 seed103 passed 0/12; 1A3N16 seed104 passed 6/12 (four unique endpoints);
1FR2_16 seed104 passed 0/12. Thus only one of three starting proposals was rescued.
The successful geometry endpoints still have approximately **10 Å** hidden-region
error: geometry repair is not reconstruction accuracy.

The strong-restraint direct A29 case moved backbone RMS 0.1346 Å, Cα RMS 0.0929 Å
and maximum heavy atom 0.2040 Å within the experimental 0.25 Å envelope. This is
case-specific bounded motion, not calibrated physical uncertainty or approval of
positive production limits. Validation therefore freezes **zero observed motion**.
Evidence: `anchor-calibration-attempt01/` and `anchor-calibration-independent-readback/`.

## Interface admission and topology limitations

All three author-declared 1BRS dimers were assessed against frozen spans/mask:
A/D refused partial occupancy at D14; B/E lacked the required interval; C/F lacked
heavy atoms at C19. Admission is 0/3, with no sampling, trimming, repair or occupancy
rewrite to create a nominal success. The original 1A3N interface refusal remains.
Evidence: `1brs-deposited-assemblies-admission-attempt01/`.

New no-sampling controls retain an analytically legitimate linked polygon pair,
rigid-transform invariance, clean endpoints with a synthetic intermediate surface
traversal, and shared SSBOND parsing preserving case-sensitive chains/insertion codes.
They pass but do not qualify protein knots, disulfide surfaces or unsaved optimizer
trajectories. Surface hits remain nonfatal diagnostics. Evidence:
`topology-limits-controls/`; tests: `tests/test_mac_topology_limits.py`.

## Independent validation and sequence evidence E

The [fixed-source validation policy](../plans/mac-reconstruction-fixed-validation.md)
admits six masks on original 1PGA A1–56 and 2CI2 I19–83 chains, ordinal intervals
`[9,13)`, `[23,31)` and `[0,4)` each. Original source entries and provenance are
tracked under `tests/fixtures/protein_context/`. Genuine terminal observed OXT
remains byte-identical; RFD3/Boltz explicitly pass it through rather than claim a
native constraint on an atom absent from their canonical axes.

After all 18 weightless model/case preflights passed, attempt03 launched the frozen
90 serial ML plus 30 MODELLER calls with seeds 211/223/227/229/233. Policy hash and
the two pre-execution-only failures are documented in the linked plan. Execution
and reserved-seed use require receipts, not inferred completion. Independent
readback subsequently completed for **all 120 outcomes**, checking policy/input
and native-output hashes, native drift, exact observations/OXT, inventory,
geometry, global chirality and heavy/backbone RMSD.

| Backend | Completed | Raw restored geometry | Repair outputs | Final accepted | Global L |
|---|---:|---:|---:|---:|---:|
| RFD3 | 30 | 26 | 4 | 30 | 30 |
| Boltz | 30 | 7 | 23 | 30 | 30 |
| Protenix | 30 | 22 | 8 | 30 | 30 |
| MODELLER | 30 | 1 | 28 | 23 | 30 |

All 120 preserve source records, including one observed OXT each. Boltz/MODELLER
native observed maxima after observation-only fit reach 12.625/10.714 Å; this is
not native fixing. RFD3/Protenix represented observations agree to numerical
precision, with the previously stated different OXT transport semantics.

All-outcome heavy RMSD medians in Å, five declared seeds per cell, no selection:

| Mask | RFD3 | Boltz | Protenix | MODELLER |
|---|---:|---:|---:|---:|
| 1PGA internal4 | 1.087 | 1.247 | 1.465 | 2.657 |
| 1PGA internal8 | 1.473 | 1.435 | 1.025 | 7.923 |
| 1PGA Nterminal4 | 1.105 | 0.814 | 1.017 | 5.010 |
| 2CI2 internal4 | 0.765 | 0.627 | 0.451 | 2.348 |
| 2CI2 internal8 | 1.619 | 0.851 | 1.266 | 7.484 |
| 2CI2 Nterminal4 | 8.308 | 3.036 | 3.663 | 6.253 |

Seven MODELLER rejections remain: all five 2CI2 internal8 seeds, 2CI2 internal4
seed229 and 1PGA Nterminal4 seed233. Seed229 internal4 retains a numerical repair
failure, raw fallback and traceback; the other six retain failed final geometry.
No replacement/resampling occurred. Rejections remain in descriptive RMSDs, not
usable-accuracy claims. RFD3's 2CI2 terminal cell passes geometry 5/5 despite median
heavy error 8.308 Å: geometry does not imply reconstruction accuracy.

Full per-seed heavy/backbone errors, native drift, failures and hashes:
`fixed-source-validation-independent-readback/{README.md,independent-results.json}`.
Two entries are not 120 independent proteins; training-set absence is unestablished.

Series E launched after validation completed: 12 Boltz/Protenix development
calls remove structural templates with no MSA in either arm. Read-only consumer
guards require actual zero template masks; Protenix's source-reinjection callback
remains unchanged. Existing source-template arms require CPU-matched repair before
comparison. No frozen local homolog/MSA assets exist for the other E cells; remote
sequence upload remains unauthorized. ML jobs do not overlap. Evidence:
`sequence-evidence-E-attempt01/`. All 12 new calls subsequently completed, as did
the 12 retained-proposal CPU matching jobs in `sequence-evidence-template-CPU-attempt03/`.
No template-arm resampling or historical Reference output replacement occurred.
Attempt01 retains an unused-import-only pre-execution amendment. Attempt02 stopped
at its first worker's duplicate directory reservation, before refinement; the
corrected worker explicitly requires the fresh executor-owned reservation and
refuses existing scientific payloads. Both failures and original scripts remain.

Independent readback passed all **24/24** final custody/inventory/geometry/global-L
outcomes. Actual consumed no-template masks are zero in all 12 new calls.

| Backend | Structural templates | Raw passes | Final passes | Repairs | Median heavy RMSD | Median backbone RMSD |
|---|---|---:|---:|---:|---:|---:|
| Boltz | None | 0/6 | 6/6 | 6 | 1.030 Å | 0.399 Å |
| Boltz | Source template | 0/6 | 6/6 | 6 | 1.035 Å | 0.393 Å |
| Protenix | None | 5/6 | 6/6 | 1 | 1.320 Å | 0.369 Å |
| Protenix | Source template | 4/6 | 6/6 | 2 | 1.041 Å | 0.313 Å |

Both arms have no MSA. Protenix still receives observations through its unchanged
reinjection callback; “no-template” does not mean sequence-only Protenix. Boltz
receives no source coordinates in its new sequence-only proposal arm, then uses
the unchanged observed-only fit/restoration. This retrospective single 1UBQ mask
does not establish isolated universal conditioning effects or significance of the
small Boltz difference. No reference-based selection was performed. Evidence:
`sequence-evidence-independent-readback/{README.md,independent-results.json}`.

## Remaining gates

The declared E subset and frozen validation execution/readback are complete.
Clean interface execution, broad independent-family
accuracy, calibrated positive anchor limits and native topology/trajectory
qualification remain unestablished. Whole-1HE8 acceptance and visual approval of
the new joint RFD3 donor ensemble remain separate unresolved gates; no limits are
relaxed to accept them. No public CLI/Homology/GUI capability was enabled, no P3
implemented, and no commit/push performed. Further sampling requires a separately
bounded cohort/asset plan; manually approving the new 1HE8 donors and any public
integration remain user decisions, not automatic consequences of these studies.

## Engineering verification at this checkpoint

The focused comparison/mask/transport/resource/refinement and portfolio/diagnostic
regression set passed **61 tests**, including independent readback mutation,
summary cohort-size/recovery metadata, retained failures, incomplete/repeated seeds
and unfiltered descriptive RMSD. Ruff checks of the new comparison scripts and
tests passed, agent documentation maps passed, and generated CLI reference is in
sync. Subsequently the full fast Python lane passed **2460 passed, 3 skipped,
33 slow deselected** (`fast-suite-attempt01.log`). The three new topology-limit
tests and two native Protenix MPS hook controls passed separately. Full slow and
GUI suites have not been rerun; no GUI code changed.

The first final fast-suite rerun exposed failures and exceeded its 240-second
outer limit. A retained fail-fast diagnosis (`final-fast-suite-attempt02.log`)
identified absent `dvbfixer` on the child process PATH (350 passed before the
first failure). A new full run explicitly prepends the selected environment's
`bin` directory and retains `final-fast-suite-attempt03.log`: **2472 passed,
3 optional Torch/Biotite skips, 33 slow deselected**, 291.68 seconds. No broad
skip, source modification or environment mutation was needed. The interrupted
first run and failed PATH diagnosis remain retained. Full slow and GUI suites
were not rerun. New auditor, read-only
consumer guard and baseline-directory regression checks remain distinct from
scientific execution/readback evidence.
