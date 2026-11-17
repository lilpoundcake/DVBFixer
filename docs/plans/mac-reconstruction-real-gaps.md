# Reconstruction continuation handoff: real gaps, joint repair and remaining gates

**Current entry point (2026-11-14):**
[Experimental structure coordinate completion](experimental-structure-completion-handoff.md).
Use it for current scope, status and continuation instructions. This document
retains historical evidence and proposals; contradictory older next-action or
pending-decision statements do not override the current handoff.

## Current user direction: existing list, consistent restoration validation

The user's continuation after the 1KTZ summary removes a separate structure-
admission phase as a prerequisite. **Work with the existing retained list**;
do not turn the 177 pending review labels into another gate before restoration.
Existing exclusions/deferred submodels remain unchanged, with no new recruitment
or model-size calibration. Review metadata remains provenance, not a new queue.

Receptor terminal chemistry is deferred to future user-facing preparation, not
a current research blocker. Preserve both OXT diagnostics and do not silently add
atoms, caps or charges. The candidate is not MD-ready, but its original 710-atom
restoration task need not wait for that separate choice.

The [read-only validation consistency reproduction](../../tmp/mac-real-gap-continuation-2026-11-14/validation-consistency-review-attempt01/report.json)
applies the same full chemistry, Ramachandran, chi12 and deterministic protein
steric diagnostics to source, first narrow-pass export and current export. It
reproduces the missed THR β warning in the former and zero new geometry findings
in the latter. The narrow validator omits stereocentre/full S–S checks and uses
severe-only sterics; its pass must not be called a clean full-diagnostic result.
Report custody, restored-mask completeness, new geometry, unchanged source
findings and later preparation separately. No thresholds, coordinates, frozen
requests or historical reports were changed by this reproduction.

<a id="current-2026-11-14-outcome-terminal-chemistry-pending"></a>

## 2026-11-14 retained outcome: terminal chemistry deferred

The full four-chain 1KTZ candidate now contains all **710 restored heavy atoms**
in the frozen mask, 468 protein residues and all 326 original waters. Independent
[full-context readback](../../tmp/mac-real-gap-continuation-2026-11-14/1KTZ-restored-heavy-export-attempt02/readback-receipt.json)
passes custody: 2,986 original protein records and 326 water records unchanged,
454 L-Cα centres, reviewed THR C60 β-R and 21 target S–S pairs. Full CLI reports
**2 ERROR / 8 WARNING / 2 INFO**: no new warnings, but genuine missing terminal
OXT at receptor ASP B136/D136. These two atoms are outside the original frozen
710-atom mask, are not baseline errors and have not been published.

**Future preparation choice, not a current blocker:** free, uncapped B/D C-termini,
permitting a separately declared two-OXT completion, or specified amidation/capping. Native oxidized
TGF-β3 and its 21-pair target were already confirmed; do not ask that again.
Until terminal chemistry is confirmed this candidate remains incomplete and
not MD-ready. No new H, protonation, charge, cap, public profile or frozen-request
change is authorized by the geometry result.

See the [candidate/provenance README](../../tmp/mac-real-gap-continuation-2026-11-14/1KTZ-restored-heavy-export-attempt02/README.md)
and [dated continuation](../research/mac-real-gap-continuation-2026-11-12.md#2026-11-14-autonomous-terminal-geometry-continuation).
The current projection is `f6360681cd9d0e9c98367b0dbd8e5604fd69ae75114d45950bc025a7c1eb9b13`;
full-water PDB is `d19713a3cb0034c05d54dc4ef3e0d695ff0980289a47674888793dcf27c6bfef`.
The final fast suite passes 2,848 tests (three optional skips, 33 slow deselected,
two SWIG warnings, 299.81 s). This is one retained case, not broad qualification,
experimental missing-coordinate accuracy, a knot/trajectory guarantee or cohort
completion. The 177 pending entries, four deferred submodels, 12 controls and
excluded 7N8X remain unchanged; recruitment and size calibration stay paused.

## Historical continuation: approved native terminal study

**Current state (2026-11-14): target native oxidized TGF-beta3 confirmed by the
user; no redox decision remained pending at that stage.** The separate21-pair native
target is approved, without changing the18 original explicit source annotations,
frozen request or old native graph failure. Native-topology preflight and its
[independent readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-topology-preflight-attempt01/readback-receipt.json)
pass: actual21pairs/42singly partnered CYS,exact202active A/C1–12atoms,3,494fixed
nonactive protein atoms and326fixed water oxygens. Native residual is7.97e-6Å;
written fixed-coordinate residual is zero. This is not coordinate restoration.

One separately declared deterministic MODELLER stereochemical terminal procedure
then executes250/750/1000CG-iteration caps with only those202atoms active, native
chirality restraints and all21DISU patches. No neural sample, randomization or
unrestricted initial automodel optimization occurs. Every native stage has zero
nonactive motion, followed by native global-L and fidelity checks **before**
materialization. Its [independent readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-terminal-restoration-attempt01/readback-receipt.json)
confirms2,986byte-identical source observations,3,494unchanged inactive protein
records,326fixed waters and454L centres. The newly declared60Å terminal scope is
not an increase of the failed3Å local-cleanup declaration; maximum active motion
is16.117Å. A/C SG7–SG16 distances improve to2.064/2.321Å, but terminal geometry
still fails. All-generated water contacts fall from13ERROR/10WARNING to8ERROR/
7WARNING; immutable baseline remains0ERROR/8WARNING. No accepted restoration.

A separate one-thread/no-solvent CPU FF derivative then uses the **actual approved
21-pair graph**, moving only710missing-coordinate heavy atoms and freezing all
2,986observations/326waters. Its3Å movement bound is explicitly relative to the
terminal parent, not original RFD3 coordinates. The
[independent FF readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-context-refinement-attempt01/readback-receipt.json)
confirms exact observation custody,global L,all21pairs and0ERROR/3WARNING restored
water contacts. Nevertheless **reject**: fifteen atoms exceed3Å (maximum5.237Å,
B24CZ), while C3–C4 amide planarity and11severe protein overlaps remain. A7–A16
twist and C7–C16 length/twist warnings also remain. Do not silently raise the bound,
adopt this derivative as accepted or erase its lineage. Source A77–C77 warnings
are unchanged. Next research is terminal conformational escape/geometry repair
and a prospectively enforced movement scope, not another blind minimization or
an unresolved redox question. At that historical stage no new user decision had
been established; the current receptor-terminal decision is described above.
The expanded focused set passes95tests; full fast suite passes2,791, with three
optional skips,33slow deselected and two SWIG warnings in299.77s.

### Historical decision boundary, resolved by the user's confirmation

Renewed autonomous continuation added private derivative
custody/fixed-water controls and declared one generated-only FF refinement.
That attempt stopped during FF preparation, before missing-heavy coordinate
minimization: its exact18-pair guard was too narrow for OpenMM's source-only
19-pair topology. No derivative PDB or new model sample was produced. The extra
A77–C77 pair already exists by distance in the immutable observed source; do not
delete it to satisfy the guard or mislabel it as a model-created edge.

The [source/biochemical audit](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-disulfide-semantics-attempt01/readback-receipt.json)
confirms original ASU9explicit S–S records render to18 in this assembly; no77pair
annotation was lost by that renderer. One retained-source DBREF-directed UniProt
P10600 metadata query, not a new structure, confirms exact mature-sequence matches
for A/C and experimentally supported native disulfide annotations. They support
an explicit **21-pair native hypothesis**, including A77–C77 and the two absent
CYS7–CYS16 links. This does not establish the redox state of this particular
crystal. A77–C77 has two immutable source geometry warnings (2.233Å length and
153° twist), which cannot be repaired while preserving all observations.

The two raw CYS7–CYS16 SG distances are9.263/17.549Å. Even reaching a generous
2.5Å envelope requires at least6.763/15.049Å motion of the unobserved endpoint,
with its observed partner fixed. Thus the prospectively bounded3Å local cleanup
cannot restore that native topology. Do not retry it, raise its bound or add the
three hypothesized links silently. The genuine decision is whether the target
is native oxidized TGF-beta3, permitting a separately planned disulfide-aware
reconstruction of the missing A/C termini, rather than assuming unresolved source
annotations define a free-CYS molecule. Original atoms/waters must remain fixed.

The [private derivative helper](../../src/dvbfixer/model/diffusion/rfd3_refinement.py)
keeps native raw fidelity separate from intentional missing-atom-only motion and
creates FF-only water serial maps without rewriting original context. Fifteen
new tests pass, including actual frozen water-O/disposable-H kernel behavior;
the focused set passes75 and the full fast suite passes2,771 (three optional
skips,33slow deselected,two SWIG warnings,297.22s). This does not close the target
chemistry decision or change the rejected raw candidate/public profile.

**Latest result: complete 1KTZ executes on Mac/MPS within the original budget;
the resulting candidate is rejected scientifically.** The exact-pair/cache/
Transition-row arm completed in443.09s (7min23s), sampled peak MPS driver
5,152,391,168bytes (4.80GiB), peak aggregate RSS2,906,308,608bytes, without CPU
fallback or a guard stop. The 200-point configured schedule made199 native
diffusion-module calls with two recycles. This is one complete-target engineering
smoke, not general Mac capacity or terminal-gap qualification.

The [independent final readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-mac-transition-attempt01/readback-receipt.json)
passes actual identity/masks and raw all-observation common-frame fidelity
(maximum7.96e-6Å), then confirms all2,986 observed records byte-for-byte,710
restored atoms and18 unchanged S–S annotation records. Canonical declared-link
metrics emit no findings and global L-chirality passes. However, generated/junction
geometry fails (including bond length, steric and four peptide-junction reports),
and the independent all-water cross-set readback confirms13ERROR and10WARNING
contacts involving restored atoms. The original8source/water warnings remain
separate. The candidate remains `candidate-UNQUALIFIED.pdb`, not a public output
or accepted structure. No refinement, source/water change or further sampling
was performed. Mac execution and acceptable restoration quality remain distinct
outcomes; further quality work needs its own bounded scope.

**User correction: Mac execution is the objective, not transfer to Linux/CUDA.**
The resource-stop/host-decision wording below records the first attempt and is
superseded as a next-action recommendation. Continue explicit Mac memory-path
engineering without cropping the complete target or increasing limits.

Private exact-expression query-row tiling is implemented in
[`mps_pairwise.py`](../../deploy/rfd3/mps_pairwise.py) and
[`mps_pairwise_extended.py`](../../deploy/rfd3/mps_pairwise_extended.py).
It preserves full dense pairs, native projections/MLP layers and masks; it does
not enable native `low_memory_mode`, whose initializer path bypasses token
feedback and the optional atom transformer. Fifty-six isolated CPU/MPS controls
cover expressions and the [explicit cache hooks](../../deploy/rfd3/mps_cache_boundary.py),
including unchanged outputs/inputs/parameter keys and CPU/MPS RNG states.
Single fixed-axis operator controls are not new structures or model-size sweeps.

The first sinusoidal-only Mac arm completed all pair embeddings but stopped near
6.09GiB. The complete-pair arm completed the native initializer and reached the
first denoiser encoder, then stopped at7.13GiB; no final coordinates exist in
either arm. These are technical arms of the same attempted seed101, not
independent scientific seeds. All original failures remain retained. The
[cache-boundary arm](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-mac-cache-attempt01/declaration.json)
adds synchronization/release of unused MPS allocator cache at selected native
module boundaries, without deleting live tensors or changing watermarks,
weights, attention, precision or model settings. Its actual outcome must be
read separately; operator tests alone do not establish a working full Mac run.

That cache-boundary arm also stopped at6.15GiB after three completed native
diffusion-module calls; the fourth call was partial. This is not a completed
scientific seed or established timestep count. Retained boundary telemetry locates
the next peak at the token encoder Pairformer gated `Transition`, whose expanded
pointwise temporaries still covered all query rows. The
[Transition-row arm](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-mac-transition-attempt01/declaration.json)
therefore tiles this native pointwise expression too, preserving the entire last
feature dimension, RMSNorm, gating and all projections. Dense-vs-tiled native
Transition controls pass on CPU/MPS; attention is not tiled or approximated.
The new bounded arm was launched separately; its measured full-model result is
recorded above, not inferred from the 56 operator/cache controls. No engine source file, checkpoint,
source input, complete target or safety limit has changed.

**Historical first attempt: resource-blocked; no completed coordinates.** The one declared
seed101 attempt stopped after 19.20s for `mps-driver-limit`. Sampled peak MPS
driver allocation reached 14.35GiB against the unchanged 6GiB stop threshold
(periodic monitoring can overshoot; this was not an allocation cap). The child
exit is confirmed, all partial evidence retained, and no retry/refinement or
limit increase occurred. The [independent readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-run-attempt01/readback-receipt.json)
confirms actual native identity/masks and the input common frame, not final
coordinate fidelity: no final native output or candidate PDB exists. Fifty-eight
MPS-bearing submodule calls are not completed samples or diffusion steps.
The 200 timesteps remain a configuration, not an executed-step claim. Further
full-target execution needs a separate resource/host decision; no silent crop or
new structure is authorized. Final geometry/water acceptance is unassessed.

Latest software checks: 75 focused tests and the full fast suite pass
**2,756 tests**, with three optional skips, 33 slow deselected and two SWIG
warnings in 300.92s. This is software evidence, not structural qualification.

The user approved the separate bounded full-1KTZ study after the offline custody
work. This supersedes earlier pending-authorization statements below, not the
historical profile refusal or failed native 18-S–S-edge transport verdict.
The [study declaration](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-run-attempt01/declaration.json)
retains the complete original four-chain target, eight gaps and 468 residues.
One prospective seed101 is declared, with the existing MPS seed amendment and
unchanged serial pressure-v2 limits; no automatic retry, refinement or expansion.

The [source/water preflight](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-preflight-attempt01/readback-receipt.json)
passes with 2,986 unchanged observations, observed global L-chirality, 18 declared
link checks and 326 original oxygen waters. Eight baseline protein–water warnings
and zero errors remain visible. The worker's lowercase severity counter was
wrong; detailed findings and independent cross-set recomputation agree on eight
warnings. Original worker/result remain retained. No source or water was repaired.

Native conditioning remains protein-only (2,984 fixed atoms); two observed OXT
are preserved-only and all waters are immutable final-assessment obstacles, not
native conditioning. Actual native identity, masks, reference coordinates and
float32 final tensors are captured losslessly before CIF/PDB rounding. One
all-observation proper common-frame fit must meet 1e-4Å maximum residual before
generated-only materialization. Public profiles/frozen requests are unchanged.
See [the custody plan](rfd3-source-link-custody.md) for P3/P4 and evidence limits.

Date: 2026-11-12. User direction: small artificially hidden reference panels are
insufficient; prioritize real incomplete structures and make difficult cases such
as 1HE8 work. This records the next scientific scope, not completed new sampling.

Status: **design for the next continuation**, incorporating the unfinished work
from the original handoff. The user asks to retain very difficult real cases and
avoid an arbitrary maximum gap/tail length. This document does not authorize an
unbounded sampling campaign, a new public capability or a silent policy change.

## Current user amendment: loop completion, not domain-scale prediction

The user subsequently narrowed the current evaluation: **do not use experimental
entries with enormous missing regions such as7N8X** when completion becomes
prediction of substantial missing domains rather than local loop reconstruction.
This is a scientific task-scope decision, not a restriction on public structural
data or an assertion about biological risk.

Terminology clarification: the current task is **restoration of missing coordinates
at known sequences**, not designing new amino-acid sequences. Use coordinate
restoration/completion in new explanations. Frozen technical fields such as
`generated_atoms` retain their compatibility spelling and denote restored-atom
ownership; changing wording does not change scientific scope or qualification.

This amendment takes precedence over the historical unrestricted-length coverage
and compulsory100+ quotas below.7N8X is deferred from the active loop benchmark;
analogous domain-scale or focused-submodel omissions are deferred as well. Retain
their original bytes, requests, preparation outcomes and provenance as a separate
research reserve. Do not count them as admitted loop cases or launch new proposals,
donor selection, refinement or evaluation for them under the current plan.

Admission now requires an explicit loop-scale interpretation in a sufficiently
observed scaffold. Gap length alone neither establishes a loop nor justifies
automatic admission; no new universal numeric cutoff is asserted. Unclassified
large spans remain pending/deferred rather than implicitly allowed. Multiple
shorter gaps do not authorize silently ignoring a deferred domain-scale gap in the
same whole-entry request. Any future narrower task needs its own explicit mask,
context and provenance, without mutating existing/frozen requests.

The200-entry ambition and group-safe80/40/80 design remain planning targets, **not
authority to fill slots with deferred huge spans**. The original30/50/50/40/30
length quota table and mandatory long-tail coverage are historical, not current
execution requirements. A replacement loop-focused quota table is not yet frozen;
report the shortfall instead of inventing new counts or a length threshold.
Continue preserving failures and distinguishing input readiness from admission.
Existing bounded source-preparation runs are inventory evidence only: no new neural
sampling has begun, and their outputs are not grandfathered into the loop cohort.

The user additionally paused recruitment: **work with the retained structures for
now**. Do not declare or start additional metadata/source acquisition to fill old
quotas. Previously retained originals, completed outcomes and existing bounded
preparation inventories remain available; preserve shortfalls and pending scope
reviews rather than silently enrolling domain-scale cases. Continue the existing
plan with durable state updates, not a session-ending handoff unless requested.

Backend-specific input budgets are a separate proposed calibration, not a shared
model ceiling. The user subsequently deferred that work: **do other correctness
tests first; revisit limits later**. Keep the production512-residue contract and
frozen requests unchanged for now. Do not launch larger-input or checkpoint size
probes under the correctness-test work below.

## Completion tracker

### Authoritative current state — retained-data consolidation

All declared calculations are settled; **no new restoration runs are queued**.
The current task is missing-coordinate restoration at fixed known sequences in
retained experimental structures. Source recruitment and model-size calibration
remain paused. The two-case series is not being expanded.

The [verified execution ledger](../../tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/ledger-attempt01/execution-ledger.json)
and [case-arm outcomes](../../tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/ledger-attempt01/case-arm-outcomes.tsv)
bind 1253 verified payloads. All 233 source tasks are accounted for: 194 ready,
31 refused, eight time-blocked, zero unexecuted. Read-only review of all 194 inputs
and 1280 declared gaps initially preserved 181 pending reviews, 12 small controls
and the explicit 7N8X exclusion. The newer case-specific scope overlay below
defers four focused-submodel whole-entry requests, leaving177 pending reviews.
Readiness is not scientific admission: **zero final test
admissions and no frozen split**. The 200-entry and 80/40/80 figures remain planning
targets, not authority to fill quotas with unsuitable cases.

Recent results are two structures × three arms, not six independent structures.
All six source-restored raw endpoints fail geometry. One final case-arm passes:
MODELLER seed42 2BE1 after common repair (optimizer seed11). Seed101 fixation masks
were verified, but observed drift remains and its repaired endpoints pass 0/2.
Protenix retains four full-native resource refusals, including two saved raw
returns. No backend-superiority claim follows from unmatched contexts/seeds.
The private library has 18 stage records and zero manual donor approvals.
Latest software verification passes2745 fast tests (three optional skips,
33 slow deselected), with64 focused custody/materialization/metadata controls and the
documentation checks. The earlier consolidation receipt retains its2573-test
result. These are software checks, not completion of the scientific work packages.

| Work package | Settled evidence | Remaining gate; not marked complete |
|---|---|---|
| WP0: reproducibility | Original publication retained; payloads checked; latest2745-test fast suite passes | Continuation publication/runtime handoff remains separate from local verification |
| WP1: admission | All 194 deposited gap axes and original projection custody checked; context omissions counted exactly; 106 entries have explicit artifact-overlapping gaps; four family flags and 26 known exposures remain | Decide task semantics for construct artifacts and terminals; loop/whole-context admission, complete exposure accounting and group-safe split |
| WP2: joint repair | Original 1HE8 0/8 and buffered 2/4 joint passes under unchanged bounds | Causal force/defect explanation and a qualified general policy; passes alone are not causality |
| WP3: restoration/comparison | Six recent case-arms, six CPU jobs, 12 stages, one final pass | Broader admitted, matched-context task matrix and complete-context acceptance |
| WP4: conditioning | Native masks/20 peptide junctions checked; MODELLER drift localized in saved stages | Matched conditioning/MSA/comparator cells; native observation fixation remains unqualified |
| WP5: topology | All12 retained 1HE8 arms/467 instances recomputed; zero segment alarms. Focused TYR972 ring:455 assessable states without hits; historical controls reproduce | Other rings/closure, independent detector qualification and unsaved-path evidence; no fatal topology guarantee |
| WP6: review/delivery | Digest-checked private library, filters, overlays, exact generated-only sidecars | Manual donor approval, selection/conflict workflow and operational delivery; no public integration |
| WP7: evaluation | No final cohort, policy or split frozen | Freeze policy/budgets and untouched group-safe evaluation before execution |

Immediate order: advance case-level admission reviews using the verified all-input
matrix and worksheet; retained-data/software consolidation is now verified. Do not resume
the two-case experiment chain, infer an untouched final component from incomplete
exposure data, widen bounds or silently repeat consumed seeds. See the
[settled chronology](../research/mac-real-gap-continuation-2026-11-12.md#settled-seed101-results-and-all-input-consolidation).
The [independently checked review worksheet](../../tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/review-handoff-attempt01/adjudication-worksheet.tsv)
now carries all 194 source records with bound evidence paths; its
[receipt](../../tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/review-handoff-attempt01/receipt.json)
verifies sequences, gap axes and TSV consistency. Controls/exclusion remain visible;
all admission and split decisions are unset. Group ordering is not a scientific
ranking, and known exposure flags do not establish complete leakage accounting.

The newer [per-case review table](../../tmp/mac-real-gap-continuation-2026-11-12/retained-case-admission-reviews-attempt02/case-review.tsv)
advances the blank worksheet with concrete evidence, without admitting any case.
All 1280 gap axes match deposited observations for every compatible source-label
candidate. Actual full/projection comparisons confirm 48 entries with no omitted
context, 34 with only HOH/WAT omitted and 112 with other components; no hidden ATOM
loss or original projection coordinate changes were found. Explicit deposited
artifact annotations overlap 229 gaps in 106 entries (33 complete gaps in 23 entries).
Fifty-six database-relative deletion annotations lack experimental construct
positions and are not missing-coordinate targets. No mask is silently trimmed,
and absence of artifact annotations does not establish native-sequence identity.
The [bound review data](../../tmp/mac-real-gap-continuation-2026-11-12/retained-case-admission-reviews-attempt02/case-reviews.json)
separates these verified factual subchecks from scientific decisions. Sequence
semantics, context treatment and final admission/split remain open.

The [gap-level review table](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-reviews-attempt01/gap-review.tsv)
now separates 469 internal gaps from 419 N-terminal and 392 C-terminal gaps.
All 229 explicitly artifact-overlapping gaps are terminal, including 33 fully
annotated gaps. No mapped artifact annotation on an internal gap proves neither
native-sequence identity nor local-loop suitability. Exact observed sequence-axis
runs on both sides are reported without a new length or occupancy threshold.
Of the 181 active-scope pending entries, ten have only internal gaps, 101 mix
internal and terminal gaps, and 70 have only terminal gaps. These are review
categories, not a selected cohort or authorization to discard terminal masks.
The [independent readback](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-reviews-attempt01/readback-receipt.json)
checks 973 bound payloads, all 1280 rows/TSV cells, request axes, scaffold runs and
artifact coverage. Whole-entry context/construct questions remain in force.
The subsequent complete fast suite passes 2610 tests, with three optional
Torch/Biotite skips and 33 slow tests deselected; this is software verification,
not scientific admission or completion of the remaining work packages.
Next: resolve task semantics at case level using these separated facts; do not
silently turn mixed entries into internal-only requests or equate a two-sided
observed boundary with a scientifically admitted loop.

The subsequent [context-connection audit](../../tmp/mac-real-gap-continuation-2026-11-12/retained-context-connections-attempt02/connection-audit.json)
checks all 194 original full assemblies/projections and deposited connection
annotations. All 315 resolved explicit pairs between retained atoms (54 entries)
survive projection, with zero added pairs or resolved projection-boundary pairs.
Three source metal-coordination annotations in 3U8X/3UQ5 belong outside assembly1;
they are not projection losses. Source polymer membership, not sequence-compatible
chain candidates, determines protein/nonprotein endpoint classification.
The [readback](../../tmp/mac-real-gap-continuation-2026-11-12/retained-context-connections-attempt02/readback-receipt.json)
also retains a non-resolving SSBOND in small-control 5XZ0 to a crystallographic
image of the same residue. This remains an explicit symmetry-context review, not a deleted self-bond,
an automatic refusal or evidence that crystal context is represented. Absent
annotations, nonbonded context, operator transport and force-field coverage remain
unqualified; the 112 nonwater-context reviews are not cleared by these results.
The subsequent full fast suite passes 2627 tests, with three optional skips and
33 slow tests deselected. Focused connection/context/gap/admission checks pass
49 tests; these remain software checks, not scientific admission.

**Counts are not lengths or a single structure.** The1280 gap instances are summed
over194 entries and include repeated chain copies in biological assemblies. The
[length report](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/length-report.json)
and [per-entry table](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/entry-gap-lengths.tsv)
now keep gap instances, absent residue instances and experimental entries separate.
The876 within-entry exact-sequence/interval classes are descriptive equivalences,
not independent structures, identical conditioning or collapsed reconstruction tasks.

The [versioned case reviews](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/case-scope-reviews.json)
defer5MSC (A-domain study),7QL0 (filament-fragment study),8YNM (death-effector-domain
assembly) and8ZJK (core assembly) under the existing focused-submodel amendment.
Deposited study scope and exact whole-entry missing spans, not a new length cap
or generated outcomes, support these individual decisions. Original sequences,
masks and earlier tables remain unchanged; no internal-only substitute was created.
The [effective scope overlay](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/effective-scope-review.tsv)
therefore contains177 pending, four deferred focused submodels,12 controls and
one explicit user exclusion. Earlier181-pending counts describe the retained
inventory before this overlay. None of the177 is scientifically admitted.

Within those177 pending entries, internal gaps have median9 residues and range1–141
(450 gap instances); N-terminal gaps have median12/range1–672 (393 instances),
and C-terminal gaps median12/range1–1069 (363 instances). Large unclassified spans
are still pending/deferred, never silently launched as loop tasks. The existing
7N8X exclusion, including its812-residue internal span, is reported separately.
The [readback](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/readback-receipt.json)
recomputes distributions and verifies all194 entries/1280 gap instances and399
bound inputs. It verifies decision custody, not a second independent scientific
adjudication. Continue remaining case reviews without shrinking masks, recruiting
new sources, calibrating input limits or restarting the settled two-case series.
The new complete fast suite passes2648 tests (three optional skips,33 slow
deselected);44 focused length/gap tests pass. Software checks do not admit any case.

The [remaining ten-arm segment audit](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-remaining-segment-frames-attempt01/completed.json)
is now complete: all367 retained optimizer frames plus30 source/proposal/final
stage instances were checked using the actual serialized covalent force graph.
Combined with the previous two-arm audit, the
[summary](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-remaining-segment-frames-attempt01/combined-arm-summary.tsv)
covers all12 arms/467 instances (431 optimizer frames and36 stage instances), with
zero segment alarms. A [separate readback](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-remaining-segment-frames-attempt01/readback-receipt.json)
reconstructs each force graph and recomputes every screen without importing either
audit worker. It uses the same distance primitive, not a second independent
topology detector. The0.75Å diagnostic threshold, prior acceptance results and
movement bounds are unchanged. No coordinates were restored or optimized anew.
This closes the saved-frame coverage backlog, not WP5 as a whole: rings, closure
ambiguity and unsaved intervals remain unqualified. All ten earlier arm rejections
remain rejections; no public gate or case admission follows from zero alarms.
The full fast suite passes2662 tests (three optional skips,33 slow deselected),
and23 focused segment/continuation controls pass. Ring/closure and continuous-path
questions remain the next WP5 boundary; no new restoration task is queued.

The [focused ring readback](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-ring-frames-attempt01/readback-receipt.json)
now checks the previously implicated TYR972 cycle against explicit heavy-atom
covalent paths in all12 saved arms. Of467 stage/frame instances,455 have an
assessable ring and no piercings; the12 original source instances lack this ring
and are explicitly unassessed, not counted as negative findings. The fixed0.10Å
planarity scope was not tuned. Four separate historical controls reproduce two
known positive endpoints, one defined negative and one undefined nonplanar raw
surface. These are different retained artifacts, not evidence that an affected
endpoint was repaired by the twelve-arm experiment. The
[per-frame table](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-ring-frames-attempt01/ring-frame-findings.tsv)
keeps statuses, witness-bond lengths and topology provenance separate.
Retrospective controls have no restored-atom ownership mask: their raw ownership
strings are not provenance findings and are marked unqualified in the readback.
All results were independently recomputed using the shared canonical graph helpers
and ring primitive, not a second independent detector. Only this one ring is
covered; other rings, macrocycles, open-path closure and unsaved intervals remain
unqualified. No new coordinate restoration, optimization, source acquisition,
acceptance-policy change or public gate was introduced.
The full fast suite passes2671 tests (three optional skips,33 slow deselected),
with28 focused ring/continuation controls passing. Broader WP5 qualification and
the177 pending source reviews remain open.

Three [case-specific construct-target reviews](../../tmp/mac-real-gap-continuation-2026-11-12/construct-target-reviews-attempt01/construct-target-reviews.json)
now resolve one source-semantics question for1AHW,1JPS and1KTZ: every compatible
deposited source candidate has the exact original requested construct sequence and
heavy-coordinate absence mask. Their targets are already the deposited Fab/
extracellular constructs, not silently expanded full-native proteins or precursors.
Domain wording in a title is therefore not an automatic deferral rule: these
individual construct-target reviews differ from the earlier four focused-submodel
requests that extend substantially outside the observed core. No universal length
cutoff or new sequence alignment was introduced.

The [summary](../../tmp/mac-real-gap-continuation-2026-11-12/construct-target-reviews-attempt01/construct-target-summary.tsv)
covers three previously exposed development candidates and20 gap instances,
including all terminal intervals and chain copies. This is a supported **partial
task-semantics review**, not scientific/backend admission: all177 cases still await
complete review. Recorded preparation-profile terminal refusals remain visible,
and no new profile bypass or cross-backend admission test was performed. Original
water context in1JPS/1KTZ remains retained with endpoint treatment undecided;
1AHW has no omitted coordinate context. The
[readback](../../tmp/mac-real-gap-continuation-2026-11-12/construct-target-reviews-attempt01/readback-receipt.json)
independently recomputes exact source sequences/absence axes and preserves source-
label ambiguity. No mask trimming, native-sequence expansion, new smaller task,
source recruitment, structural-model call, optimizer call or final split follows.
The full fast suite passes2683 tests (three optional skips,33 slow deselected),
with23 focused construct/annotation/difference/admission tests passing. Source
semantics support does not override a recorded profile refusal or freeze a cohort.

A [weightless native RFD3 check](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-full-target-data-attempt02/native-stage-summary.tsv)
now executes the pinned input builder and pre-model pipeline on the entire1KTZ
target:468 residues/four chains/all eight gaps, with3694 canonical atoms. All2984
represented observations retain exact masks and relative coordinates within8e-6Å
float32 rounding after one common centering translation; all710 missing-coordinate
atoms remain zero placeholders and every sequence position stays fixed. Two
observed true-terminal OXT records remain explicitly preserved outside native
conditioning, not claimed as model constraints. No source was cropped or repaired.

The pipeline restores all464 peptide edges, including all10 observed/missing
junctions, but its atom-array graph lacks all18 declared disulfide edges. The
[independent readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-full-target-data-attempt02/readback-receipt.json)
confirms that both SG coordinates remain fixed and their source distances are
preserved. Missing adjacency is therefore **not physical rupture**, nor is this a
measurement of final neural features or model output. The declared strict graph-
transport check remains failed/unresolved; do not patch bonds by proximity, widen
checks, infer general RFD3 incapability or launch a model call to force a pass.
The downstream bond-feature contract and source-config feature prefix have now
been checked separately below; checkpoint-specific execution remains unqualified.
The first attempt's parent-only Gemmi import failure is retained; a lazy-import
correction reran identical conditioning bytes without changing the environment,
engine or input. All public/frozen profile refusals, water treatment and scientific
admission gates remain unchanged. The full fast suite passes2700 tests (three
optional skips,33 slow deselected), with21 focused RFD3/conditioning tests passing.

The [separate synthetic feature audit](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-token-bond-semantics-attempt01/readback-receipt.json)
confirms the installed `AddAF3TokenBondFeatures` intentionally suppresses every
edge between two non-atomized polymer tokens. Five two-node native controls and
independent saved-array readback pass: an explicit polymer–polymer edge and no edge
produce identical zero `token_bonds`; polymer–atomized and atomized–atomized edges
survive. This feature is **not the full covalent graph**. Merely inserting declared
protein S–S edges upstream would not make them visible in this channel when both
endpoints remain non-atomized. No real-input edge was inserted in this audit.
The strict1KTZ graph-transport result remains false, not waived or relabelled a
pass. Native atomization and coordinate-conditioning transport are checked by the
source-config prefix below; sampler behavior remains pending. No weights or model
were loaded.

A [native source-config prefix readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-feature-prefix-attempt01/readback-receipt.json)
now covers six saved stages through feature aggregation, stopping before diffusion
batching, noising and coordinate augmentation. Actual axes are468 non-atomized
tokens/3694 atoms, with zero added virtual atoms because all sequences are fixed.
The association scheme renames1368 atoms internally; `gt_atom_name` preserves
the original identities. Every2984 represented observation retains its mask and
coordinates in `motif_pos`; all710 missing-coordinate entries are masked to zero.
All468 one-hot sequence tokens decode correctly and the atom-to-token map matches
source identities. `ref_pos` is zero on this protein-only route: it is not the
observed-coordinate channel. Restricted `token_bonds` is468×468 and zero, while
all18 declared S–S graph edges remain absent. Six in-memory mutation controls are
refused by the independent readback. This uses resolved pinned **source defaults,
not authenticated checkpoint configuration**; it does not qualify sampling-time
fixation, output topology, MPS or scientific admission. OXT custody, source/water
decisions and existing profile refusals are unchanged. Next inspect fixation and
output-frame semantics separately without conflating feature transport with a
hard coordinate guarantee or overriding the strict graph result.

Autonomous continuation now separates three fixation layers. The
[eight synthetic driver controls](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-fixation-contract-attempt01/readback-receipt.json)
show that the numerical driver does not clamp arbitrary denoiser outputs; final
reinsertion in that isolated test is conditional. This is **not a failure of full
native EDM fixation**: native denoiser code sets fixed-atom time to zero, and its
EDM output scaling is identity for finite updates. Six direct formula controls
and six three-update analytic driver compositions across ordinary/padded shapes
are independently checked in the
[EDM receipt](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-edm-fixation-formula-attempt01/readback-receipt.json).
Both EDM compositions preserve synthetic fixed points bit-for-bit without
realignment. These are non-neural controls, not new structural sampling.

The [static metadata extraction](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-checkpoint-metadata-attempt01/selected-config.json)
identifies stored `train_cfg.model.net.diffusion_module.f_pred=edm` and
`inference_sampler.allow_realignment=False`. It reads only815458 bytes of bounded
ZIP `data.pkl` metadata and treats pickle callables/storage references as inert
symbols: no unpickling, tensor-member reads or model loading. This binds stored
configuration paths, not weights or arbitrary runtime overrides; absent sampler
fields are not invented. Nineteen static-reader controls pass. Three stored
transform fields differ from prior source defaults, so a separately declared
[stored-config prefix](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-checkpoint-prefix-attempt01/readback-receipt.json)
was executed before noising. All96 saved arrays match the prior prefix exactly;
the stored absent cache path is retained, not patched or populated. All1KTZ
identity/mask/sequence/motif findings and the18-edge strict failure are unchanged.

A [historical1UBQ readback](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-retained-fixation-readback-attempt01/readback-receipt.json)
independently reproduces near-zero native drift for546 observations and exact
source-restored record bytes. This is one previously sampled smoke, not1KTZ
qualification. A verifier byte-indexing failure is preserved with its explicit
technical continuation. Native output formatting clears its atom graph in source;
final source-link custody must not be inferred from native feature transport.
No weights, real-input sampling, optimizer, public profile or admission changed.
The current1KTZ request still has recorded terminal-profile refusal, undecided
water-context treatment and unclosed explicit-link custody; a scientific run must
not be silently substituted for the completed weightless checks.

The next1KTZ execution boundary is an explicit contract decision, recorded in
`tmp/mac-real-gap-continuation-2026-11-12/rfd3-conditioning-contract-consolidation-attempt01/decision-boundary.json`:
do not replace the failed native full-graph gate with source-preserving output
materialization/final topology checks without approving a separate experimental
contract. That decision alone would not waive terminal, water, scientific-admission
or bounded-run gates. Current1KTZ execution stays blocked; the177 other pending
reviews and broader work packages are not declared complete.

The user subsequently approved **developing that separate plan**. The
[experimental source-link custody plan](rfd3-source-link-custody.md) now defines
distinct native-representation, source-link-custody and scientific verdicts, a
generated-only materialization boundary, independent final-byte/endpoint checks
and negative controls. Approval is for planning: no adapter is implemented or
model job authorized by this update. The original18-edge failure remains failed.
Source inspection confirmed that the shared materializer can overwrite fixed
coordinates if they are supplied on the candidate axis; a private wrapper must
therefore audit raw fidelity first and materialize only requested missing atoms,
then independently require original observation bytes. Offline implementation
and qualification, terminal/water admission, raw precision and a separate bounded
execution declaration remain pending. No public schema/profile or frozen request
has changed, and the177 pending case reviews are not scientific admissions.

Subsequent user authorization now implements the private offline adapter in
`rfd3_custody_records.py`, `rfd3_custody.py`, `rfd3_custody_validation.py` and
`rfd3_custody_geometry.py` below `src/dvbfixer/model/diffusion/`. No runner or
public capability is registered. Source/request/context bindings, exact canonical
ownership, raw drift/one proper frame, generated-only materialization and final
observation/link readback are enforced. Existing generated/junction/global-L
diagnostics are separate from request-declared fixed-link length metrics;canonical
warnings on immutable links require review, not automatic private acceptance.

The [corrected full1KTZ offline readback](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-custody-offline-attempt02/readback-receipt.json)
independently preserves2,986 source observation records and18 explicit S–S pairs.
Its710 inserted atoms remain **old zero-coordinate placeholders**, not model
output;the PDB is labelled `NOT-A-PREDICTION` and has no whole-candidate geometry
acceptance. Fixed-endpoint length metrics emit zero findings for the18 source
pairs, not a pass for the710 placeholders. The first supplemental valid CONECT
control exposed a serial-map integration error;source snapshots/error note and
the original SS-only control remain retained. Explicit attempt02 corrects that
canonical-parser call and adds a valid reciprocal-CONECT regression without
changing invalid-record checks or original1KTZ inputs.

A [historical1UBQ geometry replay](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-custody-geometry-replay-attempt01/readback-receipt.json)
re-materializes only the already retained raw proposal. The complete final PDB
is byte-identical to the old source-restored output;546 observations and28 old
proposed atoms reproduce existing geometry and global-L acceptance. This is not
new sampling,an untouched case or1KTZ qualification. All native graph failures,
terminal/water decisions,future raw precision and scientific-admission gates
remain unchanged. The offline test matrix has45 new controls and64 focused tests
including existing materialization/metadata cases;P3/P4/P5 remain unclosed.
The final fast suite after declared-link geometry controls passes2,745 tests
(three optional skips,33 slow deselected,two SWIG warnings) in301.67s. Earlier
attempt logs remain historical,not substitutes for this corrected-code result.

### Historical execution notes — superseded by the state above

The launch/queue descriptions below record earlier boundaries, not current jobs.

Renewed autonomous continuation: the user explicitly requested continuing the
whole handoff after the bounded two-call tranche, not stopping at its bookkeeping
boundary. A separate seed11 sampler-return transport study is now declared for
1NON/2BE1, one new call each, same200 updates/one cycle/context and unchanged
RSS8GiB/MPS6GiB/1800s guards. Seed7 refusals remain immutable, not retried.
`scripts/retained_sampler_return_capture.py` passively saves the actual sampler
return/full native identity axis/200-step trace and returns the same tensor to
unchanged downstream heads. No confidence stage is skipped. Full-inference resource
failure and independently assessable raw transport are separate outcomes. Both
new calls are now accounted: full inference resource-refused, but2 actual raw
endpoints/200-step traces/full native axes survived and were independently
materialized. Observations remain byte-exact and global L holds; raw geometry
passes0/2 (junction/angle/amide/clash defects). These are transport successes, not
native inference or geometry successes.

Separate bounded CPU jobs now refine the saved outputs, with zero new sampling:
all-generated-only first, one sequence±2 anchor release only if geometry still
fails,1000 iterations/stage, no restarts/solvent and unchanged source references/
outer0.5/0.5/1.25Å bounds. Independent repair readback is prepared. A data-only
actual native peptide-graph audit is declared to test continuity at intended
junctions; mask coverage alone is not evidence of that continuity. A separately
declared narrow topology readback covers the two accepted1HE8 endpoints and all
their saved frames. The attempted12-arm read-only declaration hit outer120s limits
before analysis; ten other arms remain pending, not silently qualified.
Recruitment, limit calibration and public integration remain paused. Continue
these independent tests; do not call an exhausted small development budget
completion of WP1–WP7. Historical boundary status/docs are snapshotted before
their live continuation updates.

Historical seed7-tranche execution snapshot (superseded by the seed11 continuation
above, with original outcomes retained): the autonomous correctness
queue is complete (2546 fast tests passed, three optional skips,33 slow deselected).
All233 source tasks are accounted for:194 ready and independently axis-checked,
31 refused, eight retained time blocks, zero pending execution. These are inventory
outcomes, not cohort admission:12 small controls and the explicit7N8X exclusion
leave181 loop/construct/context reviews pending. Recruitment and input-limit
calibration remain paused. The historical mandatory long-gap quotas stay superseded.

Two local-context Protenix controls (1NON/2BE1, seed7,200 updates, one cycle,
354/365 tokens, no repair) are separately declared technical development calls.
Native data/mask and synthetic MPS reinjection checks passed. Attempts01/02 failed
before runner construction (namespace collision, then omitted torch LayerNorm
setting); both retain zero model calls and no predictions. Attempt03 restored the
established MPS environment and consumed1NON seed7, then exceeded the unchanged
6GiB MPS limit (peak6740672512 bytes). It retained200 diffusion-module forwards,
but no candidate or reinjection trace; this is a resource refusal, not completed
geometry evidence. The child self-terminated; a parent group-signal permission
race prevented its final receipt. Independent failure readback and executed script
snapshots are retained. The guard now reaps a confirmed-dead child on that race,
without ignoring permission failures for a live child or changing thresholds.
Attempt04 executed only the original unexecuted2BE1 control;1NON was not repeated.
It also exceeded6GiB (peak6874939392 bytes), with confirmed child exit and a
retained resource-failure receipt. Both actual calls are consumed, zero endpoints
exist and geometry outcomes cannot be assessed. The independent final readback
binds225 evidence files under `retained-loop-context-sampling-final-readback-attempt01/`.
At that original boundary no further model calls were declared. The later seed11
study has its own new declaration, not a retry or increased resource envelope. The final
fast lane passed **2556 tests**, three optional skips,33 slow deselected in295.44s;
82 focused correctness tests also pass. Earlier final fast failures are retained:
attempt01 omitted the CLI environment bin from PATH; attempt02 hit an inherited-
pipe fixture timing failure. The fixture now avoids Python subprocess import
latency using shell/sleep/printf, with unchanged timeout/assertions and unchanged
production preflight. See the detailed chronology; an isolated diagnostic does not
prove the exact historical timeout cause.

That bounded technical tranche is settled, **not the scientific handoff**. Its
zero-active-job boundary is historical; current continuation jobs are listed above.
Every next model-executing stage needs its own explicit method/budget declaration;
do not silently replay
seed7, infer usable raw structures from200 forward-hook counts, increase resources
or resume recruitment/limit calibration. Cohort admission, raw endpoint geometry,
matched scientific comparison, topology qualification and manual donor review
remain separate unfinished gates.
Do not assume success from launch or replay a consumed seed. Independent raw
readback is required before declaring any refinement. See the chronological
`docs/research/mac-real-gap-continuation-2026-11-12.md` and
`tmp/mac-real-gap-continuation-2026-11-12/status.json` for receipts/current state.

- [x] WP0: close publication verification and validate retained evidence.
- [ ] WP1: freeze a large natural-gap cohort, original sequences and actual tasks.
- [ ] WP2: explain 1HE8 residual failures and compare joint mobility policies.
- [ ] WP3: run matched real-gap proposals, repair and compatible assembly.
- [ ] WP4: finish staged conditioning and unmatched comparator cells.
- [ ] WP5: qualify topology diagnostics and inspect optimization passage.
- [ ] WP6: complete private donor-review/selection workflow and operational handoff.
- [ ] WP7: freeze and assess the completion policy on new natural-gap cases.

These boxes refer to the new deliverables below, not to a claim that all their
underlying code is absent. Existing work and negative results must be reused.

## What success means

The task is usable gap completion, **not recovery of an unknown crystal pose**.
Primary endpoints are complete intended sequence/heavy atoms, correct junctions,
no severe clashes or D-Cα, exact nonmobile observations, bounded explicitly mobile
observations, compatible simultaneous gaps and retained topology warnings.
Reference RMSD is optional descriptive evidence for artificial-mask controls only:
it is not an acceptance gate, donor-selection objective or reason to reject a
geometrically acceptable reconstruction. Historical scoring/output is unchanged.

The small panels and independent 1PGA/2CI2 study establish transport and bounded
workflow evidence; they do not constitute a representative production benchmark.

Do not replace crystal-RMSD ranking with a different arbitrary "realism" target.
Secondary structure, compactness, model confidence and agreement between seeds
are useful annotations, not universal rewards: a flexible loop or tail may be
appropriate. An ensemble of different acceptable completions is a valid outcome.

## Historical length policy — superseded by the loop-focused amendment above

The initial policy is **no new fixed scientific upper cutoff** on the number of
missing residues, separately for internal gaps and terminal tails. Long spans are
not discarded from the denominator merely because short loops are easier.
Existing frozen requests/capabilities and historical limits are not silently
rewritten; any new route is private and separately declared.

Distinguish three different limitations:

1. **Resources:** token/atom/context limits, memory, wall time and disk. Keep these
   real and fail closed. A resource refusal is not evidence that a given residue
   length is scientifically unreconstructable. Length alone does not determine
   cost: source context, joint task size and atom count matter.
2. **Method:** a direct sampler may fail while local-context, template-assisted or
   explicitly staged construction could work. Evaluate a labeled alternative,
   rather than truncate the target or pretend the unsupported arm succeeded.
3. **Completion quality:** full inventory and safe geometry may fail after the
   declared search/repair budget. Report the attempted route and residual defects;
   do not turn that result into a universal length cutoff without evidence.

Stratify reports by internal vs N/C-terminal, length (for example1–8,9–24,25–49,
50–99 and100+ residues), anchor separation, number/contact graph of gaps, source
size and environment. These bins are reporting labels, **not admission thresholds**.
The current inventory reaches92 residues; longer authentic spans/tails should be
added prospectively if the admitted panel lacks them, not synthesized by hiding
coordinates or chosen after inspecting favorable outputs.

First prefer a direct whole-gap proposal with explicit conditioning. If it cannot
fit or consistently fails, a separately designed local-context or staged route
may be evaluated. A staged route must preserve both anchors for internal gaps,
retain every target residue, record temporary interfaces/overlaps and audit every
join. Previously generated segments remain generated provenance, never become
measured observations or an excuse to freeze all but the newest fragment. Final
assembly/refinement keeps **all generated spans mobile together**. Progressive
construction is a new method to test, not an already working feature or permission
to implement P3 under another name; a route needing excluded P3 work requires a
separate scope decision.

If data eventually justify an operational default limit, propose it separately
for gaps/tails with failure curves, supported routes, resource costs and explicit
experimental override semantics. Do not add that cap now or automatically raise
resource watermarks to avoid having one. No unlimited-time/unlimited-memory promise.

## Main cohort: 200 real incomplete structures, not 200 masks

User expansion: the main benchmark must be **substantially larger** than the
initial nine-entry pool. Target **200 distinct experimental entries/declared
assemblies with genuine natural missing regions**, screened from approximately
500–1000 source entries. This is a prospective collection target, not a claim
that200 eligible entries have already been obtained or admitted. Expand source
screening if necessary; do not fill hard strata with artificial masks or silently
reduce the final benchmark to the easiest admitted cases.

An entry counts once even if it has many chains, gaps, seeds, assemblies, derived
projections or repeated uploads. Prefer representative sequence/architecture
diversity; cluster related constructs/near-duplicate entries and keep all their
derivatives on the same split. Report entry counts and cluster counts separately:
200 entries do not automatically mean200 independent families. The existing nine
large cases belong to development; old tiny artificial controls and small1BRS do
not contribute to the200-entry primary target.

### Historical prospective collection quotas — not current execution requirements

| Longest genuine unresolved span in an entry | Target entries |
|---|---:|
| 1–8 residues | 30 |
| 9–24 residues | 50 |
| 25–49 residues | 50 |
| 50–99 residues | 40 |
| 100+ residues | 30 |
| **Total** | **200** |

Count each entry in exactly one length row using its longest genuine span, while
evaluating **all** declared spans in its final reconstruction task. These are
collection strata, not method admission limits;100+ has no finite upper bound.
Do not discard an entry's difficult gap and reclassify it by its easiest gap.
If a stratum cannot be populated with authentic sequence/observation evidence,
retain the shortfall and review the sampling frame explicitly before changing
quotas. Unknown chemistry or partial occupancy refusals are not hidden successes.

Additional **overlapping** diversity targets, not another sum to200:

- At least150 entries with at least300 observed protein residues; at least50
  with1000+ target residues. Long unobserved tails alone must not make a tiny
  observed core count as a large observed-context benchmark.
- At least100 multichain assembly tasks and100 entries with three or more natural
  missing spans; preserve complete partner context and actual task counts.
- At least60 entries with an authentic N- or C-terminal missing segment; at
  least30 of these have a tail of25+ residues. Report N/C/internal separately.
- At least40 crowded/interface-adjacent gap cases, classified from observed
  anchor/flank/partner context rather than inaccessible missing coordinates.
- Include disulfide-rich and mixed short/long-gap architectures, and retain
  unsupported cofactor/metal/glycan/covalent-context cases in the screening ledger
  without pretending the protein-only routes handle omitted chemistry.

Aim for at least400 authentic missing spans in the200 entries, but do not count
gap number as the biological sample size. Publish per-entry actual sizes and
length distributions, source preparation/refusals, and achieved versus requested
quota counts. The hard cases are required coverage, not optional demonstrations.

### Development, mobility calibration and final evaluation

Split before inspecting new generated outcomes, grouping related structures:

| Split | Entries | Use |
|---|---:|---|
| Development | 80 | Backend/context experiments, donor/repair development; includes inspected1HE8 and the old pool |
| Anchor-policy calibration | 40 | Mobility/strength/neighborhood trade-offs on separate real structures |
| Frozen final evaluation | 80 | Only after the completion policy and budgets are frozen; no outcome-driven tuning |

Historically, length, tails, complexity and architecture were to be distributed
across splits, including100+ cases. Under the current user amendment, only admitted
loop-scale cases may enter the active evaluation; no100+ reservation is mandatory.
Case admission/sequence checks before splitting are allowed; generation outcomes
from final entries are not. Previously inspected entries cannot be moved to final
evaluation because their particular gap configuration has not yet been sampled.
Maintain a prospectively ordered same-stratum reserve list. Pre-execution source
ineligibility may trigger a documented reserve replacement, preserving the original
refusal ledger; execution/geometry/resource failures may **not** trigger replacement
to improve reported success. Final failure cases stay in final denominators.

## Existing real sources: the initial development/stress pool

Retained experimental downloads already include 9CQD, 1N8Z, 1JPS, 1AHW, 1G9M,
1GC1, 1BRS, 1I2M, 1KTZ and 1HE8. Inventory the original biological assemblies,
deposited sequence/observations, actual missing spans, partner context, partial
occupancy, missing observed heavy atoms, disulfides and heterogens before choosing
new backend cells. Earlier cc94 runs on this cohort are historical evidence,
**not** new RFD3/Boltz/Protenix execution.

The retained assembly-1 preparation inventory has substantial natural-gap tasks:

| Entry | Target residues, including gaps | Chains | Natural missing spans (residue counts) |
|---|---:|---:|---|
| 1HE8 | 1131 | 2 | 31,4,21,9,92,7,6,13,9,24 |
| 1AHW | 1294 | 6 | 3,8,8,3,8,8 |
| 1N8Z | 1041 | 3 | 9,3,4,10 |
| 1G9M | 949 | 4 | 4,12,4 |
| 1GC1 | 948 | 4 | 11,13,4 |
| 1JPS | 658 | 3 | 1,4,8,4,7,8 |
| 1I2M | 618 | 2 | 7,5,39,4,6,4 |
| 9CQD | 530 | 3 | 8,14,3,3,3,2,40 |
| 1KTZ | 468 | 4 | 12,18,10,6,12,18,10,6 |

These are inventory counts from the retained preparation summary, not fresh
admission or successful reconstruction claims. In particular, 1N8Z's initial
preparation refused severe remaining geometry defects; retain the separate
follow-up lineage rather than call that initial refusal a pass. Small 1BRS
(199 target residues, two two-residue gaps) can remain an interface/control case,
not the main evidence for large-structure readiness.

Keep original-source hashes and refusals. No trimming to a tiny protein merely to
make an inference shape fit, no occupancy rewriting and no missingness invented
by reference masking. An explicitly scoped protein-only projection remains a
projection, not completion of omitted chemistry. Compare the same deposited
missing-region task across backends; where whole-context inference cannot fit,
declare source-only local conditioning separately and assess the complete assembly.

Freeze the admitted cases, actual sequence placements, generated masks, source
preparation, backend/context arms and serial pressure-v2 call budget before the
next samples. Keep six development seeds, failures and native/restored/repair
stages; do not silently expand to a full factorial or reuse validation seeds for
tuning. Reuse existing donors for the first 1HE8 relaxation investigation.

## What the anchor experiments actually did

`scripts/refine_patchwork_all_gaps.py` constructs **one simultaneous optimization**:
all generated atoms in all gaps, plus selected observed anchors, are mobile.
Other observed heavy atoms are massless/fixed but retain their interactions in
the OpenMM force-field system. Hydrogen placement is temporary; output ownership
and source-record audits are independent. Anchors have positional restraints and
explicit movement checks relative to the original source, not ratcheted references.
Bonded/soft-steric preconditioning and full-force-field minimization are separate
recorded phases. This is local energy minimization, not global conformational
sampling or a demonstrated untangling algorithm.

In the retained full 1HE8 patchwork experiment:

- Ten gaps, **1712 generated atoms +148 atoms in19 observed anchor residues**.
- Initial bond/clash/amide findings were removed, but a Ramachandran finding
  remained (also represented by an aggregate generated-region finding).
- Three observed residues exceeded the historical movement gate:
  A321 backbone RMS0.5291/Cα0.5003 Å; A353 Cα0.5020 Å; A895 backbone RMS0.5447,
  Cα0.5312 and heavy maximum1.2612 Å. Historical limits were0.5/0.5/1.25 Å.
- Therefore the rejected result does **not** mean that anchor relaxation had no
  effect. It improved the proposal but did not pass all recorded gates.

The separate richer-context two-anchor control moved **all generated atoms plus
A895/A902**, not only those two residues. It retained an amide at A901–A902 with
ω≈−144.6° and failed its separately declared tighter anchor envelope. Its limits
must not be confused with the earlier full-patchwork envelope. No RMSD-to-crystal
gate caused either rejection.

Evidence paths under the retained `tmp/` tree:

- `1HE8-patchwork-all-generated-anchors-v1/attempt-01/`:
  `protocol.json`, `movements.json`, `summary.json`, all trial structures.
- `mac-reconstruction-2026-11-11/1he8-two-anchor-engineering-feasibility/`:
  `declaration.json`, `summary.json`, independent audits and final structure.

## Next diagnostic sequence

1. Re-audit exact movable atoms, peptide/disulfide topology and fixed neighboring
   atoms for each retained 1HE8 attempt; distinguish geometry and movement failures.
2. For A901–A902 and the remaining Ramachandran center, inspect force components,
   gradients, convergence and constraint/positional-wall activity. Energy reduction
   alone does not demonstrate removal of a diagnosed junction defect.
3. Compare direct anchors against a predeclared sequence-expanded window and
   spatial neighboring-residue shell while keeping **all generated regions mobile
   together**. Existing `spatial_relaxation.py` is implementation, not proof that
   this additional controlled 1HE8 matrix has executed successfully.
4. Retain initial/final and useful intermediate geometry/topology diagnostics.
   A clean endpoint cannot certify an unsaved trajectory. Do not solve a ring
   piercing by permitting atoms to pass through the protein during optimization.
5. If a starting donor needs a different conformation, compare another compatible
   donor or separately declared local reconstruction. More minimizer iterations
   and looser anchor limits are not assumed solutions.

Do not silently relax scientific gates to pass a favorite case. If a movement or
geometry policy is changed, declare it explicitly, justify it against the actual
completion objective and independently report old/new outcomes. New whole-1HE8
donor approval remains a separate manual decision; no public capability or P3 is
enabled by this plan.

## WP0 — publication and reproducibility prerequisite

The published branch is `research/mac-reconstruction-real-gaps`, initial study
commit `b17e7af`. Archives preserve ten series,12489 files, native/restored/final
outputs and refusals. Start by validating archive parts and extracted payloads
against their SHA256 manifest; reconstruct only in an empty owned work directory.
Local checkpoint/env trees still require separate provisioning and valid licenses.

The post-publication fast rerun finished **2473 passed, one failure, three optional
skips,33 slow deselected**: `test_repository_inventory_is_current` reports stale
`docs/research/license-inventory.json`. Regenerate the inventory from the current
tracked tree, inspect the change and rerun its check plus the complete fast lane.
Do not relabel this run as all passing or discard the older passing2472-test log.
Ruff/new archive controls, agent maps and CLI reference checks are separate evidence.

Deliverable: verification receipt, synchronized generated inventory and explicit
source/runtime/checkpoint/asset provenance for the next execution manifest.

Execution update: the three newly tracked patch records have been added to the
generated inventory; its check and12 focused tests pass. The full fast lane now
passes2474 tests (three optional skips,33 slow deselected). Archive payload hashes
were independently checked at publication. See the
[continuation log](../research/mac-real-gap-continuation-2026-11-12.md) for receipts,
the completed1000-entry census/grouping,412 original-source inventories and
1HE8 temporary-force inspections. Assembly preflight produced380 inventories and
32 retained refusals. All118 initial source tasks are accounted for:101 are privately
ready,13 retain source-preparation refusals and four retain1800s timeouts; none
remain unexecuted. All101 ready placements match
deposited axes. Eight original traced1HE8 controls completed (0/8 joint passes,
3/8 geometry passes); four separately declared inward-wall controls completed
(2/4 joint passes,3/4 geometry passes). The passing richer-context sequence±2 and
spatial policies retain the unchanged outer0.5/0.5/1.25Å acceptance. This is private
protein-only endpoint evidence, not whole-entry chemistry or broader method
qualification. Additional181 originals from unused existing census entries are
inventoried;41 new source preparations are fully accounted:32 ready, eight refused
and one retained timeout. All133 ready placements across the two settled programs
match deposited axes, with seven small controls explicit. This replenishes short/medium
reserves, not long strata. The500-new-entry public metadata screen and union
grouping completed:1500 records,1016 conservative components, largest176. A new
140-original acquisition completed for63 annotated50–99 and77 annotated100+
cases without repeating prior downloads; assembly preflight completed120 inventories
and20 refusals.39 hard-source preparations are active. The141-original recruitment
for25–49 completed download/label-axis inventory, excludes all733 prior downloads
and leaves both long quotas unchanged. Its assembly preflight is declared and
waiting for the CPU executor. Scientific admission/splits remain open.
Context/construct/format/occupancy reviews and targeted long-stratum recruitment
remain open. These steps
are not final200-entry admission or new sampler execution.

## WP1 — natural-gap cohort and admission

Use the nine substantial entries in the table as the initial development/stress pool,
not the final benchmark size; screen/curate the200-entry cohort specified above.
Keep1BRS as a small interface control, not the principal production benchmark.
Audit original deposited sequence against observations, REMARK465 and assembly
metadata. Distinguish true unresolved sequence from expression tags, different
constructs, unobserved chains and ambiguous placement. An invented tail from an
incorrect FASTA/SEQRES alignment is not a reconstruction task.

Freeze per case: original source/sequence hashes, assembly/altloc policy,
case-sensitive chain/residue/insertion-code axis, every missing span, terminal vs
internal identity, anchor availability, lengths, observed completeness/occupancy,
disulfides/covalent links, heterogens, source geometry, prepared-source provenance,
input size and source-only gap interaction graph. Existing measured-atom repair
must be a separately retained preparation derivative; never silently repair frozen
replays or call template-completed observations diffusion-generated.

Require real partner/crowded context. Preserve chemistry/occupancy refusals and
the full original source; a protein-only route cannot claim whole-entry chemistry
support. Omitted heterogens/obstacles must be explicit. No chain or environmental
obstacle may be silently discarded merely to admit a model input.

Output a **task-level** matrix. An entry may contain several local sampler tasks
and multiple interacting groups, but its final endpoint remains the complete
declared assembly. Keep requested/admitted/executed/accepted counts separate at
entry, gap, task, seed and complete-assembly levels. Do not count ten easy gaps in
one successful seed as ten independent structures or ignore a failed eleventh gap.

## WP2 — causal investigation of 1HE8 joint relaxation

Begin without neural resampling, using immutable saved proposals and the five-step
diagnostic sequence above. Produce an issue-to-atom/topology/force map for A901–902,
the remaining Ramachandran center, source clashes and the known SER1044/TYR972 ring
alarm in failed donors. Separate source-only findings, generated-internal findings,
junction defects and defects introduced by fitting/restoration/refinement.

Verify movable identities against the **final** temporary topology after atom/H
rebuilds; inspect massless particles, peptide/disulfide bonds, hydrogen geometry,
positional restraints and movement walls. Report physical force-field terms,
auxiliary geometry terms, gradients on problematic atoms, displacement and whether
the run converged or merely reached its iteration ceiling. Do not assume stale
indices, force imbalance, restraints or local trapping as an established cause
before measurement.

Freeze a small comparison on the same starting proposals:

- fixed observations (control);
- all generated atoms plus direct observed anchor residues;
- all generated atoms plus a prospectively sized sequence-expanded anchor window;
- all generated atoms plus an explicitly bounded spatial neighboring-residue shell.

Use the same original observation references and unchanged common heavy-atom
ownership. Compare a small preregistered restraint-strength/movement-envelope set;
do not repeatedly widen windows until a favored structure passes. Full generated
regions remain jointly mobile in every arm, including otherwise clean gaps.
Mobility selection is a physical neighborhood policy, not an ERROR-only atom list.

If tiny gate exceedances explain a rejection, show actual values and classify
them separately from severe geometry; a new operational tolerance needs an
explicit reviewed policy. If changing the mobility neighborhood removes the
defect, report its effect on all observed contacts and other gaps, not only the
chosen junction. If no local minimum solves it, test a different compatible donor
or a separately labeled existing reconstruction route, not limitless minimization.

Deliverable: an explanation supported by controls, per-arm outcomes and a bounded
joint repair policy, **or an explicit unsolved defect**. A passing outcome alone
does not establish causality, and successful repair is not physical uncertainty
calibration.

## WP3 — matched proposals, shared repair and assembly on real gaps

Compare RFD3, Protenix v1, stock Boltz proposal mode, the explicit MODELLER protocol
and cc94 **where each admits the same task**. Do not silently swap weights or give
one backend a much smaller/easier gap. Retain model-specific constraints and
observed terminal-OXT semantics; Boltz/MODELLER drift before restoration remains
separate from final observation custody. cc94 refusal is an outcome, not grounds
for cropping the public/frozen request to force admission.

Every cell retains native proposal, identity/frame transport, source restoration,
raw findings, necessary common repair and independent final readback. Clean
single proposals need no arbitrary rebuilding; selected multi-donor assemblies
undergo the mandatory shared all-generated complex operation. Anchor release is
necessity-gated under the frozen WP2 policy. Final geometry/custody is measured in
the complete declared context, including obstacles omitted from a local model input.

Compare independent gaps, interacting joint pairs/small groups and large connected
groups. Keep each jointly sampled group atomic during donor selection. A resource-
driven split changes the method and needs its own declaration/lineage. Bounded
compatible search must retain locks/allow-lists, search budget, ties, donor/repair
derivatives and whole-assembly rejection reasons; local passes never imply a valid
whole1HE8. Compare homogeneous single-seed assembly against compatible mixtures
without using hidden-reference scores. Preserve the old manual selections and
require separate review for new1HE8 donors.

Deliverable: per-entry/per-gap/per-seed raw and final geometry, source movement,
topology flags, resource cost and complete-assembly success. Diversity and
secondary-structure diagnostics are descriptive. Crystal RMSD is not required for
these genuinely unknown regions.

## WP4 — finish the original comparator and conditioning backlog

| Outstanding cell | Required continuation |
|---|---|
| A: context | Compare flanks, spatial neighbors, local domain and full complex on actual large-gap tasks; retain complete assembly for validation |
| B: fixation | Extend template vs per-update enforcement beyond pilots; study soft coordinate restraints only where the pinned interface actually supports them; quantify drift and restoration damage |
| C: atom information | Repeat backbone-only vs all-heavy context on crowded real cases, with matched sequence/residue coverage |
| D: coupling | Compare independent, pair, small group and difficult large group; distinguish changed task size from coupling itself |
| E: MSA/templates | Add frozen local MSA and additional homolog templates for Protenix/Boltz; keep no-MSA/source-template/no-template controls and actual consumer receipts |
| cc94 comparator | Include matched admitted cells with its original500-update profile; no inference from unrelated historical runs |
| Long synthetic controls | Complete missing five ML seeds for the existing16-residue panels if those historical comparisons remain in scope; MODELLER already has six, but these controls must not displace the natural-gap panel |
| Interfaces | Admit and execute real complete partner-context cases, or retain explicit source/chemistry/resource refusals; earlier refusal alone is not a finished comparison |

Select informative cases/arms prospectively; do not run the full factorial. Local
MSA/homolog provenance must bind sequence, search/database version, alignment,
coverage and template coordinates. Do not fetch a structure that directly supplies
an artificial control's hidden reference. A homolog may legitimately guide unknown
real gaps, but this is template-assisted completion, not de novo reconstruction.
No remote sequence upload without explicit authorization. Unsupported native mask
or soft-restraint interfaces remain unsupported, not renamed successful arms.

## WP5 — topology and optimization-path evidence

Extend the existing explicit-topology segment/ring/surface diagnostics with
positive/negative, ambiguous, rigid-transform, legitimate linked/native-topology,
disulfide and chain-break controls. Preserve actual covalent connectivity, chain
case and insertion codes; never infer peptide links merely from spatial proximity.
Run diagnostics on source, raw, restored, assembled and final structures.

For controlled repair trajectories retain meaningful intermediate states with
exact optimizer-stage metadata. Instrumentation must not silently change the
method: repeatedly restarting/chunking a minimizer is a separate arm unless
equivalence is demonstrated. Saved frames can expose passage but sparse clean
frames cannot prove its absence. Open-path closure ambiguity and unqualified
alarms remain visible; no blanket ban on naturally knotted/linked configurations.

Deliverable: reviewed diagnostic scope, reproducible implicated bonds/residues and
viewer selections, measured false alarms and explicit limits. A fatal policy
requires separate qualification; until then do not claim a topology guarantee.

## WP6 — private donor library, review and operational delivery

Unify the existing library/search scripts into a reproducible private workflow:
overlay candidates on the same source, filter by entry/gap/joint-group/backend/
seed/stage, inspect junctions/clashes/topology/movement, set allowed donors or locks,
show conflicts and export a versioned selection recipe. Clearly distinguish
native, restored, repaired and manually accepted derivatives. Register generated
ownership using canonical `.dat` semantics where applicable; known atoms and
prepared-source additions must not become diffusion-generated.

Scientific materialization stays in Python. Manual preferences are review records,
not scientific acceptance. Preserve source/request/asset hashes, joint lineage,
exclusive output reservations, failure/cancellation evidence and validated resume.
No crash-left directory adoption or overwrite of original trials. Finish the
private recipe/viewer surface before proposing a GUI slice; no automatic public
CLI/Homology/GUI/default change. If product integration is later approved, its
workspace isolation, autosaves, provenance and actual browser/GUI tests are a
separate deliverable, not inferred from this research.

## WP7 — freeze the completion policy and assess new real cases

Freeze input preparation, routing/context, donor ranking/tie-breaks, joint repair,
anchor limits, topology treatment and per-task budgets before final evaluation.
Use the80 prospectively reserved natural-gap entries not inspected during tuning,
spanning internal/tail, long/multiple/crowded cases. All nine current pool
entries are already inspected research candidates, so relabeling them independent
validation would be incorrect. Previously used validation seeds211/223/227/229/233
are no longer untouched. Declare a fresh disjoint final-evaluation seed set, record
first actual consumption and do not tune from its outcomes.

Assess completion geometry/custody, movement, unresolved topology evidence,
complete-assembly success and resource cost—not unknowable true missing-region
coordinates. Report every requested case and every failure; distinguish source
admission, model input transport, resources, numerical repair and geometry failures.
Keep seeds/gaps nested under entries, not independent biological replicates.
Qualification may remain conditional; no length cap or public enablement follows
automatically from a small successful cohort.

## Execution order, budgets and stop conditions

1. WP0 and weightless WP1 admission/inventory first.
2. WP2 uses retained1HE8 structures before new neural samples; freeze its exact
   CPU arm/proposal/optimizer count and resource ceilings in a new declaration.
3. Select a frozen10-entry technical tranche **within development**, including
   difficult long/multigap cases and1HE8, not just the easiest shapes. Execute one
   declared development seed per admitted backend/task. This is an operational
   pilot, not an acceptable endpoint for the200-entry benchmark.
4. Freeze amendments after technical pilots, then cover all80 development and40
   calibration entries in declared batches, preserving the80-entry final holdout.
   Proposed base candidate counts: six development seeds7/42/101/102/103/104,
   three fixed calibration seeds7/42/101 and five freshly declared disjoint final
   seeds. Reuse valid pilot outputs rather than resampling seed7. Keep failures;
   no adaptive best-seed/backend choice based on reference coordinates.
5. Stage WP4 subsets and WP5 controls separately, with their own incremental costs;
   implement WP6 reusable custody/review while outcomes accumulate.
6. Freeze WP7 before new evaluation cases/seeds are executed.

For cost visibility, the proposed base matrix at **one native task per entry per
backend** is `80×6 +40×3 +80×5 =1000` proposals per backend. RFD3/Protenix/Boltz/
cc94/MODELLER would therefore reach **up to5000 proposal jobs** if all cells admit,
before conditioning ablations, multiple local tasks, donor search and CPU repair.
This is not a small Mac smoke and **not blanket authorization to launch5000 jobs**.
Even a uniform30-minute per-job ceiling would sum to2500 serial hours, not an
expected runtime. Measure actual costs on the10-entry pilot, enumerate admitted
native tasks and publish per-batch time/storage/compute proposals before launch.
If the host cannot support this matrix, explicitly choose scheduling/hardware or
amend per-cell replication with the user; do not quietly shrink200 entries to ten,
remove long-gap strata, increase concurrency beyond safe limits or spend an
unbounded cloud/GPU budget. Changes of hardware/runtime need their own acceptance.

For every stage compute `calls = sum(native tasks per cell × new seeds)` and list
CPU repair/search trials separately. Preserve one-ML-process concurrency,
pressure-v2 RSS8GiB/MPS6GiB, host headroom/pressure checks and explicit wall/disk
limits. Existing ML per-call1800s is the starting declared ceiling, not proof that
a long-gap method is impossible if it times out. Any justified budget change is
explicit and prospective; do not disable watermarks, silently crop context, use a
generic CPU fallback or spend an unbounded Cartesian product. Retain the pinned
model numerical methods and per-backend update counts.

Stop remaining launches on technical/identity/resource failures and retain valid
completed native outputs for an explicit recovery amendment without resampling.
Scientific geometry failures stay in denominators; continuing the remaining
unchanged declared scientific cells is not a retry. Manual donor approvals,
remote services/new assets, excluded chemistry/P3 and product integration remain
separate boundaries. No new checkpoint sampling has been performed for this plan.

## Publication and navigation

Code/tests/reports and checksum-bound research archives are published on a new
branch from `main`; local environments, weights, browser profiles and unrelated
force-field directories are not part of that publication. See the
[comparison report](../research/mac-reconstruction-comparison-2026-11-11.md),
[original handoff](mac-reconstruction-portfolios-handoff.md) and
[evidence archive instructions](../research/evidence/mac-reconstruction-2026-11-12/README.md).
