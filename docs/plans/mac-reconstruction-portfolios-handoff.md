# Mac reconstruction backends and local multi-seed portfolios: handoff

Date: 2026-11-11. **Design/documentation approved; implementation, installation,
weight downloads and experimental execution require separate user permission.**
Do not interpret approval to write this handoff, or previous 1HE8 runs, as
permission to execute the new program below. No work package is marked complete
merely because an upstream project advertises a capability.

Execution update (2026-11-11): the user subsequently authorized autonomous
execution. Three native-MPS backends, single-case six-seed controls, explicit RFD3
seed-transport correction, source-only local 1HE8 conditioning, constrained donor
selection and unqualified surface diagnostics are recorded in the
[feasibility log](../research/mac-reconstruction-feasibility-2026-11-11.md).
This is not full-plan completion or approval to waive the resource/capability
gates. The authorization statements below describe the original design boundary;
new runs remain bounded by their separately retained declarations.

## 1. Objective and boundaries

Evaluate RFdiffusion3 (RFD3), Protenix v1 and Boltz-2 on Apple Silicon for a
source-preserving reconstruction workflow:

```text
source + known sequence + declared missing regions
  -> independent local multi-seed proposals (or explicitly joint small groups)
  -> candidate library with immutable provenance
  -> manual/automatic compatible donor selection
  -> assembled patchwork and junction checks
  -> joint minimization of ALL generated atoms
  -> bounded observed-anchor release only if necessary
  -> independent final geometry, motion and interaction assessment
```

Initial scientific scope: canonical, protein-only inputs; whole missing residues
with known sequence; all canonical heavy atoms including sidechains. Hydrogen
placement/protonation is separate. Do not silently change sequence or treat
backbone-only completion as success. Partial-residue atom completion, ligands,
glycans, PTMs, metals, covalent heterogen context and production-default changes
are out of this first scope. P3 alternative-reconstruction work remains excluded.

A complete heavy-atom prediction is not by itself an inpainting interface.
Distinguish strict source-preserving reconstruction from a new domain prediction
used as a proposal. MODELLER remains the production baseline. Protpardelle cc94
local-motif donors remain a comparison arm, not deprecated by this plan.

Related evidence: [1HE8 experiment record](../research/1he8-motif-portfolios-and-joint-generation.md),
[geometry repair plan](diffusion-geometry-repair.md),
[anchor-relaxation handoff](geometry-gated-anchor-relaxation-handoff.md),
[modeling backend research](../research/reconstruction-and-modeling-backends.md).

## 2. Starting state and authorization

- Work from the repository root, not the historical worktree. The inspected host
  is an Apple M3 Pro with **18 GB unified memory**. Recheck host/runtime identity
  at execution time; do not copy hardware identifiers into public reports.
- Existing runtime assets include `dvbfixer/engine`, `dvbfixer/sampling`,
  `dvbfixer/refinement` and `.artifacts/apple-cc94/model_params/`. They are local
  assets, not guaranteed to arrive through Git. Do not modify pinned engines,
  existing environments or repair damaged historical prefixes in place.
- Important experimental code/evidence is currently uncommitted/untracked.
  Inspect Git status and preserve it. Do not assume a fresh clone contains it.
- Original five-seed outcomes, original 1HE8 inputs and negative joint-generation
  results must remain intact. Store new run artifacts inside the project.
- Before execution, present the chosen revisions/checkpoints, input manifest,
  concrete matrix, total sample counts, resource caps and storage estimate to the
  user. Record the authorization and its scope. Stage expansions need their own
  approved manifest; do not expand a budget merely because a preceding run fits.
- Until that authorization: documentation and design only. Do not install,
  download, patch engines, sample, refine or launch background experiments.

## 3. Backend admission: test the actual scientific interface

| Route | Intended experiment | Required admission proof |
|---|---|---|
| RFD3 | Fixed-sequence reconstruction with observed atom constraints | Exact target sequence including unknown-coordinate residues; actual fixed-atom mask; canonical heavy-atom output; bounded source motion; no inadvertent sequence design |
| Protenix v1 | Existing coordinate-reinjection idea ported to MPS | Compatible checkpoint and callback; source-coordinate enforcement during every update; native vs restored-frame evidence; no silent numerical/backend substitutions |
| Boltz-2 | Template/constraint-conditioned alternatives, and possible reconstruction adapter | Establish whether observations can actually be fixed. If not, label domain-prediction/proposal mode separately; templates alone do not prove immutability |

Candidate upstream entry points, to resolve and pin before installation:

- [Foundry](https://github.com/RosettaCommons/foundry) documents community MPS
  support through [fnachon/foundry](https://github.com/fnachon/foundry).
  [RFD3 input specification](https://rosettacommons.github.io/foundry/models/rfd3/input.html)
  describes fixed atoms and sequence controls. Confirm behavior in the exact
  checkpoint/runtime rather than inferring it from field names.
- [Protenix](https://github.com/bytedance/Protenix),
  [model inventory](https://github.com/bytedance/Protenix/blob/main/docs/supported_models.md),
  and [community MPS fork](https://github.com/DeepFoldProtein/Protenix-MPS).
  Check fork/checkpoint compatibility, actual inference code and per-update
  callback availability. Mini/Tiny are separate checkpoints/arms, not transparent
  memory optimizations of Base; do not assume their template features match.
- [Boltz-community](https://github.com/Novel-Therapeutics/boltz-community)
  documents MPS inference. Pin the fork separately from upstream Boltz and record
  numerical/behavioral changes.

This is not permission to install all forks. Review code, licensing, checkpoint
terms, native dependencies and version compatibility first. Protenix v2 is not
an assumed freely distributable alternative: upstream explicitly distinguishes
its restricted weight terms. Do not let a repository-wide license override
checkpoint-specific restrictions.

Each admitted backend gets an isolated environment and an owned clone/patchset.
Use native macOS execution for MPS, not an assumption of Metal access in Docker.
Retain existing Protenix CUDA process separation and frozen historical engines.
Unsupported chemistry or unavailable strict conditioning is an explicit outcome,
not an opportunity for silent cropping or source restoration disguised as fixing.

## 4. Seeds, cohorts and cost accounting

### Fixed seeds

- **Development/comparison: six seeds `7, 42, 101, 102, 103, 104`.** The original
  five are preserved and 42 is the user-requested addition.
- **Post-freeze validation: five seeds `211, 223, 227, 229, 233`.** These are
  disjoint from development seeds. Do not use them to tune selection or thresholds.
- Six candidates per local gap and backend/conditioning arm are the initial
  portfolio, not six accepted structures or evidence of exhaustive coverage.
- Same integer seeds do not imply matched noise across different architectures.
  Compare controlled inputs/budgets and distributions, not seed-paired causality.
- Optimizer restarts are not new diffusion samples. Log their own fixed seed,
  budget and perturbation policy separately. Failed seeds stay in the denominator;
  technical retries have new attempt IDs, retained failures and unchanged settings
  unless a protocol amendment explicitly states otherwise.

### Cohorts

1. A small native-MPS functional smoke with exact masks and an inspectable gap.
2. A development panel covering short internal loops, longer structured regions,
   termini, close loop pairs and interface context. Include 1HE8 only as a known,
   retrospectively examined stress case; it is not independent validation.
3. Source-complete structures with synthetically hidden coordinates for calibration
   and reference scoring. Reference coordinates must never enter conditioning,
   assembly selection or refinement. Preserve original observations/occupancies.
4. Separate structures and the validation seed set after freezing the pipeline.

Before sampling, freeze structure IDs, checksums, assembly/altloc choices,
case-sensitive chain/residue identities, sequence placements, masks, context
selection and train/validation split. Consider upstream training cutoffs: a
held-out workflow case is not necessarily absent from a model's training data.
Do not remove hard cases post hoc or rewrite occupancies to bypass admission.

Cost formula for independent gaps is `backends × conditioning arms × gaps × 6`
development calls. A joint group counts once per seed, but every member and its
shared lineage must be recorded. Freeze the actual number of admitted cells
after small smokes; do not launch an unbounded Cartesian product of all studies.

## 5. Mac feasibility and memory gates

The earlier 512-token public budget and private 1131-token ceiling apply to
specific Protpardelle experiments, not automatically to these models. Begin with
a small input, then a few increasing context sizes; freeze each tested shape.
Fit the same scientific input or report unsupported/over-budget, rather than
silently discarding observations to produce a success.

Before launch, specify per-process and whole-host memory policies, wall timeout,
disk quota, concurrency and graceful interruption. First use one ML process at
a time. On 18 GB unified memory, leave explicit headroom for macOS and the GUI;
RSS and MPS allocation overlap and must not be naively added as independent RAM.
Monitor swap/memory pressure as well as allocator peaks. Disable automatic CPU
fallback for the MPS qualification arm; a deliberate mixed/CPU arm needs a
separate label and accounting. Do not disable safety watermarks to force a pass.

Record actual device placement of representative tensors/kernels, precision,
finite outputs, measured peaks, runtime and all warnings. Plain Torch kernels
may be needed instead of CUDA-only extensions. Precision selection is pinned
per backend; no unvalidated global FP16/BF16 conversion. OOM or missing MPS
operators are retained technical failures, not scientific geometry rejections.

MSA and template searches are a separate provisioning/data step. Default to local
or cached assets; no sequence upload to a remote MSA service without explicit
permission. Record whether a nominally local workflow used a remote service.

## 6. Conditioning studies

Run staged, one-factor comparisons before testing selected interactions. Use the
same source/mask, candidate count and downstream assessment within a comparison.
Report architecture-specific mechanisms instead of pretending all backends
implement the same interventions.

| Series | Arms | Question and measurements |
|---|---|---|
| A: context | Sequence flanks; flanks plus spatial neighbors; local domain; full complex where feasible | Which observations matter? Record coverage, omitted nearby obstacles, token/atom counts, native quality, source drift and final quality. Full context is an arm, not the new default. |
| B: fixation | Ordinary template/conditioning; per-update coordinate enforcement; soft coordinate restraints where supported | Measure drift before restoration, restoration-induced junction damage and fixed-mask fidelity. Do not equate a template with a hard constraint. |
| C: atom information | Observed backbone vs all observed heavy atoms | Is sidechain obstacle information helpful? Keep sequence/selected residues otherwise matched; retain full source for final independent assessment. |
| D: joint scope | Independent gaps; a declared interacting pair; a small group; a large connected group as a stress arm | Separate benefits of coupled generation from increased task size. Use joint group lineage; never call sequential independent sampling a joint call. |
| E: sequence evidence | No MSA; fixed local/cached MSA; observed-structure template; additional homolog template | Compare within Protenix/Boltz capabilities. Hold MSA/template assets fixed; record unsupported cells rather than substituting checkpoints without disclosure. |

Define flank lengths, neighbor radius/ranking, maximum shape and physical sequence
offsets prospectively. Inspect pinned consumer masks/placement before weights are
loaded. Local projections must be explicit and bind the original structure;
never crop a frozen replay silently. Do not feed a prior generated prediction as
observed conditioning; any such future arm would need separately labeled ownership.

For D, retain the historical 8-observed-residue flank / ≤6 Å heavy-contact graph
as a comparator. Its transitive closure grouped nine 1HE8 gaps into one 704-token
problem and did not solve the problem. Investigate smaller interaction tasks,
but define them as an explicit new policy, not an undocumented split of a frozen
component to satisfy a memory budget. If a component cannot be handled, report it.

Save native output before fitting/restoration and the materialized source-frame
proposal separately. This is necessary to distinguish model defects, context
drift and final source-restoration artifacts.

## 7. Candidate library, ownership and manual selection

One library entry binds source/request digests, gap or joint-group identities,
backend/checkpoint/patches, seed, conditioning assets, native arrays, raw PDB,
all generated atoms, assessments and optional refined derivatives. Stages must
not be mixed: a raw prediction and a generated-only-refined donor are distinct
entries/lineage nodes. Preserve failed candidates for inspection.

Reuse the principles in `gap_portfolio.py`, `portfolio_evidence.py`,
`portfolio_publication.py` and `DatRecord`; inspect exact interfaces before
implementation. Do not silently extend frozen public schema-4 or relabel private
portfolio bundles as ordinary single-seed results. Preserve the matching `.dat`
generated ownership, template provenance, chain case and insertion codes.

Manual workflow requirements:

- Switch/overlay candidates against the same observed source; show gap, seed,
  stage, backend, secondary-structure diagnostics, junctions and clashes.
- Permit one selected donor, a set of allowed donors, or an explicit lock.
  Show conflicts with other selections and which repairs moved a fragment.
- Keep the user's preferences as versioned review records. Existing 1HE8 visual
  choices (36/50) are development labels, not acceptance or reference truth.
- Export a reproducible selection recipe. Scientific materialization and
  validation belong in Python, not duplicated in Node/React.
- Inspect the existing persistent Homology/template-plan patterns, but do not
  confuse a generated-candidate library with the template mosaic implementation.
  Preserve workspace isolation, revision-safe autosaves and artifact custody.

For a jointly generated group, initially choose the **whole group from one run**.
Mixing its members across seeds changes their coupled geometry and must not be
represented as the original joint proposal. A separately enabled recombination
would require a new recipe and complete compatibility assessment.

## 8. Automatic compatible-combination selection

Define separate local, pairwise and whole-assembly assessment layers:

1. Exact sequence/chemistry/inventory, finite coordinates, generated ownership,
   chirality, backbone continuity, local stereochemistry and junction eligibility.
2. Candidate-vs-source and candidate-vs-candidate sterics, backbone near-crossings,
   suspected threading and compatibility of feasible junction/anchor motions.
3. Full assembled-complex validation: pairwise passes are not proof of collective
   compatibility or correct topology.

Build a conflict/compatibility graph over donor alternatives. For small spaces,
enumerate all combinations; otherwise use a deterministic bounded beam or a
documented constraint solver. Report explored combinations, budget exhaustion,
ties and whether optimality is established. User locks/filter sets constrain the
same search, not a separate unvalidated materializer. Compare with the existing
64-evaluation / width-4 portfolio baseline and homogeneous single-seed controls.

Do not optimize a single unexplained sum of incomparable scores. Keep admission,
geometric defects, motion eligibility, ambiguity and structural-quality evidence
separate. Secondary-structure assignment/coherence and ensemble stability are
diagnostics, not a universal reward for helices: loops/disordered regions may be
correct. A model's own confidence is not independent validation. Freeze ranking
and tie-breaks before held-out evaluation.

## 9. Intersection/threading assessment

Implement the separately requested
[detector plan](diffusion-geometry-repair.md#generated-loop-intersectionthreading-detector-requested-follow-up-2026-11-10).
Distinguish close contacts, 3D backbone-segment intersection/near-crossing and
suspected topological threading. Open paths lack a unique knot/link invariant
without a closure convention. Native knots/threading exist: do not declare every
visually intertwined pair impossible or enforce a blanket unknotted-protein rule.

Check generated–generated and generated–observed paths, actual peptide continuity,
chain breaks, disulfides and case-sensitive identities. Treat inherited source-only
findings separately. Require positive/negative controls, rigid-transform invariance,
projected-but-not-3D crossings, legitimate topology and explicit ambiguous cases.
Expose exact segment pairs and viewer selections, method/version and uncertainty.

Run at raw, assembled and final stages. Investigate intermediate minimization
states for chain passage; endpoint-only clash checks cannot prove an optimization
trajectory never crossed. Initial experimental alarms are not silently new fatal
public gates. A versioned reviewed policy and validated controls are required.

## 10. Junction checks, optional repair and shared minimization

**Checking junctions is mandatory; rebuilding/closing them is not.** Coordinates
inserted into a PDB can still have wrong C–N distances, angles or peptide
planarity, but an already valid junction requires no separate closure procedure.
Closure also does not untangle mutually threaded loops.

Protocol:

1. Materialize the selected assembly in the original observed frame and assess
   every junction. Preserve this pre-minimization assembly.
2. Perform shared generated-only minimization: **every generated atom of every
   selected gap/group participates in the same complex**, with observations fixed.
   Do not freeze clean generated gaps merely because only one loop has a finding.
3. If residual junction defects need observed motion, retry with an explicitly
   selected bounded anchor window while keeping ALL generated atoms mobile.
   Reference coordinates always come from the original source.
4. Independently validate whole-complex geometry, outside-mask records,
   calibrated per-anchor motion and interaction/threading findings. Retain both
   successful and rejected attempts; no publication on independent failure.
5. If a junction remains unclosable, choose a different donor combination or
   perform a separately labeled optional local-repair experiment. First evaluate
   the existing licensed MODELLER route rather than automatically implementing
   CCD/KIC. Verify its actual movable mask and whether it rebuilds an entire
   otherwise useful fragment. Record any changed fragment as a derivative.

Only investigate another constraint-based closure algorithm if existing repair
is inadequate and a separate scope decision authorizes it. No mandatory MODELLER
replacement/reconstruction of already valid generated segments is intended.

Use the existing `boundary_refinement.py`, `geometry_conditioning.py`,
`runner_refinement.py`, `anchor_policy.py` and `staged_runtime.py` where appropriate.
Inspect their current masks: the ordinary whole-candidate refiner and the older
issue-local **post-portfolio** repair are not the same stage. Explicitly integrate
the shared post-assembly operation; a private experiment script alone does not
make it mandatory public behavior. If invoking `minimize`/`zbs`, pass
`--no-solvent`; preserve `keepIds=True` and all existing chirality safeguards.

## 11. Calibrate anchor mobility rather than ratifying arbitrary limits

The historical `(backbone RMSD 0.5 Å, Cα 0.5 Å, heavy maximum 1.25 Å)` limits
are engineering controls, **not established physical uncertainty bounds**.
Keep them unchanged for historical replay; new calibration is a distinct arm.

Use complete reference structures with preregistered hidden regions. Include
short/long loops, termini, secondary-structure boundaries and crowded interfaces.
Split structures/sequence families as well as seeds between calibration and
validation. Hide reference coordinates from generation, selection and refinement;
score against them only in the independent evaluation stage.

Evaluate incrementally, not an unrestricted parameter sweep:

- Fixed observations as a zero-motion control.
- Direct flanking residues versus a prospectively bounded expanded window;
  do not recursively release an entire chain or pick mobility only from ERRORs.
- A small declared range of positional restraint strengths, followed by a
  separate movement-envelope comparison. Candidate engineering grid for review:
  `(BB RMSD, Cα displacement, heavy maximum)` in Å of
  `(0.10, 0.10, 0.50)`, `(0.25, 0.25, 0.75)`, `(0.50, 0.50, 1.25)`.
  These are **experimental grid points, not approved production limits**. Freeze
  the final grid, strengths and sample budget before running calibration; the
  historical setting is a comparator, not the desired winner.
- Initially use uniform restraint weights to isolate motion policy. Study
  B-factor weighting separately only with provenance checks. B factors include
  disorder/occupancy/model effects and cannot directly define allowed motion.

Measure per-anchor and whole-input displacement, held-out generated-atom/backbone
error, junction quality, secondary-structure/contact damage, acceptance rate and
repair-induced crossing alarms. Report Pareto trade-offs and uncertainty. Do not
pick the loosest setting because it yields more geometry passes. Account for
coordinate precision/resolution and errors in the observed structure; a deposited
reference is not exact ground truth. Finalize a conservative policy on calibration
data, then assess it unchanged on the validation cohort. An unresolved result is
acceptable; retrospective threshold tuning to accept 1HE8 is not.

## 12. Deliverables, gates and execution order

| Gate | Deliverable / acceptance | Stop or defer when |
|---|---|---|
| G0: authorization | User approves exact scope, environment/checkpoint choices, cohort and resource manifest | This document alone is not execution approval |
| G1: isolated runtimes | Pinned license-reviewed installs, real MPS checkpoint smoke, measured resources, preserved existing environments | Unsupported operators, unsafe memory pressure, unavailable terms or incorrect device placement |
| G2: reconstruction admission | Exact sequence/masks/atom output; native/restored-frame audit; observed preservation or explicit proposal-only classification | Merely producing an all-atom PDB is not enough |
| G3: baseline library | Six declared seeds per admitted development task; failures and native/raw/derivative stages retained | No silent retries, sample-count expansion or test-set tuning |
| G4: conditioning comparison | Staged A–E subset with matched inputs and explicit unsupported cells | Do not launch full factorial without an approved count/cost |
| G5: compatible assembly | Manual recipes and deterministic auto search, whole-group lineage, interaction controls, whole-complex audit | Local cleanliness/manual preference is not acceptance |
| G6: shared refinement and calibration | All-generated joint minimization, necessary-only anchor release, reference-blind calibration and unchanged held-out test | No arbitrary motion-limit relaxation or mandatory closure of valid stubs |
| G7: product decision | Independent reports, retained negatives, model-specific readiness; separately approved CLI/Homology/GUI integration | Technical smoke success does not qualify broad public support |

Implement/test pure contracts, masks and geometry controls before costly sampling.
Use no-sampling fixtures/fake adapters for library/search/recipe tests. Run narrow
tests first, then the fast Python suite in a complete scientific environment.
GUI work must also pass `npm run typecheck` and `npm test -- --run`; exercise actual
workspace artifacts and rejected results without auto-loading them as accepted.

Each executed cell needs a prospective declaration, input/engine/checkpoint hashes,
commands and environment, progress, seeds/sample counts, raw/native outputs,
fit/restoration diagnostics, selection lineage, all refinement attempts, resource
measurements and independent final report. Retain cancellation/failure evidence;
resume requires validated completion records, not directory existence. Do not
delete rejected attempts, rename them accepted or overwrite historical bundles.

Start with G1/G2 and a small comparable local reconstruction panel after user
permission. Defer broad conditioning sweeps and GUI product implementation until
the backend interfaces actually work. At every gate distinguish implementation,
technical execution, scientific performance and production qualification.

## 13. Open decisions to freeze before execution

- Exact backend forks/commits/checkpoints and per-checkpoint license review.
- Concrete smoke/development/calibration/held-out structure and mask manifest.
- Safe per-runtime context sizes, memory/timeout/disk budgets on the 18 GB host.
- Admitted conditioning cells, candidate sample counts and optimizer budgets.
- Algorithm/tolerances for crossing alarms and classification of ambiguous paths.
- Final calibration grid, anchor windows, restraint strengths and reference scores.
- Minimal manual selection surface (private recipe/viewer first or GUI slice),
  followed by explicit approval for product integration.

Do not fill these decisions by copying old private Protpardelle ceilings or by
inferring authorization from a previous experiment. Preparation of this handoff
does not close any execution gate.

## 14. Executed bounded-study disposition

Continuation evidence and complete tables are in the
[comparison log](../research/mac-reconstruction-comparison-2026-11-11.md) and
[fixed validation policy](mac-reconstruction-fixed-validation.md). Retained artifact
catalog: `tmp/mac-reconstruction-comparison-2026-11-11/README.md`.

| Gate | Current disposition | Unclosed boundary |
|---|---|---|
| G0 | User authorized bounded autonomous execution; per-series declarations retained | No automatic new cohort/assets/budget expansion |
| G1 | Three isolated checkpoint-backed native-MPS runtimes executed under pressure-v2 | Host/case-specific feasibility, not universal size support |
| G2 | Exact mask/sequence/transport audits executed; stock Boltz explicitly proposal-only | Native OXT representation differs; no universal strict-inpainting claim |
| G3 | Declared admitted development cells completed, including terminal/pair and private MODELLER comparator | Interface cells refused original source data; historical 16-residue ML panels remain seed7-only |
| G4 | Bounded A/C/D and no-template/no-MSA E subset executed and independently audited | Other MSA/homolog cells lack frozen local assets; no remote upload permission |
| G5 | Compatible planner and manual locks/group lineage tested; all-generated joint assembly/refinement exercised | No accepted whole1HE8 assembly; new six-MODEL RFD3 donor visualization unapproved |
| G6 | 36-attempt anchor experiment retained; zero-motion frozen validation completed on two entries, all120 independently audited | Positive motion bounds not physically qualified; geometry does not imply accuracy; topology/trajectories unqualified |
| G7 | Evidence supports a private bounded-study report, not a public capability | Public CLI/Homology/GUI integration and production-default decisions require separate approval |

The frozen validation yielded 90/90 ML and 23/30 MODELLER final passes. Every
candidate retains exact observations/OXT and global L; seven MODELLER rejections
remain. E yielded 24/24 combined final passes (12 new calls, 12 reused CPU-matched
proposals), but only on one retrospective 1UBQ mask. These counts do not establish
general backend superiority or permit changing scientific gates.

Next decisions are distinct: manual review/approval of the new 1HE8 RFD3 MODELs;
a separately frozen complete interface/broader-family cohort if further sampling
is desired; local MSA/homolog asset provenance if those cells are desired; and a
separate product scope if integration is desired. None is silently inferred.
No public capability, P3, commit or push was enabled by this continuation.
