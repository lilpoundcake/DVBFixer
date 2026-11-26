# WP2R1: prospective joint cycle constructor, starts and comparison

Date: 2026-11-26; baseline `4565ac8`; persistent
`feature/chemical-reconstruction`; version remains 0.10.2. Freeze this protocol
and the [cohort/settings ledger](chemical-reconstruction-cycle-solver-cohort-v1.json)
before ANY new R1 initializer coordinates, optimizer call or scientific control.
Identity-only enrollment is allowed before freezing; it is not source admission.

## Scope and documentation plan

The user approved proceeding from WP2R0 to the prospective R1 protocol. This
document defines the first PRIVATE numerical experiment, not a shipped capability.
Its freeze does not mean an initializer, solver, DTO or validator is implemented.
Documentation scope: protocol and frozen cohort/settings ledger, handoff/task
discovery, and a local plan/quality note. Implementation validation/results require
a separate record after implementation; no future-result link is added now.

Reuse the standalone CycleConstructionProblem and its complete source-bound graph,
metrics, constraints, variables, blocks and retained context. R0 is unchanged.
The new constructor solves ALL generated variables together, including ring branches
and islands separated by fixed atoms. No PRO/CD rule or per-span independent fit.
The R1 initializer is NOT a Candidate, GeometryPlan or passing projection parent.
The old projector/staged-H gates and WP2P/Q results must remain unchanged.

No public API/CLI/GUI/Homology/DatRecord, publication/materialization, H placement,
default/dependency/backend replacement, version bump, main merge, release or
production/data-license approval. No licensed rotamer integration or new dependency.
No PDB structure output from this study. Retain private JSON, masked input bytes,
source proofs, raw/rounded additions, logs and every failure under
`tmp/chemical-reconstruction-cycle-solver-v1/` with exclusive writes.

## Explicit initializer provenance and refusal contract

Implement separate immutable CycleInitializerRecord, CycleSolvePlan,
CycleBranchResult and CycleConstructionResult roles (names are prospective).
Initializer inputs: ONLY the admitted CycleConstructionProblem. No external start
argument, old rejected/passing candidate, baseline, PDBFixer output, reference,
rotamer, fitted ideal pose, native score, learned embedding or random seed. Access
to authority ideal positions is for existing scalar angle/stereo/plane authority
only, never fitting/copying them into the initializer. Separate input preparation
and withheld-reference access from construction; construction receives the bound
problem, not complete original target coordinates.

Initializer provenance must bind algorithm revision, problem/source/observation/
graph/constraint digests, exact variable and anchor order, derived harmonic system,
length scale, branch index/sign/rotation, local scaffold and proposed additions.
It must also record refusal/failed starts with no fabricated coordinates or success.
All properties indicating scientific/geometry/MD approval remain false at this role.

Algorithm `boundary-harmonic-orbit-v1`:

1. Construct the right-handed observed frame from the three R0 anchors `(a,b,c)`:
   origin `a`, e1 = unit(b-a), e3 = unit((b-a) cross (c-a)), e2 = e3 cross e1.
   Local coordinates are dot products with e1/e2/e3. Fail on nonfinite/degenerate
   arithmetic; do not replace anchors or choose a different frame.
2. In exact R0 variable order, assemble the heavy-graph Dirichlet Laplacian with
   edge weights `1 / ideal_bond_length_angstrom**2`. Generated/generated edges
   contribute to both diagonals and negative off-diagonals; generated/fixed-target
   edges contribute diagonal weight and weighted observed local position to RHS.
   Do not use external/context coordinates as fictitious bonded boundary atoms.
   Solve `L * h = rhs` in float64 with NumPy's direct linear solve, all variables
   simultaneously. Every generated block must reach an observed boundary per R0.
   A singular/nonfinite solve or maximum absolute equation residual >1e-10 A^-1
   refuses initialization. No pseudoinverse, jitter, ideal overlay or alternate start.
3. Let ell be the ordinary median of the lengths of ALL generated-involving target
   heavy bonds (even count: arithmetic mean of central two). It must be finite and
   positive. For zero-based variable j of n, theta = 2*pi*(j+1)/(n+1),
   v_j = (cos(theta), sin(theta), (-1)**j / 2).
4. Declare EIGHT starts in order `(sign,k)` = `(+,0)..(+,3),(-,0)..(-,3)`.
   phi = k*pi/2; rotate v about local e1 by phi, using
   `(x, y*cos(phi)-z*sin(phi), y*sin(phi)+z*cos(phi))`.
   Start j = h_j + sign*ell/2*rotated(v_j). Only then transform to source Cartesian
   coordinates and create generated Observation records with authoritative element,
   occupancy 1 and angstrom units. Preserve every measured record unchanged.
5. This scaffold is artificial, not claimed physically valid. Measure and retain
   initial geometry/findings but do NOT require a passing parent to attempt R1.
   Do not select/deduplicate starts, stop after first success, retry failed branches
   with another direction, or treat the eight starts as eight independent proteins.
   Any nonfinite/degenerate start/residual fails its labelled branch explicitly.

This initializer is a bounded hypothesis, not an established superior distance-
geometry method. Its perturbation is not a conformer library. No selected pose or
ring pucker is taken from either comparator or withheld source atoms.

## Frozen numerical plan and residual definition

Solve Cartesian displacements in the observed local frame relative to EACH own
initializer. There are exactly 3*n variables, no observed-coordinate variables.
Every coordinate displacement bound is +/-8/sqrt(3) A, giving a raw Euclidean
per-atom cap 8 A. Independently enforce the SAME Euclidean cap after PDB rounding;
rounding that crosses it rejects, not clips. This is an artificial start-relative
search budget, NOT evidence of small physical motion. It is distinct from the old
projector's unchanged 2 A passing-parent cap. Measured atom motion is exactly zero.

One SciPy `least_squares` call per initialized branch: method trf, jac 2-point,
tr_solver exact, x_scale 1.0, diff_step 1e-6, loss linear, ftol/xtol/gtol 1e-10,
max_nfev 600, no callback restarts or robust-loss reweighting. Count ALL residual
calls separately, including finite-difference and audit/objective calls outside
SciPy nfev; optimizer residual-call ceiling is 50,000 per branch. Ceiling breach
is a retained budget failure. Record package/runtime/BLAS/platform versions and
all stated/effective settings. Reproducibility is scoped to that environment,
not promised bit-identical across CPU/library implementations.

Objective: sum of squared dimensionless residuals, all weights 1, in this order:

- All R0 heavy bond/angle metrics, INCLUDING fixed-only constants, in coverage
  order: `(measured-target)/tolerance` (existing 0.10 A and 12 degrees).
- All source then extra constraints in original order. Distance/ring-closure/angle:
  `(value-(lower+upper)/2)/((upper-lower)/2)` in the declared units. Stereo:
  `max(0, minimum_volume-orientation*ordered_volume)/minimum_volume`, using existing
  ordered-frame/ligand semantics, NOT independently inferred CIP or a torsion wall.
  Plane: one signed distance/upper per point after its first three ordered atoms,
  with existing normal and converted distance units. R1 explicitly refuses a
  nonzero lower plane bound, which this profile does not model. Exclusion:
  `max(0,lower-value)/lower + max(0,value-upper)/upper`; both declared bounds are
  respected, unlike silently dropping a finite upper bound. Unsupported periodic
  torsion semantics remain refused. No constraint is silently omitted or reordered.
- Eligible generated/generated and generated/observed severe-separation pairs,
  canonical identity order: `max(0,0.7-distance_angstrom)/0.7`. Exclude only self,
  1-2 and 1-3 graph neighbours from resolved chemistry and declared covalent links;
  never invent opaque-context internal bonds or missing chemical authority.

Include complete full-complex context. For an observed pair, omission from the
numerical vector is allowed ONLY if its distance from that branch's initializer
is strictly >8+0.7 A, proving it cannot enter the cutoff anywhere in the movement
ball. Always retain all eligible generated/generated pairs. Enumerate and bind
this fixed per-branch pair list BEFORE its solve, not a changing-length neighbour
list at each evaluation. Keep the proof/count of omitted far pairs; independently
audit final geometry against ALL retained observations with no spatial crop.
Max 65,536 numerical residual scalars including context; larger cases refuse, never
truncate. R0's separate 512 chemical-term planning limit is not changed.

No weak positional regularizer, energy/minimization/force field, native objective,
chirality reflection, graph/variant/mask/threshold adjustment or candidate-derived
conditioning. A degenerate angle/plane/stereo evaluation, nonfinite residual or
linear-algebra/optimizer numerical failure fails that branch with retained evidence;
do not replace it by a favourable zero residual or an undocumented finite penalty.
Missing NumPy/SciPy/OpenMM/required source resources is environmental unavailable,
not chemical rejection. Unexpected programming errors abort with evidence, not
blanket convert to scientific refusal. No fallback constructor.

## Independent acceptance, selection, traces and saved readback

Retain every initialized branch even if no geometry passes. Record solver returned
x/raw additions, raw residuals, initial/final sums of squares, nfev/njev/call counts,
success/status/message/exception, movement, rounding and all validation findings.
If the solver cannot return a complete finite x, retain the initializer and last
finite evaluation as explicitly PARTIAL diagnostics (never select/score as a full
proposal). Solver success alone does not approve a structure or source chemistry.

Independent measurement must NOT call/reuse the numerical residual objective as
its sole validator. Reuse existing exact-identity/contracts and scalar measurement
owners, comparing all R0 metrics, all declared constraints, complete generated mask,
authoritative elements/units/occupancy and full-context severe separation separately.
No fabricated GeometryPlan or old strategy merely to gain access to its validator.
Check raw and float(format(coord,'.3f')) rounded additions independently, retaining
both. Observations, source records/links and masks are immutable and excluded from
generated additions. Byte/hash proofs bind original/prepared/masked source and
unmodified observed record inventory; no input or public output writer is added.

Use the existing `assert_all_l` detector on a metadata-only topology containing
all observed plus generated protein atoms, for BOTH raw and rounded states. Preserve
exact chain/resid/icode/occurrence identities through the adapter. Reuse/factor the
existing comparator's pure global-L readback rather than a second chirality rule;
do not call PDBFixer, rebuild, Modeller, reflect or regularize to perform this check.
Pre-existing source D also refuses a proposed whole-complex pass. Retain findings,
not a hidden input repair. Chemical stereo constraints remain additional checks.

A branch passes only with solver success, complete finite additions, both independent
geometry/global-L audits, exact transport/rounding and raw/rounded movement caps.
All eight labelled starts are attempted unless initialization/source admission is
unavailable/refused or an unexpected error aborts. Select among passing branches
lexicographically by `(raw final squared residual, raw max movement, branch index)`.
If none pass, return rejected/no selection with all attempts, not infeasible-chemistry.
Do not choose a lower held-out RMSD branch or score partial additions.

Distinct DTOs bind problem/plan/start digests and all attempts/selection/refusals;
no new Strategy/GeometryPlan variant. If strict private transport is extended,
register standalone roots only and preserve old parser limits/variants. Canonical
digests/JSON, not Python numeric equality, bind cached records. Planned, initialized,
converged, geometry-passed, native quality, source coverage and MD readiness remain
separate states. Even a selected R1 result is private heavy declared-geometry-only,
not an H parent, qualified general reconstruction or production/public capability.

Saved readback must reopen pinned sources/authorities and reconstruct exact masks,
problem/plan/pair lists and enrollment; measure saved starts/outputs/movement/rounding/
global L/selection/native scores without any coordinate constructor, optimizer or
PDBFixer. For harmonic scaffold, check the saved equation and offset formula with
the declared equation tolerance, retaining the distinction from bit-exact replay.
A SEPARATE initializer-only deterministic replay (optimizer/baseline/reference
disabled) must reproduce all start records byte-for-byte on the same environment.
That replay is not an independent scientific sample or the constructor-free audit.
Test tampering of source/state/mask/constraints/coverage/start/pair proof/status/
solver metadata/raw/rounded/global-L/selection, including numeric representation.

## Independently declared 30-cell cohort and comparator arms

The ledger enrolls three tracked systems, NOT 30 independent proteins. All were
used in older research and are known, not blinded/external/unseen qualification.
Identity-only discovery: first file-order ATOM chain, then first interior target
of the declared residue name via the existing selector, no geometry/native filter.
No fallback chain/site/residue, altloc choice, occupancy rewrite, context crop or
natural-missing-mask enlargement. Source/proof/planning failures stay in denominator.

- Primary: all seven nonempty subsets of PRO CB/CG/CD on 1TIM A/2, 1CRN A/5,
  2CI2 I/25: 21 cells. Includes all nine old WP2P cells plus noncontiguous subsets.
  Those nine are retrospective development regressions, NOT new independent data.
- Secondary: complete heavy PHE/TYR/TRP sidechains at first interior targets in the
  SAME first chains: nine cells, including the metadata-absent 1CRN TRP slot.
  Keep that absence without replacement; other cells have no source-admission
  guarantee merely because an identity exists. These test ring/fused/branch scope,
  not all-chemistry coverage. No HIS microspecies is silently selected.

Exact literal masks/state/identity/path/SHA256 are frozen in the ledger. Use existing
source preparation/observed probes/explicit peptide authority and opaque/fractional
proof semantics. PRO uses unchanged R0 incoming-type/byte proof and v2 port extras;
PHE/TYR/TRP use their declared state (same parent name) with original authority
constraints only, no fictitious PRO port. Real peptide partners remain observed;
ALA chemical witnesses are not actual-neighbour/full-polymer qualification.

Three heavy arms, independently attempted for each source-admitted cell:

1. New eight-start joint R1 constructor (only this arm uses R1 starts/solver).
2. Existing deterministic original graph-internal proposal, original request,
   no projector, search or H; retain refusals/incomplete/raw/rounded results.
3. Existing separate seeded-PDBFixer baseline helper, one call per cell, internal
   seeds 1..5, OpenMM Reference, exact source-bound masked input and target atoms,
   no H/reflection/minimize or other residue completion. Retain temporary topology,
   fixed movement/resource checks and global L per unchanged comparator contract.

Do NOT gate the new constructor on either comparator's geometry or copy their
coordinates as starts. Attempt all three available arms even if one fails. Repeated
old-nine comparator calls here are labelled repeats, not new independent evidence;
do not overwrite or rerun the historical WP2P/Q driver to change its result.
Read/verify original retained old evidence separately without scientific reruns.

After ALL available proposals/branches for a cell return, extract the exact masked
native references for measurement ONLY. Measure all complete raw/rounded branches
and selected result, plus both comparator outputs, under the SAME effective
original+extra request and frozen heavy policy. Use measurement-only DTOs, not
constructor relabelling. Report exact-name fixed-frame RMSD/max/count, no fitting,
symmetry/name permutation or reference-based branch selection. Common geometry
includes generated-context unique severe-overlap pairs below 0.7 A/minimum eligible
separation, not a vdW/MolProbity/energetic score. Native geometry is a separate
post-hoc diagnostic; do not remove a cell or relax bounds because it fails them.

No required pass rate/RMSD win count. Report all 30 cells, phase/arm/refusal/reason,
start/solver/geometry states, selected and complete-rejected native scores, explicit
paired denominators and native worsenings. Never pool branches, nested masks,
raw/rounded, replays/repeats or same-system sites as independent trials. No H or
native-H/naming/full-polymer superiority claim. Incomplete branches have no partial
RMSD. Preserve failures as failures; no retuning on this cohort.

## Implementation and execution gates

After this freeze, implement separate start/solver/independent audit/trace modules
with narrow synthetic and generic ring tests BEFORE tracked-source numerical runs.
Include full-sidechain, fused/branched rings, noncontiguous subsets, linked/opaque/
fractional context, fixed conflicts, exact identities/case/insertion/occurrence,
degenerate/budget/environment/refusal controls and old projector/H gate regressions.
Metadata absence/readiness is not a geometry pass. Qualification results do not
exist at protocol freeze; do not add prospective counts as completed experiments.

Before scientific execution, record source-code/test/environment/protocol/ledger
hashes and test outcomes without inspecting held-out reference scores. No new
R1 coordinate runs until the separate initializer/solver/audit/trace slice and its
nonleakage tests are ready. Run focused/boundary/full-fast, Ruff, scoped/source-bound
typing, maps/generated/version/fixture/link/whitespace checks as appropriate for
implementation. Keep old 814 ambient errors, GUI dependency state and WP2F timeout
visible; these are historical baselines, not assumed current passes.

If this declared algorithm fails, retain it and its outcomes; any change in starts,
weights, solver, caps, masks, states, thresholds, source selection or retries needs
a prospective amendment before another numerical cohort, never retroactive tuning.
H-parent integration/publication and broader native/retained-complex scientific
qualification remain later separate gates. Production/data-license/public/default/
dependency replacement gates remain closed regardless of this bounded study.
