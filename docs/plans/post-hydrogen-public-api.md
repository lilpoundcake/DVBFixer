# Private public-API post-H adapter and 1VCU continuation

## Authorization and documentation plan

After `1c8140a` the user instructed continued work until input is needed.
Continue privately on `research/post-hydrogen-public-api`, preserving variants
1/2 and their immutable evidence. No production/default registration, OpenMM
source modification, main merge or release. The post-H geometry rule remains
unchanged; this stage qualifies ownership/handoff and one independently activated
case, not a broader public capability.

Documentation order: this prospective contract, pure/private adapter tests,
frozen two-arm experiment, complete scalar outcomes/audit/validation, result
note and canonical handoff. No future links before creation; pending links none.
Validate links, source hashes and no public coordinate/RNG/State payloads.

## Public API adapter, before scientific execution

Use only public Modeller topology/positions and the public `addHydrogens` call.
The private adapter prepares a deep-copied working model, snapshots occurrence-
aware identities/graph before and after, runs the unchanged selective rule on
the raw returned geometry and commits topology/positions only after admission.
On API failure or guard/ownership refusal, the caller's model is untouched.
Nonempty residueTemplates mappings are outside this first adapter and refuse
before calling the API; do not silently lose object-keyed residue mappings.

Derive new H from exact before/after identity difference, cross-check element,
single heavy parent and canonical alias chemistry. Forbid new heavy atoms,
missing existing atoms, changed observed bond graphs and ambiguous occurrences.
Every caller-protected original H must be present before and after the call;
a stripped observed H is never relabelled new. The caller must supply that
lineage explicitly. Newly regenerated H from an upstream source with zero
observed H may be generated, but H already present at the current boundary are
always immutable to the guard. Preserve chain case, insertion codes and repeated
component occurrences using shared ExactAtomRef/variant knowledge.

Separate API transport from guard motion: original→raw API observation coordinates
must obey the already fixed half-larger-adjacent-float32-spacing +1e-15 nm bound
(Reference explicitly requested must be exact). Do not silently restore or move
observations to hide native rounding. Raw API return→guarded output preserves
every atom outside the exact new-H move mask byte-for-byte. Record both stages;
do not claim the whole native API is byte-preserving on single-precision devices.

Guard/ownership refusal uses a distinct geometry exception compatible with the
existing prepare HydrogenGeometryError escape and not ValueError/KeyError.
It must not trigger native FF fallback/retry. Preserve the raw returned model,
refusal reason and zero publication on failure; do not weaken thresholds.
Test atomic model commit/failure, return-value preservation, removed/protected H,
alias/graph/occurrence validation, transport bounds, no-op, nonempty templates
and RNG neutrality of the adapter (the underlying API may sample normally).

## Paired native 1VCU enrollment (freeze before execution)

Source: retained fixed-cohort case-08 prepared.pdb + prepared.dat, original tracked
1VCU.pdb, and the admitted static post-H case14 proposal at A/SER0/HA. Verify the
source really has zero H/D observations. Validate the .dat via DatRecord only.
Verify target/parent/heavy-neighbor identity and nm coordinates across the saved
return→prepared PDB handoff within the fixed PDB rounding bound (0.0005 Å plus
1e-12 Å/component); no heavy-coordinate modification or guessed atom matching.

Two arms, each exactly once: (1) unmodified prepared baseline, (2) prepared HA
proposal plus post-H guard at every subsequent successful public H-add call.
Copy both source sidecars byte-identically. The candidate prepared file may
change only the exact A/SER0/HA coordinate columns to the rounded static proposal;
all other bytes remain unchanged. This is a new controlled whole-complex
experiment, **not** an exact historical captured-System replay or new preparation.

Both invoke current minimize with --ff amber --ph 7 --no-solvent --keep-water
--atom-naming standard; normal 100/5 restraint tiers, max-iter1000, automatic
platform, automatic non-strict GAFF2 handling and all standard chirality/variant
policies remain unchanged. Do not force keep-H, disable variant stripping,
remove ligands, change parameters, or add H-only relaxation. A production
strip-and-splice fallback may execute normally, but is a visible failure of
whole-complex qualification, never counted as a full-system pass.

In each isolated research worker, a scoped public-method wrapper delegates to
the original addHydrogens through the adapter. Baseline observes/evaluates the
same candidate rule without applying it; guarded commits it. This is research
instrumentation only, not a patch installed into production or a global
LocalEnergyMinimizer interception. Seed Python random at 1 before each public
H call in both arms, then restore its caller RNG state. Preserve failures,
selected variants, actual FF provenance, raw/committed outputs and all masks.

Use an owned project ligand-cache directory. Baseline starts with an empty cache;
copy its resulting cache to a separate candidate cache before that arm. This
shares chemistry parameters, not target-H coordinates. Retain hashes before/
after and compare actual whole System particle/force parameters; any parameter
drift blocks paired-physics interpretation. No external user cache mutation.

Observe each Simulation.minimizeEnergy with a nonterminating reporter and save
exact Systems, atom/bond identity, requested/applied start coordinates, native
States/checkpoints and all callbacks privately. Do not intercept the global
minimizer or change tolerance/maxIterations/forces. Allow normal emergency
chirality recovery if production invokes it, but record every extra phase/
reflection; it is not evidence of an unassisted guard success.

Budgets: two workers, 1,800 s each / 3,600 s aggregate; at most eight public H
calls and eight whole-system phases per worker, 4,096 callbacks per phase.
Timeouts kill process groups. Durable per-path records and one final index;
retain partial results, no retries, seeds, alternate FFs or budget extension.
If candidate refuses or baseline fails, finish accounting without replacement.

## Independent audit and acceptance boundaries

Audit raw API→guarded exact masks/identity/radius/sign, original-source H lineage,
the removal/regeneration of the prepared HA, later production postflight and
final-System handoff. Report if normal H rebuilding eliminates the original
conflict before the guard, rather than claiming protection of a conflict that
did not occur. Check final outputs/DatRecord and complete ligand membership.
Re-evaluate saved whole-phase energies/mobile gradients, native stopping tests,
target HA spatial compatibility, Cα geometry and final heavy chirality using
shared domain primitives. Preserve starting defects, warnings and any emergency
repairs. Compare Systems across arms with restraint-anchor differences isolated
and reported, not silently ignored or retuned.

This is one independently activated structure and one declared sampling label,
not an accuracy/rate estimate or public readiness. If an unrelated scientific
blocker prevents whole-system construction, preserve it and stop that arm; do
not repair chemistry or fall back to an easier scientific case to get green.
