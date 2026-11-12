# Authorized O–H postflight implementation and qualification

## Authorization and scope

The user approved development and qualification of a behavior-changing H
postcheck on a **separate branch**, without changing the native force field or
weakening existing tests. Branch: `fix/hydrogen-postflight-contract`, starting
from `347efe6`. This is not authorization to merge/release an unqualified change.
The [SER study](../research/ser-hg-capture.md) and lost original failure remain
historical evidence; do not rewrite their sources, hashes or expected outcomes.

The [frozen qualification plan](hydrogen-postflight-qualification-plan.json)
binds the baseline, complete local-graph grid and six saved native inputs. Its
protocol-note input is a preserved pre-implementation copy, not this live handoff.

Documentation scope: this implementation/qualification record, the shared
hydrogen invariant/task map, preparation context, and a handoff update. Validate
all links and source references. One owner coordinates detector, repair and
callers. No GUI/CLI flags, version bump, naming, `.dat` or parameterization change.

## Frozen behavioral proposal (before implementation)

- Detect an O–H distance **at or above 1.2 Å** as `parent_distance`, independently
  of the existing 1.6× gross-distance heuristic. This upper sanity limit is
  oxygen-specific, not a blanket absolute limit for C/N/S/P bonds. Existing
  covalent target lengths, coincident and sibling-H criteria stay unchanged.
- Keep already-valid H byte-for-byte unchanged. For newly detected moderate
  O–H elongation (below the old 1.552 Å gross cutoff), first shorten **only that
  H along its existing bond direction** to the existing 0.97 Å target. Do not
  replace a valid hydroxyl/phenol torsion with an arbitrary new direction.
- Accept radial shortening only if it introduces no same-residue coincidence
  or sibling-H clash under the existing guards. Otherwise use the existing
  tetrahedral repair path. Gross protein/ligand O–H collapse keeps that path.
- For oxygen with no heavy neighbor (water), radial shortening is also allowed
  for gross elongation, when existing H directions and sibling checks permit it;
  no arbitrary water orientation is invented.
- If an affected O–H remains at/above the limit after bounded repair, raise a
  typed geometry error rather than silently publishing it or resampling H.
  Prepare must not reinterpret this error as an FF-template fallback. Existing
  fallback for genuine template/ligand errors remains intact.
- Heavy atoms, atom/residue identity, topology, insertion codes, valid sibling H,
  nonoxygen distance thresholds and force-field/minimization settings stay fixed.
  This is not a solution to the separate THR2 HB handedness mechanism.

## Predeclared qualification

1. Preserve the baseline geometry module bytes from `347efe6` and its digest in
   a project-local private evidence directory before editing production.
2. Keep the original SER regression assertions unchanged. Add a deterministic
   noncoincident 1.245813 Å O–H regression that fails against the baseline.
3. Run the complete Cartesian grid of **eight chemistry contexts** (SER, THR,
   TYR, ASH, GLH, water, ligand hydroxyl and sugar hydroxyl) × **six distances**
   (0.97, 1.199, 1.2, 1.201, 1.2458133889150464, 1.7 Å): **48 local graphs**.
   Water includes a well-placed sibling H. These are deterministic local graph
   controls, not 48 independent experimental structures or checkpoint samples.
4. Preserve C/N/S/P behavior in nominal/stretched controls, explicit chain-case
   and insertion-code controls, nearby-atom fallback, degenerate/unrepairable
   failure and fixed-heavy/valid-H identity checks.
5. Compare old/new repair on **all six saved native SER H returns**, without
   new H sampling or optimization. Keep every changed/unchanged/error outcome.
6. Pin old research guard-gap assertions to their recorded historical evidence;
   do not make them pretend current production still has the old bug. Preserve
   their old expected values and add separate executable tests for the fix.
7. Run focused tests, strict mypy, Ruff, generated specs, agent-doc validation
   and one full fast suite. Retain failures without weakening assertions or
   adding broad skips. The old 3,082-pass/one-failure run remains a failure.

Output PDB quantization can perturb a distance just below 1.2 Å across that
boundary; this in-memory upper guard is not a proof about every rounded boundary
case. Qualification must report that limitation rather than silently move the
test's threshold or claim bit-identical CPU/GPU sampling.

## Candidate v1: retained unsuccessful qualification

Prospective checkpoint: `eeae9dc`. The [first matched qualification](../research/hydrogen-postflight-qualification-evidence.json)
passed its narrow distance contract on all **48 local graphs + six saved native
states** (baseline passes 29/54). Heavy coordinates/identity and valid siblings
were preserved; all six old gross SER repair outputs were exactly unchanged.
This was insufficient for the full existing regression contract.

The first full fast run retained **3,144 passed / two failed**, three optional
skips and 33 slow deselections:

- The AST-based source inventory still points at historical source hashes and
  line numbers. Preserve those historical inventories and add a new current
  source amendment rather than overwrite old scientific evidence.
- The original unchanged SER test found OG–HG **0.969632 Å** but HG–OC1
  **0.718741 Å**, failing its existing **>1.0 Å** contact assertion. A correct
  covalent length alone is not enough. The original failed final PDB/Dat and
  candidate source are retained in-project. Its native intermediate path was
  not captured; do not claim an exact unobserved sequence from the final PDB.

[Candidate-v1 failure record](../research/hydrogen-postflight-candidate-v1-validation.json)
binds those sources, outputs and logs. This is an unsuccessful implementation
candidate, not a green run, and is not authorized for merge/release. Additional
validation requires an explicit pre-execution amendment, not a silent repeat.

## Prospective amendment v2: include the existing contact assertion

This amendment is frozen before candidate-v2 editing or execution. It is prompted
by the retained v1 failure, not by a search for favorable random seeds.
The [v2 plan](hydrogen-postflight-qualification-plan-v2.json) binds its exact
inputs, retained failed candidate, canonicalized failure and contact grid.

- Keep the oxygen-specific 1.2 Å upper bound and every original assertion.
- Also detect an oxygen-bound H **within or at 1.0 Å of a same-residue nonparent
  heavy atom**, with a distinct `oh_nonparent_contact` finding and measured
  contact distance. This generalizes the existing HG–OXT assertion without
  hard-coding SER or OXT. The 0.5 Å general coincidence guard and all nonoxygen
  policies stay unchanged.
- Require a radially shortened H to clear such a neighbor by **1.01 Å**: the
  existing 1.0 Å predicate plus 0.01 Å serialization clearance (larger than the
  maximum 0.001733 Å relative-distance perturbation from two PDB coordinate
  roundings). This is an admission check for a new candidate, not a widened test.
- Fall back to the existing tetrahedral placement when radial admission fails.
  If bounded repair still violates either O–H predicate, raise the typed error;
  do not search orientations, change native forces or resample H. Diagnose must
  report the nonparent distance, not falsely describe a correct O–H bond length.
- Preserve and bind the actual v1 failed PDB/Dat, canonicalize names only using
  the shared naming service, and require that its unchanged geometry is either
  repaired successfully or blocks qualification. The real retained SER failure
  must be repaired, not classified as an allowed refusal.
- Re-run the fixed 48+6 distance/saved-state comparison under the new candidate,
  add **12 deterministic contact graphs** (O–H 0.97/1.245813 Å × contacts
  0.719/0.999/1.0/1.01/1.1/1.5 Å), and the **one actual retained v1 failure**:
  **67 total records**. Intentionally blocked toy geometries may refuse only
  with the typed error when the independent pre-state violates a predicate.
  Valid contact inputs may not refuse or move. Keep corrected and refused
  outcomes distinct; do not call them independent structural samples.
- Add a separate current AST source inventory amendment, proving unchanged
  PDBFixer import/use semantics and totals apart from hashes/line numbers.
  Never overwrite either historical research inventory or re-run its probe.
- Run targeted checks and one new full fast suite **after actual code changes**,
  with a new project-local basetemp. Retain v1's two failures and all v2 outcomes;
  no retries of an unchanged candidate merely to obtain green.

## Candidate v2: matched qualification results

Amendment checkpoint **`693f54b`** precedes the v2 implementation and execution.
The [v2 evidence](../research/hydrogen-postflight-qualification-v2-evidence.json)
reports all **67 records**:

- Original **48 local graphs + six native SER states** satisfy both independent
  distance and contact checks, with no candidate errors. All heavy coordinates
  and identities remain fixed; valid graph inputs and unselected siblings stay
  unchanged.
- Of 12 intentionally stressed contact graphs, eight finish safely and four
  explicitly refuse with `HydrogenGeometryError`. Those four blocked toy
  geometries were invalid before repair; they are not discarded failures or
  successfully reconstructed structures.
- The **actual retained candidate-v1 SER failure** is corrected, not refused:
  OG–HG **0.969632 → 0.970000 Å** and HG–OXT **0.718741 → 1.788592 Å**.
  Historical baseline and failed v1 both leave that final state uncorrected.
  The other SER hydroxyl, heavy coordinates and atom identity stay unchanged.
- Independent qualification readback is byte-identical. Neither qualification
  invokes H sampling, optimization, FF changes or a fresh original-trajectory
  reconstruction. PDB-name normalization uses the existing shared service and
  preserves the observed coordinate columns exactly.

The [current source inventory amendment](../research/pdbfixer-usage-inventory-hydrogen-postflight.json)
updates geometry/prepare/diagnostic source hashes and line numbers while proving
unchanged PDBFixer import/use semantics and totals. Both older inventory files
and their template evidence remain unchanged.

The expanded focused check passed **94 tests**, including the original unchanged
SER assertions, nonoxygen nominal/moderate/gross controls, fail-closed preparation,
identity and diagnostic-contact cases. The single full candidate-v2 fast run
passed **3,174 tests**, with three optional Torch/Biotite skips and 33 slow
deselections. The [validation record](../research/hydrogen-postflight-validation.json)
binds logs/XML and final source hashes; candidate-v1's two failures remain retained.
Ruff, strict-target mypy, generated CLI/GUI specs and agent-documentation checks
also passed. No original SER assertion, force-field parameter or minimizer setting
was changed. All pytest failure/success artifacts now use project-local basetemps.

### Limits and publication gate

These are deterministic local graph and saved-state controls, not 67 independent
physical structures. They establish the stated coarse O–H contract, not a native
minimizer convergence guarantee, correct protonation or torsion, all-atom steric
validity, or a fix for the separate THR2 HB handedness issue. Short O–H bonds
outside the unchanged coincidence guard remain outside this upper-bound patch;
for example, a retained SER125 O–H around 0.684 Å is not silently normalized here.
Near-boundary PDB quantization remains a separate limitation.

Implementation stays on `fix/hydrogen-postflight-contract`. Merge/release remains
a user decision after review of qualification and complete validation, not an
automatic consequence of the research authorization.

## Real-input follow-up

The next bounded stage is the [prospective real-input protocol](oh-native-qualification.md):
six tracked structures, two FF requests each, exact matched saved-H replay,
explicit actual-route coverage and retained failures. This is not a main merge.

[Real-input results](../research/oh-native-qualification.md): all twelve attempts
executed once, eleven preparations completed, and all twelve postflight states
replay exactly with no new refusals or heavy/valid-H movement. Twelve additional
O–H violations were corrected, all on forcefield-free fallback. The unchanged
GLYCAM/CHARMM fallback case fails naming with duplicate sugar H1 created before
postflight. Full-cohort qualification remains false; no native rerun or main
merge. The separate offline chain-run adapter fixes only auditor transport.
