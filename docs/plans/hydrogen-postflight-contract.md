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
