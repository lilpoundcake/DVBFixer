# Mac reconstruction: controlled comparison continuation

Date: 2026-11-11. User authorized continuation after reviewing the incomplete
handoff gates. This is a staged private study, not public capability approval.

## Execution order and bounded first slice

1. Freeze comparison rules and add the missing full MODELLER proposal baseline
   on the five already inspected masks below. This is retrospective development
   evidence, **not** a new independent accuracy benchmark.
2. Prepare a prospective, source-bound expansion covering termini, interacting
   gaps and protein interface context. Preserve admission refusals. Freeze inputs,
   shapes and per-backend support before proposing the next sampling manifest.
3. Execute matched development proposals and common repair; assess context,
   backbone/all-heavy information and independent/joint generation in separate
   bounded comparisons. cc94 participates only where its current interface admits
   the exact task. Unsupported cells must not be silently cropped or relabeled.
4. Calibrate necessary-only anchor movement on reference cases separately; then
   freeze the policy and execute independent validation. Do not tune on 1HE8.
5. Consider product integration only after these gates; no GUI implementation yet.

### Slice 1: declared MODELLER baseline

| Projection | Hidden residues | Seeds | Full proposals |
|---|---|---|---:|
| 1UBQ A1–72 | A10–13 | 7, 42, 101, 102, 103, 104 | 6 |
| 1A3N A1–59 | A30–37 | same six | 6 |
| 1FR2 A3–85 | A24–31 | same six | 6 |
| 1A3N A1–59 | A30–45 | same six | 6 |
| 1FR2 A3–85 | A24–39 | same six | 6 |

Thirty serial CPU jobs; no new ML checkpoint sampling in this first slice.
Each job allows 1800 s, aggregate process RSS 8 GiB, pressure-v2 host guards and
a 10 GiB study evidence cap. Worst-case serial wall budget is 15 h, not a runtime
estimate. Existing MPS/CPU runtimes remain unchanged. No remote sequence services.
All jobs reserve the existing study process lock. Outputs remain inside project.

The private seeded adapter uses MODELLER 10.8 full LoopModel construction, not
the earlier complete-proposal conjugate-gradient repair. It reuses production
chirality restraints, initial whole-model optimization and gap-only loop MD.
Explicit seed is `-development_seed`. One initial model and **one** fast loop
derivative are declared per job (production default is two). There is no
best-of-many selection, fallback to an initial model, empty-selection all-model
fallback or post-hoc raw chirality reflection. The exact PIR mask is transported
from the established request, not re-inferred by another alignment.

## Frozen assessment rules

- Keep unmodified native, identity-mapped native, observed-only fitted native,
  source-restored and final repaired outputs distinct. Report native drift and
  restoration-induced junction/geometry changes; MODELLER is not falsely called
  strict source-preserving because its final observations were restored.
- The generator and identity transport read only source/request. Held-out native
  coordinates are read only by post-generation evaluation. This procedural
  separation does not claim an OS sandbox or training-set disjointness.
- Require exact requested chemistry/heavy-atom inventory, finite coordinates,
  independently verified source custody and global L-chirality. Native terminal
  OXT additions may be explicitly reported and excluded during identity transport;
  the original native file remains retained. No other extra atoms are accepted.
- Use existing versioned assessment thresholds unchanged. Copy/bind scientific
  implementation hashes in the execution plan. Geometry failures stay failures;
  errors in raw data/transport/resources remain separate technical failures.
- Clean proposals are not repaired. Geometry-failing source-restored proposals
  receive the same existing generated-only CPU repair as the earlier reference
  study: 1000 iterations, seed 7, zero restarts, no geometry preconditioning or
  anchor movement. All generated regions are mobile together; no solvent.
- Report per-candidate native/source-restored/final geometry, junctions, heavy and
  backbone held-out RMSD, observed drift, repair movement, runtime and process RSS.
  Do not pool seeds as independent structure-level replicates. Report all six
  seeds with failures in denominators, plus per-case medians and ranges.
- Do not choose candidates or tune protocols using held-out RMSD. The first slice
  evaluates the full six-candidate distribution; it does not establish an automatic
  one-candidate selector. Model confidence is not independent acceptance evidence.
- The previous longer-mask ML arms have only seed 7. Their comparison with six
  MODELLER seeds is descriptive and unmatched until the ML cohorts are extended.
- Segment/surface/ring diagnostics remain experimental. No topology or optimizer
  trajectory safety is inferred from clean endpoint alarms.

Validation seeds 211/223/227/229/233 remain unused. 1HE8 remains a known stress
case, not an independent validation target. No arbitrary anchor-limit changes,
production defaults, public schema changes, P3, commits or pushes are authorized
by this continuation.

## Evidence and implementation

- `scripts/mac_reconstruction_comparison.py`: reference-blind alignment/identity
  contracts, stage distinctions and disjoint seed sets.
- `scripts/run_mac_modeller_proposal.py`: one full seeded native proposal.
- `scripts/run_mac_modeller_comparison.py`: hash-bound declaration, serial guards,
  post-generation audit/common repair and retained failures.
- `tmp/mac-reconstruction-comparison-2026-11-11/`: new execution evidence; historical
  `tmp/mac-reconstruction-2026-11-11/` remains unchanged.

Execution results are recorded separately; this protocol is not a completion claim.
