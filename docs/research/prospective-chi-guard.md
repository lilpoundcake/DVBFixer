# Prospective new-input comparison: mean/refinement does not beat the simple baseline

Status: **bounded new-input research, negative for the library-mean/refinement
advantage hypothesis**, not a production capability or broad scientific
qualification. Continues the [actual-endpoint guard](chi-endpoint-guard.md).
Installed PyRosetta remains [research-only with mandatory production gates](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03).

## Enrollment before download, choices before native metrics

Before fetching any new structure, `prospective_chi_guard_cohort.py` froze:

- **1ENH** (engrailed homeodomain), **1CRN** (crambin), **1QYS** (Top7);
- first coordinate MODEL, first encountered ATOM chain, first two VAL and first
  two THR identities in file order, whole-sidechain masks; no substitutions;
- unchanged observation/geometry/score/mean/χ kernels and complete-pool selector;
- six comparison arms: stock seeded PDBFixer Reference, synthetic baseline fixed
  score, mean raw-family prior, mean fixed score, unchanged mean-C, bounded mean-guard;
- **any execution blocker makes a pool incomplete; any guard budget exhaustion
  blocks the entire guarded pool**, never select only the successful starts;
- no new angles, tolerance/restart grid, weights or preferred-profile selection.

This is prospective **input/policy/core-kernel enrollment**, not a fully
preimplemented, blinded, randomized or training-independent benchmark. Transport
adapters were completed after download without changing those scientific kernels
or selection rules. Top7 is Rosetta-designed, not independent evidence about
Rosetta design/training; its projected requests were not executed here.

All six arms' choices committed in the
[native-free ledger](prospective-chi-guard-evidence-selection.json) **before native
accuracy reporting**. A separate plan then froze **all 15 unordered arm pairs**
before the accuracy report existed. No favorable pair or native-error subset is
chosen for the [matched summaries](prospective-chi-guard-pairs.json).

## Scope, refusals and implementation evidence

All three official RCSB downloads returned HTTP 200; original bytes, URLs, UTC
times and SHA-256 hashes remain inside the project. The frozen rule yields:

| Input | Enrolled masks | Outcome |
|---|---:|---|
| 1ENH | 2 THR, 0 VAL | Both admitted; two-VAL shortfall retained |
| 1CRN | 2 THR, 2 VAL | Four admitted; THR1 lacks the required observed φ/ψ context for mean queries |
| 1QYS | 2 THR, 2 VAL | All four share the preserved missing/mismatched LINK/SSBOND-endpoint refusal |

Thus **10 masks**, two count shortfalls, **four correlated shared refusals**, six
admitted inputs and **five mean pools**. The Top7 refusal concerns the frozen
projection retaining touching metadata while excluding an endpoint: do not call
the original deposition invalid or silently sanitize the projection into admission.
THR1's absent peptide neighbor is not filled with an invented library angle.
Only **two proteins** contribute executed scientific cases.

Current PDBFixer/Reference preflight is separately pinned (PDBFixer 1.12.0,
OpenMM 8.5.2). It establishes required-resource usability here, not the cause or
repair of earlier damaged installations. No environment was modified. The initial
library-query executable failed **before any query** because it eagerly imported
an OpenMM-dependent adapter into the isolated PyRosetta environment. Original
source/log retained; a separately frozen stdlib/PyRosetta-only entry point fixes
dependency transport, not scientific settings. No dependencies were installed.

The installed binary remains the previously hashed `StpDwn_0-0-0` library at
PyRosetta revision `1a56185`; no equivalence to original recommended upstream
tables or normalized probabilities is assumed. Five complete three-well queries
and their periodic controls passed; individual means/weights/sigmas stay private.

Independent readback verifies **30 NCAC stage snapshots**, source-bound proper
fits/rotations, every candidate materialization and **45 exact physical/topology
score imports** (including the auxiliary relative-NCAC control). Scorer extras,
virtual coordinates and non-target H match within and across fixed-score arms.
Synthetic proposals pass the common gate **14/18**; means **15/15**. Both C and
guard return **15/15** final common passes, unchanged protected coordinates and
no target/global chirality findings. The guard uses **16 calls for 15 starts**,
achieves the actual gradient criterion on all 15 and exhausts no budget.

## Accuracy: report unequal coverage and matched comparisons separately

[Accuracy evidence](prospective-chi-guard-evidence.json), rounded generated-heavy
atom RMSD, Å. These rows have **different selected denominators**, so their
unmatched medians are not direct superiority comparisons:

| Arm | Selected / enrolled | Exact-name median |
|---|---:|---:|
| Stock seeded PDBFixer Reference | 2 / 10 | 0.167449 |
| Synthetic baseline + fixed score | 6 / 10 | 0.191971 |
| Mean raw-family prior | 5 / 10 | 0.283035 |
| Mean + fixed score | 5 / 10 | 0.278574 |
| Mean + C | 5 / 10 | 0.281341 |
| Mean + bounded endpoint guard | 5 / 10 | 0.281341 |

On the **same five mean-supported masks**, synthetic baseline median is
**0.135648 Å**, versus **0.278574 Å** for mean fixed score and **0.281341 Å** for
C/guard. The declared refinement and corrected numerical stopping do **not**
rescue selected-output accuracy. C and guard select the same indices and produce
**byte-identical selected rounded heavy PDBs in 5/5 cases**; the guard's additional
call concerns a nonselected start. No numerical mode is adopted as preferred.

Conversely, on the **two stock-supported masks**, stock median is **0.167449 Å**
versus synthetic **0.305992 Å**. This counterexample is retained: synthetic's
greater coverage does not establish template accuracy inferiority. Stock is the
bare shared seeded PDBFixer comparator, **not full production legacy preparation
plus restrained minimization or ZBS**.

All 15 matched arm pairs are reported. Six admitted native controls pass the
unchanged common gate. Generated-only VAL symmetry gives the same cohort medians
here; no native-based renaming or candidate selection. Sample sizes remain small,
correlated and deliberately selected; no significance, population, training-
independence or broad chemistry/production qualification claim follows.

## Retained evidence and reproduction

All paths below are under `tmp/pdbfixer-research/`:
`prospective-chi-guard-enrollment-01/` (pre-download policy/source hashes),
`prospective-chi-guard-inputs-01/` (official original bytes/provenance),
`prospective-chi-guard-masks-01/`, `prospective-chi-guard-preparation-01/`,
`prospective-chi-guard-query-plan-01.json` (failed dependency entry), explicit
`prospective-chi-guard-query-plan-02.json` and
`prospective-chi-guard-query-run-02/`, mean/score/C/guard
plans/runs, `prospective-chi-guard-audit-manifest-01.json`, original selection,
report and pair-summary plan. New raw sources are retained privately in-project,
not registered as tracked production fixtures by this experiment.

`audit_prospective_guard_cohort.py` replays immutable complete-pool choices before
parsing references. It independently recomputes declared deterministic transport
and common geometry from retained data, audits protected ownership/guard phases,
but invokes no new optimizer, external score evaluation or library query. The
separate pair summarizer reads only committed accuracy/selection reports.

```bash
PYTHONPATH=src:. python scripts/audit_prospective_guard_cohort.py select \
  --manifest tmp/pdbfixer-research/prospective-chi-guard-audit-manifest-01.json \
  --output <fresh-selection.json>
PYTHONPATH=src:. python scripts/audit_prospective_guard_cohort.py report \
  --manifest tmp/pdbfixer-research/prospective-chi-guard-audit-manifest-01.json \
  --ledger <fresh-selection.json> --output <fresh-report.json>
```

**795 focused tests passed in 27.04 s**; not the full repository/platform suite.
No production/default/dependency/installed-library changes or public completion
capability, no commit/push. Original rotamer tables, documented appropriate
provenance/license rights and explicit user approval remain mandatory before
production use of that research route, separately from scientific qualification.

The subsequent [full production legacy comparator](legacy-comparator.md) executes
all six shared-admitted masks without restoring observations. Four final heavy-
geometry passes and two THR-CB stereo failures remain visible alongside actual
observed motion; it is separate from the bare seeded stock arm above.
