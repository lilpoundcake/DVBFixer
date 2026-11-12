# GLYCAM fallback graph restoration: bounded qualification

## Outcome and authorization

The user approved this separate fix after the original
[O–H real-input cohort](oh-native-qualification.md) retained an 11/12 preparation
result. Branch: `fix/glycam-charmm-hydrogen-fallback`. The prospective contract
was committed at `f0911a3`; implementation and the
[new frozen plan](../plans/glycam-hydrogen-fallback-qualification-plan.json) were
committed at `d735083` **before** native execution. No main merge or release.

All **12 new declared preparations completed once**, with no replacement cases,
new seeds or native retries. All twelve postflight states replay exactly and
all final raw PDB atom identities are unique; every `.dat` loads through
`DatRecord`. This closes the declared **changed-candidate cohort gate**, not
broad scientific or native CHARMM glycan qualification.

- [Complete qualification evidence](glycam-fallback-qualification-evidence.json)
- [Exact saved-failure graph control](glycam-fallback-saved-state.json)
- [Current PDBFixer source-binding amendment](pdbfixer-usage-inventory-glycam-fallback.json)
- [Final validation and retained-artifact digests](glycam-fallback-validation.json)
- [Implementation contract and fixed budgets](../plans/glycam-hydrogen-fallback.md)

## Mechanism and what the implementation changes

In the original failed CHARMM-requested GLYCAM hinge/CH3 preparation, the saved
fallback input has one existing but **unbonded H1** at each of `H/UYB232` and
`K/UYB242`. C1 is also missing its intra-residue O5 and C2 links. OpenMM matches
existing H by bonded parent; the forcefield-free pass therefore creates a second
H1. The old postflight performs zero repairs, and final naming refuses duplicate
identities. The complete original failure and its false cohort status remain
unchanged and retained.

`prepare/fallback.py::restore_glycam_fallback_bonds` runs only inside the existing
whole-FF-error fallback. It selects canonical GLYCAM residues with unbound
existing H, builds isolated metadata views through `build_glycam_system`, and
validates exact-name templates before committing any live bond additions.
Duplicate names, unknown templates/atoms, incompatible elements or H bonds,
missing parents and incomplete wrapper restoration produce a typed error.
They do not trigger another H attempt.

Only missing **intra-residue** bonds are added. No atom is created, deleted or
renamed; coordinates and inherited bonds, including all external links, remain
unchanged. The metadata view needs the normal AMBER14+GLYCAM XML dependencies
because GLYCAM protein templates reference AMBER atom types. It never creates
a System, samples H or supplies the fallback's minimization parameters.
Fallback remains explicitly **forcefield-free**; successful whole-FF paths and
unknown non-GLYCAM ligands retain the previous route.

The early engineering attempt with GLYCAM XML alone failed on `protein-H1`.
That source/error was retained privately before loading the standard dependency
set in the isolated metadata view. No native case was run on that engineering
attempt, and no force-field parameter was edited.

## Exact old-state control (no sampling)

The preserved original pre-fallback state receives exactly six additions:

| Residue | Added links |
|---|---|
| H/UYB232 | C1–H1, C1–O5, C1–C2 |
| K/UYB242 | C1–H1, C1–O5, C1–C2 |

Every atom, exact nm coordinate and inherited link is unchanged. There are no
new external links; the second recovery call adds zero bonds. This deterministic
control establishes the graph correction independently of native variability.
It does not recover any lost historical THR2/SER trajectory.

## All twelve changed-candidate outcomes

The same six sources/two requests, pH, H seed label 1, water retention and standard
naming were enrolled prospectively. All actual H minimizer boundaries used the
unchanged automatic OpenCL / Apple M3 Pro / single-precision route. Limits stayed
1,200 s per case, 7,200 s aggregate and eight H/native-H calls per case.

Repair counts below are old/new **O–H postflight-module** replay counts and include
nonoxygen legacy repairs. They are not baseline/candidate graph-repair runs or
paired native trajectories. The saved-old-state control above qualifies graph
restoration; the changed-candidate cohort qualifies bounded integration.

| Case | Source / requested FF | Successful H route | Graph links added | Old / current postflight repairs | Preparation |
|---|---|---|---:|---:|---|
| 0 | 1UBQ / amber | AMBER19 | not entered | 4 / 4 | completed |
| 1 | 1UBQ / charmm | CHARMM36 | not entered | 7 / 7 | completed |
| 2 | 1DQJ / amber | AMBER19 | not entered | 2 / 2 | completed |
| 3 | 1DQJ / charmm | CHARMM36 | not entered | 2 / 2 | completed |
| 4 | 8CDE pH 3 / amber | forcefield-free after terminal-GLU template failure | 0 | 0 / 5 | completed |
| 5 | 8CDE pH 3 / charmm | CHARMM36 | not entered | 11 / 11 | completed |
| 6 | 3RY6 / amber+glycam | AMBER14/GLYCAM | not entered | 31 / 31 | completed |
| 7 | 3RY6 / charmm | AMBER14/GLYCAM bridge, then CHARMM naming | not entered | 32 / 32 | completed |
| 8 | 1VCU / amber | forcefield-free after DAN template failure | 0 | 4 / 4 | completed |
| 9 | 1VCU / charmm | forcefield-free after DAN template failure | 0 | 2 / 2 | completed |
| 10 | GLYCAM hinge/CH3 / amber+glycam | AMBER14/GLYCAM | not entered | 0 / 0 | completed |
| 11 | GLYCAM hinge/CH3 / charmm | forcefield-free after NLN template failure | **6** | 0 / 0 | **completed** |

Only the intended recovery in case 11 adds live bonds, exactly matching the
six-link old-state control. All four recovery calls preserve atom identity,
exact nm coordinates, inherited typed/ordered bonds and external links. The
remaining three return zero additions.

The previously failing case now emits **4,117 unique atom records** and a valid
sidecar, without duplicate sugar H1. It does not pass native CHARMM template
matching: that expected upstream failure/fallback remains explicitly recorded.
The metadata load is not relabelled as a native AMBER H-placement success.

All twelve candidate postflight states satisfy the coarse O–H contract and
preserve heavy/valid-H coordinates and graph identity. There are no new typed
refusals, short-O–H warnings (<0.8 Å) or ≤1.0 Å cross-residue-heavy warnings in
these states. Five additional O–H violations are repaired in the pH-3 fallback;
this is a new saved sample, **not a rate comparison** with the original cohort's
twelve violations. Native OpenCL endpoints and repair counts can vary between
otherwise matched native preparations. Terminal-GLH normalization warnings and
other ordinary pipeline diagnostics remain visible.

## Observer and audit engineering

The new observer retains explicit native chain/residue occurrence indices,
exact nm floats and bond types/orders. It preserves repeated case-sensitive
chain IDs without the old auditor's inferred-run transport limitation. Original
observers and artifacts remain frozen rather than overwritten.

Seventeen focused helper/observer tests cover exact expected links, idempotence,
unrelated ligands/already-bonded no-ops, duplicate names, wrong/missing parents,
element mismatch, all-or-nothing validation, template/wrapper failures, no native
calls during recovery, failure propagation, chain occurrences/case/insertion
codes and typed-bond transport. The expanded focused set passed **112 tests**.
The initial red helper test failed on the then-missing module; the duplicate-H
mechanism itself is established by the retained real snapshots, not that import
failure alone.

An independent private toy audit detected removal of an inherited external bond
while accepting the exact recovery, and never called an incomplete toy cohort
qualified. Repeating the offline audit consumes the same saved states, not a
new native sample; its complete report is byte-identical.

The final complete fast suite passed **3,202 tests**, with three optional
Torch/Biotite skips and 33 slow deselections. Ruff, strict-target mypy, generated
CLI/GUI specs and agent-map checks passed. Original SER assertions, native FF
parameters and minimizer settings were not changed. The separate current-source
inventory updates only the prepare source binding and preserves all historical
inventories, PDBFixer import/use semantics and totals.

## Limits and retained evidence

The original cohort remains **11/12 and not qualified**; this separate changed
candidate is **12/12 under its declared gates**. Neither result is rewritten.
This is six input structures and one H draw label per request, not broad
accuracy, protonation, convergence or hardware-repeat qualification. Direct
late ligand-H coverage and native CHARMM glycan qualification remain absent;
the separate THR2 HB problem and empirical-library licensing gates are unchanged.

Private sources, failed engineering attempt, raw/native/RNG states, recovery
snapshots, original failure hashes and complete new outputs live below
`tmp/hydrogen-postflight-contract-01/glycam-fallback-01/`. Original failed case
artifacts remain at their original project-local paths. Tracked evidence carries
only identities, statuses, hashes and scalar measurements, not coordinates or
native arrays/checkpoints. Implementation stays on the fix branch, not main.
