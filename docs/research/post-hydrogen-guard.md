# Post-H guard (variant 2): source-bound qualification

## Scope and static result

User selected variant 2 after `70cbda6`. The [contract](../plans/post-hydrogen-guard.md)
was frozen at `d7fe7a3`; the complete static implementation/enrollment at
**`8930708`**, before proposal execution. Work remains private on
`research/post-hydrogen-guard`; variant 1 remains intact and available as a
separate alternative. No dependency, production behavior, default or main change.

The [static plan](post-hydrogen-guard-static-plan.json),
[complete outcomes](post-hydrogen-guard-static-evidence.json) and
[independent audit](post-hydrogen-guard-static-audit.json) retain all **18**
post-H returns: six mode-matched unchanged THR2 controls (three sources × two
modes), plus all twelve fixed real-input cases. Only new H may move; the new-H
mask is bound to the original pre-H observations and H-only System mobility,
not inferred from the relaxed geometry. All observations and topology identity
are checked against the original initializer before any proposal.

All 18 states admitted without refusal or execution error. Three states activate:

| Return | Moved atom | Current radius | Result |
|---|---|---:|---|
| original THR2 unseeded, OpenCL | A/THR2/HB | 1.036867 Å | heavy/H sign-compatible proposal |
| original THR2 unseeded, Reference | A/THR2/HB | 1.036089 Å | heavy/H sign-compatible proposal |
| 1VCU protein–ligands, requested AMBER | A/SER0/HA | 1.198827 Å | heavy/H sign-compatible proposal |

The other **15 states are exact no-ops**, including all four seeded THR2
returns. Each changed state preserves all other coordinates exactly; only the
original new H moves and its returned parent–H radius is preserved. The unchanged
candidate rule is the same private selective geometry rule as variant 1, but
applied after H-only minimization. No additional H relaxation is used.

Signs in the independent audit use source bond order and may differ from other
notes' explicit ligand ordering. Heavy versus H-centred volume compatibility is
permutation-invariant; neither sign alone is a chemical CIP label.

### The independent 1VCU activation is not end-to-end qualification

For 1VCU, the AMBER19 H invocation failed template matching and the successful
H-only return came from **forcefield-free fallback**, not native AMBER protein/
ligand physics. Its new A/SER0/HA is sign-incompatible after that fallback;
the static proposal resolves the spatial sign conflict. The corresponding
requested CHARMM case is a no-op. Keep actual routes distinct from requested FFs.

The HA radius is preserved at ~1.199 Å: this spatial guard is **not a C–H bond-
length regularizer** and does not certify that radius. It also does not correct
heavy stereochemistry, protonation, ligand chemistry or all steric contacts.
No matching whole-system Context was captured for this real-input cohort case;
it is therefore **not executed in the exact native whole-phase replay below**.
Do not claim end-to-end 1VCU success or silently construct a different FF/system
and call it the original continuation. A separately defined scientific
continuation is still required for that independent activated context.

## Native follow-up enrollment

The separately [frozen native plan](post-hydrogen-guard-native-plan.json) includes
the two changed original THR2 proposals, **two paths / four new native calls**.
Start from each mode's own saved H-only return, apply the post-H proposal and
the exact historical repair bridge, then run only the original two whole-system
phases. No new H-only call, extra relaxation, FF retuning or restraint reanchoring.
The four exact-noop seeded paths reuse their previous unchanged controls by hash.
The independent 1VCU activation remains explicitly unsupported for this exact
whole-system replay, not dropped from the results.

Six new source/return tests and two multi-path accounting regressions pass
(30 focused tests with the 22 existing guard controls). The new runner writes
each per-path record once and its final summary once, and records worker exit
codes durably; it does not reuse the earlier failed exclusive-write driver.
Private coordinates/States/traces remain under the ignored project tmp tree.
