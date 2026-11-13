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

Final fast-suite validation: **3,245 passed**, three optional Torch/Biotite skips
and 33 slow deselections. All 30 focused controls, Ruff (production plus new
research files), strict-target mypy, generated CLI/GUI specs and agent maps pass.
The [validation/quality record](post-hydrogen-guard-validation.json) binds the
complete retained evidence and distinguishes favorable local dynamics from the
independent context's unexecuted continuation and disabled production capability.

## Native result: no extra H relaxation needed in retained THR2

Native enrollment was frozen at **`47a5fcf`** before execution. Both new paths
completed once, with **four returned whole-system phases** and both worker exit
codes recorded as zero. The corrected per-path/final-index accounting succeeded;
no offline summary recovery or native retry was needed. The older variant 1
orchestration failure remains preserved and is not relabelled by this result.
See the [complete independent native audit](post-hydrogen-guard-native-audit.json).

| Own saved H-only return | Before post-H guard | After whole phase 1 | After whole phase 2 |
|---|---|---|---|
| unseeded, captured OpenCL | heavy Cβ R, incompatible HB | Cβ **R** | Cβ **R** |
| unseeded, Reference | heavy Cβ R, incompatible HB | Cβ **R** | Cβ **R** |

The reused unchanged controls invert Cβ R→S during whole phase 1 and remain S.
Here the post-H intervention changes only the returned HB point and retains its
own returned radius. The other 314 returned H and all 327 heavy coordinates
remain exact at the proposal boundary; later whole-system motion follows the
original native masks. No H-only minimization is rerun, no extra relaxation is
added, and no early OpenMM hook is used.

Independent same-platform physical readback matches all **278 callbacks** with
zero maximum mobile-gradient/objective error. All four phase Systems are
byte-identical to the captured originals, including the original same-anchor
divide-ten phase transition. Both historical bridges make zero repairs;
requested/applied coordinate bounds, copied observations, own-return continuity
and serialized return States pass the independent audit. No negative Cα node
or last-callback/returned position mismatch appears. The complete static and
native offline audits each repeat byte-identically, without new trajectories.

All **four whole-system returns meet the native gradient criterion**. The saved
upstream H-only returns still came from the unchanged 50-iteration runs that did
not meet their measured native gradient criterion: this stage does not erase or
improve that stopping limitation. It shows that extra H-only relaxation was not
necessary to retain R under these particular subsequent whole-system settings.
Final target C–H distance is ~1.091 Å in both modes; that is an observed endpoint,
not a newly imposed constraint or a general geometry-acceptance certificate.

This supports variant 2 as a **local working alternative** to the early hook.
The two modes remain one poor-heavy context, not two independent structures or
broad chirality protection. Nonparent overlaps, protonation and unsupported
chemistries are not comprehensively qualified by this finite comparison.

## What remains before integration

For THR2 there is currently no experimental reason to require a patched early
OpenMM hook: post-H correction plus the existing whole stages retained R without
extra minimization. Keep variant 1 available as a fallback research route rather
than deleting or weakening its evidence.

The independent **1VCU A/SER0/HA** activation is still only statically assessed.
A future whole-complex experiment must first declare its System construction,
actual ligand/FF/fallback route, comparison controls, source-bound handoff through
late heterogen-H placement/variant naming, budgets and failure retention. It may
not be counted as an exact historical replay or as an end-to-end success here.

Supported production integration would also need a before/after public-API
ownership adapter: snapshot observations before `addHydrogens`, resolve exact
identities afterward, and allow changes only to genuinely new H. The private
mass-based source admission relies on saved internal Systems and is **not itself
a public runtime integration**. Deletions, variant aliases, repeated chain
occurrences, partial observations and exception cleanup need explicit handling
and live API qualification. No global minimizer monkeypatch or installed OpenMM
modification is part of the proposed post-H path. Public/default capability,
main merge and release remain disabled.

In particular, minimize's required variant strip/readd path may replace H again
after prepare. A live workflow must audit each H-creation boundary and retain
observation lineage, rather than assume a prepare-time correction survives or
classify a stripped observed H as scientifically new. Do not bypass canonical
variant/bond handling or force the keep-H fast path to obtain a favorable result.
