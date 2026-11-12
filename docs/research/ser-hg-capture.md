# SER-HG: unavailable original artifacts and bounded fresh evidence

## Scope, source and retention

The user authorized access to the failed pytest directory. The read tool still
rejected external-directory access; the user's own copy command then reported
that the source did not exist. Thus the original failed `repaired.pdb` and `.dat`
are **unavailable**, not inspected, recovered or replaced. The retained original
[failure record](thr2-hswap-validation.json), project-local full-suite log and XML
remain authoritative: **3,082 passed / one failed**, OG–HG **1.245813 Å** against
the unchanged **<1.2 Å** assertion. New successful runs cannot turn that run green.

After that loss the user instructed autonomous continuation. This stage uses only
the tracked, synthetic two-SER input literal in `tests/test_prepare_broken_geom.py`.
It is not another protein cohort or a reconstruction of the lost trajectory.
Production code, scientific thresholds, test assertions and safeguards stay
unchanged. Documentation scope: this protocol/results note, a continuation link
from [H-swap](thr2-hswap.md) and the canonical
[handoff](../plans/pdbfixer-reconstruction-handoff.md); validate links and evidence.

## Prospective protocol

The [frozen plan](ser-hg-capture-plan.json) is a byte-identical public metadata
copy of `tmp/pdbfixer-research/ser-hg-capture-plan-01.json`; raw payloads stay private.

Freeze source/runtime/input hashes before scientific execution. Enroll exactly
**six fresh preparation invocations**: one unseeded and H Python seeds **1–5**,
each once. The same `prepare main([input, '-o', output])` call as the tracked test
uses otherwise unchanged defaults and automatic platform routing. Seed and
restore Python RNG only around native `Modeller.addHydrogens` for seeded labels;
the unseeded label preserves natural RNG evolution.

Observe native H minimizer boundaries without adding a reporter, energy query,
force query or changed optimizer arguments. Retain System/Integrator/State XML,
checkpoint, actual platform, topology/bonds and H before/after coordinates. Wrap
the actual production H-repair call to retain its input/output, pure detector
findings and return count, forwarding arguments/results unchanged. Preserve full
RNG state privately. ForceField wrappers retain the native signature/defaults.

All artifacts are written from the start into exclusive directories below
`tmp/pdbfixer-research/ser-hg-capture-run-01/`, not disposable pytest directories.
Retain every outcome, error, timeout and budget nonexecution. Bounds: **eight H
calls and eight native H calls per case**, **600 s per case / 3,600 s total**.
No extra seed, repeat-until-failure, iteration change, repair or acceptance filter.
Every final distance assertion is reported separately from execution success.

## Deterministic detector-contract control

Before examining fresh outcomes, freeze a static local-graph sweep at O–H
distances **0.97, 1.199, 1.2, 1.2458133889150464, 1.551, 1.552, 1.553 and
1.7 Å**. These toy graphs have one H bonded to O, one heavy neighbor and no close
nonparents. Report the existing detector's findings, actual repair count and
post-repair distance, without changing thresholds or running a minimizer.

Current production's distance guard is **0.97×1.6 = 1.552 Å**, whereas the regression
requires **<1.2 Å**. The synthetic 1.245813 Å O–H graph reproducibly fails that
regression predicate yet produces no finding or repair. This establishes a
**contract gap**, not the complete geometry or cause of the lost failure. Other
detector predicates can still repair an atom below that length if it coincides
with a nonparent; do not infer those predicates were absent in the lost output.

## Interpretation boundary

Fresh native captures can show whether an independently generated endpoint lies
in the gap and identify its observed H/repair lineage. They cannot prove the lost
run's stochastic path or qualify a general change in O–H tolerances. No native
sampling is repeated merely to obtain a pass or failure. Production integration,
upstream issue publication and threshold changes are outside this stage.
