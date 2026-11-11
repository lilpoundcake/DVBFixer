# Original-input THR2: bounded fresh native capture

## Prospective protocol and authorization

The user approved this separate bounded original-input/platform capture and
autonomous continuation on 2026-11-12. This supersedes the earlier authorization
boundary in [whole-path research](carbon-H-whole-path.md), not its findings.
Documentation scope: this note (protocol, results, limits) and the canonical
[handoff](../plans/pdbfixer-reconstruction-handoff.md); validate both links and
evidence before publishing. No other documentation expansion is planned.

Before scientific execution, `research_thr2_native_capture.py plan` binds the
original comparator plan, all 160 historical source hashes, original THR2
prepared PDB/Dat/log/output, OpenMM Python source and observer dependencies.
Execution resolves scientific imports from the preserved historical tree.

Exactly **three fresh native CLI executions**: unseeded, H seed 1 and H seed 2,
once each, automatic platform, `--no-solvent --atom-naming standard`, otherwise
historical defaults. Seed only Python random during `Modeller.addHydrogens`;
restore it for seeded arms, preserve its natural evolution in the unseeded arm.
Retain RNG states and before/after H coordinates. No initializer, force,
restraint, minimization tolerance/iteration or safeguard changes.

Observe every actual `LocalEnergyMinimizer.minimize`, including internal H
placement: preserve topology/bonds, System/Integrator/State XML, checkpoint,
actual platform properties, arguments and before/returned coordinates. No native
reporter or energy/force readback is added. Boundary serialization is itself
instrumentation; this is not proof of an identical uninstrumented trajectory.

After **all three native executions**, bind captured artifacts, then replay each
captured call independently **once on captured-checkpoint platform/properties and
once on Reference from XML State** with a read-only native callback reporter.
Never chain replay endpoints into the next phase. Native call parameters remain
unchanged. Record starts, callback arrays/objectives/gradients and returned state;
report OpenCL initial/endpoint differences, never demand cross-platform identity.

Bounds: at most eight native calls per label, 4,096 callbacks per replay,
1,800 seconds per job and 14,400 seconds total. All enrolled labels are attempted;
partial capture, failures, timeouts, missing calls and budget nonexecutions remain
evidence. No retry, adaptive extension, extra seeds or success filtering.

## Interpretation gates

The original failing THR2 run has no saved pre-minimization System/State/checkpoint.
These new captures cannot recover or replace that stochastic trajectory. Failure
to reproduce inversion is a valid result. Reference behavior cannot establish
OpenCL causality; callback observations are not every line-search evaluation.
Native error enums are not exposed by this observer. Production and upstream
issue publication remain outside this authorization.

Private payloads belong under `tmp/pdbfixer-research/thr2-native-capture-run-01/`;
only scripts, tests, bounded aggregate evidence and notes may be tracked.

The [frozen plan](thr2-native-capture-plan.json) is a byte-identical public copy
of `tmp/pdbfixer-research/thr2-native-capture-plan-01.json`; it contains hashes and
protocol metadata, not native coordinates, forces or checkpoint bytes.
