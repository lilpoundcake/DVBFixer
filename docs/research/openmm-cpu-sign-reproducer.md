# Isolated synthetic OpenMM constraint-gradient research

## Scope and first prospective series

User instructed private continuation after `330b40a` and the explicit isolated-
reproducer proposal. Work is on `research/openmm-cpu-sign-reproducer`. No installed
OpenMM patch, new environment/package installation, production/main change,
external issue or public capability. Independent child processes use exclusive
project-owned output/temp/cache directories. Native-library directory inspection
was permission-denied and not bypassed; runtime metadata is not binary-hash
attestation. A preliminary unsupported Platform.getLoadedPlugins introspection
was an engineering failure, not a scientific case.

The [first contract](../plans/openmm-cpu-sign-reproducer.md) and
[12-case plan](openmm-cpu-sign-reproducer-plan.json) were frozen at **`287eacb`**
before native execution. Three positive-mass particles, one harmonic bond,
independent quadratic stress, optional pair constraint; no biological structures
or force-field templates. OpenCL-single/CPU/Reference × constraints on/off ×
low/high stress all ran once, with public nonterminating reporters. Analytic
energy/gradients and termwise finite differences are independent of OpenMM.

## Complete first result: desired trigger not reproduced

[All 12 outcomes](openmm-cpu-sign-reproducer-evidence.json) have worker exit0 and
returned native calls, **162 callbacks** total. Finite differences pass all162.
The predeclared primary source-sign signal is **false**; all-controls gate is
also **false**. These are retained negative results, not discarded tuning runs.

- All eight CPU/Reference controls pass physical gradient/energy envelopes.
- OpenCL free-low passes. OpenCL constrained-low passes its gradient envelope
  with the predeclared constraint-arithmetic allowance, but one of43 callbacks
  fails physical energy: error~0.00148483 kJ/mol. Keep this failure, not a wider
  absolute bound. The allowance is not strict gradient equality.
- Both OpenCL high-stress cases produce one callback, with no reported gradient
  above the fixed-point capacity criterion. Their stress-particle derivatives
  are wrong even WITHOUT constraints: an analytic x gradient~6.70954e13 is
  reported as~-6.16563e8. None of the three declared gradient hypotheses matches.
  This is consistent with force-buffer aliasing at excessive forces, not an
  observed CPU-sign reproduction. No internal buffer/fallback event is directly
  witnessed.
- Positive penalty energy/gradient finite differences remain correct. Neither
  source-sign matching in old 1VCU readbacks nor new synthetic native completion
  rehabilitates physical qualification or HA286 geometry.

The huge-stress trigger is therefore not a minimal sign reproducer on this
runtime. Do not rerun these unchanged cases until green, alter their recorded
classification, or modify the old code/tolerances after results. A separately
declared representable-force threshold control is needed to distinguish the
source's overflow detector from out-of-range force aliasing. No dependency patch
or corrected native gradient has been tested.

## Separate representable-force study: observed public-reporter inconsistency

The [threshold-control amendment](../plans/openmm-representable-force-control.md)
and [new 12-case plan](openmm-representable-force-plan.json) were frozen at
**`63adfc2`** before native execution. Pinned `convertForces` source requests
fallback above2^30 force units, below the signed fixed-point capacity2^31.
The new particle2 start and stress values put the trigger inside that interval;
maxIterations1 is explicitly diagnostic, not a convergence claim. These are
new cases/conditions, not repetitions or rewritten outcomes from the first series.

[All second-series outcomes](openmm-representable-force-evidence.json) completed
once: 12 native calls, one callback each. The predeclared primary witness is
**observed**, and **all11 other controls pass** their physical/energy/finite-
difference criteria. No missing paths, platform substitutions, timeouts or retries.
The complete two-series accounting is24 native calls /174 callbacks, with the
first negative series and its failed controls preserved unchanged.

| Second-series group | Result |
|---|---|
| OpenCL constrained representable-trigger | physical gradient FAIL, source-sign PASS |
| OpenCL constrained below detector | physical PASS, both source-sign hypotheses FAIL |
| OpenCL unconstrained below/trigger | physical PASS |
| CPU/Reference constrained below/trigger | physical PASS, source-sign hypotheses FAIL |
| CPU/Reference unconstrained below/trigger | physical PASS |

The constrained OpenCL trigger reports maximum gradient1,201,164,672, strictly
between2^30 and2^31. For its constrained particle0 x component:

- physical derivative: **199.4642138481195** kJ/mol/nm;
- independent termwise finite difference: **199.46421384833772**;
- reported native gradient: **200.52517700195312**;
- source-sign prediction: **200.52517652511037**.

The physical mismatch~1.060963 exceeds BOTH the original strict bound~0.005011
and the predeclared extra OpenCL constraint-arithmetic allowance~0.762941. This
is not rescued by the existing rounding allowance. The complete source-sign
prediction passes strict comparison; after float32 transport it matches **all9
reported components exactly**, maximum residual0. Physical and penalty energy
errors at that node are both0. This wrong-sign formula is NOT the derivative of
the positive penalty energy, as the independent finite differences verify.

This establishes a bounded **public-reporter/objective-gradient inconsistency
on the installed runtime** (self-reported8.5.2.dev-36a30cb, OpenCL single,
Apple M3 Pro). It strongly supports the inspected CPU-overflow source-sign path,
but does not directly observe the internal fallback event or attest installed
binary/source identity. It does not establish the cause of HA286 inversion,
generalize to every release/platform, certify constrained stationarity, fix the
dependency or rehabilitate the earlier structural qualification.

## Portable compact reproducer and scope of its verification

`scripts/openmm_cpu_sign_minimal.py` is a compact, standalone stdlib+OpenMM
wrapper with no dvbfixer imports, private files, biological structures or NumPy
dependency. Copy it to an OpenMM environment and explicitly run:

```bash
python scripts/openmm_cpu_sign_minimal.py --run-native --platform OpenCL
# Platform controls (separate explicit invocations): --platform CPU / Reference
```

Default invocation describes only and creates no Context/minimizer. This wrapper
was packaged AFTER the factorial study and was **not run for an extra native
trial**. Its System serialization is byte-identical to the frozen model; its
scalar analysis was replayed offline against the unchanged retained native
callback and reproduces physical FAIL / exact float32 source-sign agreement.
The actual native executions are those of the frozen full drivers, not a
claimed new end-to-end wrapper execution. Three packaging controls plus10 pure
math/System controls pass. Do not treat a diagnostic native return as convergence.

No dependency replacement, installed-source edit, protected library-directory
access, upstream issue or production change was made. Further action now needs
an explicit choice between upstream reporting and separately scoped isolated
dependency-fix testing; this record does neither automatically.

## Final validation and preservation

Full fast suite: **3,305 passed**, three optional Torch/Biotite skips, 33 slow
deselections. All62 combined focused controls, Ruff, strict-target mypy,
generated CLI/GUI specs and agent maps pass. The
[validation/quality record](openmm-cpu-sign-reproducer-validation.json) retains
both native execution indexes, frozen source/evidence digests and wrapper/raw-
callback replay binding. Local links and public raw-payload absence pass; all
planned contracts/results/validation/handoff exist with no pending links.

The first synthetic signal/control failure remains negative; the second reproduces
the reporter inconsistency with controls. Neither means physical/structural
qualification passed. Old post-H local-policy, strict-guard, variant1/2 and
production source/fixtures/assertions remain unchanged. No jobs are running.
