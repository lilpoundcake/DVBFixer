# Native Reference replay: returned Context can retain a rejected line-search trial

Status: **controlled observational replay plus source-level forensics**, not a
production fix or explanation of the original unforced 1CRN THR2 Cβ event.
Follows the [saved-only whole-force attribution](carbon-H-whole-attribution.md).

## All 24 exact computational replays

The [complete independent readback](carbon-H-whole-path-evidence.json) retains
every matched path, callback and returned-versus-reported comparison below.

The same two seeds × three THR24 conditioning contexts × four original H arms
were replayed once, with the unchanged original two **Reference / no-solvent /
tolerance 10 kJ/mol/nm / maximum 1,000 iterations per phase** calls. The original
State/integrator, forces, anchors and tenfold phase-2 positional-weight reduction
remain fixed. Native callbacks only copy their coordinates, gradients and objective
metadata; they never request an early exit. The original post-H whole start is
used directly, without running H repair again or regenerating H.

All **24 paths and 48 phases** complete without blockers, timeouts, retries or
budget nonexecutions. Every original start and both original nm endpoints
reproduce exactly. All **3,991 callbacks** independently match native physical
forces and energies, maximum component gradient error **0.0 kJ/mol/nm**. All
**4,087 callback/boundary nodes** retain target Cα, target Cβ and GLN16 Cα volumes,
force projections, native gradient criterion and terminal comparisons. Constraint
penalty is zero; physical positional restraints remain part of native energy.

The largest callback count in any phase is **188**, below the frozen 1,000
iteration cap. Storage bound is 2,048 callbacks per phase, not an optimization
extension. Budgets are 1,800 seconds/path and 14,400 seconds total. Partial callback
records are retained on exceptions or endpoint mismatch, without a success trace.
These are computational replays, not new independent scientific samples.

Plan SHA256:
`748b13a82a83656dd99677084ac19b1b08e490c69f066adf8c5a0dfa8e0e8bea`.
Private plan/run/audits: `tmp/pdbfixer-research/carbon-H-whole-path-*`.
No raw coordinate/gradient arrays or candidate PDBs are redistributed.

## Returned state versus last reported iterate

In **17/48 phases**, returned Context coordinates differ from the last reported
iterate. All 17 returned states fail the native gradient convergence criterion;
the other 31 phases satisfy it. Five returned states have higher energy than the
last reported iterate; twelve have lower energy but still do not establish a
successful line search or convergence. All phases have at least one reported
iterate, so no unreported initial state is substituted in these comparisons.

All three negative-Cα boundary snapshots have a positive last reported target Cα:

| Context | Phase | Callbacks | Last reported Cα | Returned Cα | Returned minus reported energy |
|---|---|---:|---:|---:|---:|
| Seed 1, beta-only, methyl-zero | 1 | 26 | +0.021875 Å³ | **−0.002076 Å³** | **+357.415 kJ/mol** |
| Seed 2, raw, CH₂ | 2 | 7 | +0.020362 Å³ | **−0.006118 Å³** | **+943.606 kJ/mol** |
| Seed 2, raw, combined | 2 | 7 | +0.017283 Å³ | **−0.006974 Å³** | **+1093.030 kJ/mol** |

For the two seed-2/raw endpoints, maximum per-atom displacement from the last
report is only **0.005301 / 0.004996 Å**. Near an already collapsed heavy frame,
that small displacement can change the volume sign and activate a very stiff
relative guard. Last reported gradient ratios are about **40.017 / 39.864**,
returned ratios **122230.397 / 131588.604**, versus native epsilon about **3.270**.
Neither the returned nor last reported state should be called a converged
geometrically valid repair. We do not restore or publish either as a candidate.

## Pinned source explains a rejected-trial Context leak

[Source provenance](openmm-minimizer-source-36a30cb.json) records the revision
reported by this runtime: **8.5.2.dev-36a30cb**, commit
`36a30cbca54e727b216b606f3c011b67201eb8b4`, Darwin arm64. Exact upstream files and
copyright notices are retained privately; only source URLs/checksums are public.
Runtime revision metadata is **not** independent attestation of an unmodified
compiled binary. No installed runtime or production file is patched.

At that declared revision:

1. Reference `evaluate()` writes each trial position into the Context before
   computing forces/energy.
2. [L-BFGS on a failed line search](https://github.com/openmm/openmm/blob/36a30cbca54e727b216b606f3c011b67201eb8b4/libraries/lbfgs/src/lbfgs.cpp#L483-L489)
   restores its own `x` and `g` buffers to the previous point and returns an error.
   The progress reporter is called only after a successful search.
3. [The Reference wrapper](https://github.com/openmm/openmm/blob/36a30cbca54e727b216b606f3c011b67201eb8b4/platforms/reference/src/SimTKReference/ReferenceMinimize.cpp#L196-L219)
   ignores the L-BFGS return value and reads positions back from the Context,
   without restoring the accepted `x` into that Context.
4. For these unconstrained/all-mobile Systems, the constraint-error check passes
   immediately. No subsequent constraint projection restores the accepted point.

This source-level mechanism explains how the returned Context can expose the
last rejected trial rather than the optimizer's rolled-back buffer, consistent
with the exact replays and higher-energy negative endpoints. The exact numerical
L-BFGS error enum is **not observed** by the Python API here; do not assert a
specific maximum-line-search/minimum-step code or that the 1,000-iteration cap
was reached. No generalization to other platforms, builds or OpenMM versions is
established by this Reference evidence.

## Accepted iterates also have adverse transient geometry

The terminal leak is not the only stereochemical issue. Among reported iterates
there are **16 negative target-Cα measurements across eight phase-1 paths**, and
**137 negative target-Cβ measurements across 16 phase-1 paths**. All Cβ paths
are the raw/beta-only stress strata (eight each); none is full-generated.
The most negative target-Cβ volume is **−19.101650 Å³**, a severely distorted
intermediate, not evidence of a valid D-THR minimum. Counts are correlated callback
measurements, not independent D residues or continuous-time durations.

All saved whole boundary THR-Cβ signs remain positive; the original unforced
THR2 Cβ endpoint is not reproduced or explained by transferring this result.
Do not claim that every adverse sign is introduced only after return, that simply
restoring the last report guarantees valid geometry, or that a finite Cα guard
protects a Cβ centre. These intentionally poor-heavy conditions remain distinct
from newly admitted healthy production inputs.

## Validation and authorization boundary

Eleven focused reporter/readback/source-provenance tests pass, including exact
native endpoints, interrupted partial-record retention, endpoint mismatch
refusal, no-minimization audit, nonfinite input refusal and unversioned-source
rejection. Targeted Ruff passes. The full fast suite through this checkpoint
passes **3,066 tests**, three optional Torch/Biotite skips and 33 slow deselected.
Retained log: `tmp/pdbfixer-research/carbon-H-whole-path-fast-suite-01.log`.

The current fixed Reference cohort is exhausted for this replay question. To
continue toward the **original unforced 1CRN THR2 cause**, obtain its original
failing pre-minimization System/State/checkpoint, or request explicit approval
for a separate prospectively bounded fresh capture on the original input/platform.
The original comparator folder has input/prepared/output PDBs and logs but no
such original checkpoint. Do not silently expand seed/platform search until an
inversion appears, transfer the Reference cause to that uncaptured event or
rewrite any frozen request/evidence.

Production integration or publishing an upstream issue is a separate authorization
decision. No production/backend/default/threshold/safeguard has changed; all
empirical-table production provenance/license/approval gates remain closed.

Code: `scripts/research_carbon_H_whole_path.py`,
`scripts/audit_carbon_H_whole_path.py`, `scripts/capture_openmm_minimizer_source.py`.
