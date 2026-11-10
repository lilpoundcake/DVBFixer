# Generated-only target-χ relaxation: protected geometry holds, choices unchanged

Status: **bounded retrospective research**, not a production refiner or accuracy
qualification. Continues the [saved-context/pair diagnosis](mean-context-pairs.md)
and [library-mean proposals](dunbrack-means.md). All original scientific failures
and the research-only PyRosetta/original-table production gate remain intact.

## Frozen intervention and ownership

`scripts/research_generated_chi_relax.py plan` froze **all 42 geometry-feasible
mean starts in 15 pools**, retaining the other seven original masks as inherited
nonexecutions. No native-error or score-margin filtering. The same old eleven
and four executable new-input masks are already evaluated; this is not a new
prospective accuracy cohort.

Each initial fixed-heavy pose must reproduce its original numerical score and
complete compressed snapshot exactly (JSON container types normalized, no numeric
rounding). Then enable only the target residue's χ degrees of freedom in a
MoveMap. All backbone/jump/other-residue χ degrees remain disabled.

Only the generated **distal** VAL/THR heavy atoms and target sidechain H may move.
Generated **CB stays protected**, as do all observed atoms, backbone H, non-target
context H, scorer extras and virtual atoms. No packing, Cartesian minimization,
reflection, source/context repair or coordinate restoration. Active installed
library hash must match the previous queries. Weights remain unchanged.

Predeclared minimizer: `lbfgs_armijo_nonmonotone`, tolerance 0.01, maximum 200
iterations, torsion-space mode, neighbor list disabled. The run does **not** claim
full convergence or global minimization. The protected-coordinate numerical bound
is 1e−9 Å; actual protected motion is **0.0 Å in every start**.

## Important score semantics

The final score belongs to the **raw relaxed pose including its target sidechain
H**, not a fresh import of rounded heavy atoms with newly rebuilt hydrogens.
This is an explicit change from the earlier fixed-heavy score-only protocol, not
a hidden replacement for it. The final heavy proposal must pass independent
raw **and** rounded common geometry before ranking. Source PDB records are kept
unchanged in the diagnostic materialization.

The independent auditor blocks selection if any expected start in a pool has an
execution/mobility blocker. It ranks all accepted refined candidates by final raw
score, stable original index tie. All pool selections are written to the
[selection ledger](generated-chi-relax-evidence-selection.json) before native
reporting; [accuracy evidence](generated-chi-relax-evidence.json) binds that ledger.
This is a coupled distal-heavy/sidechain-H intervention, not an isolated test of
whether hydrogen optimization alone explains the previous ranking errors.

## Results

| Measure | Result |
|---|---:|
| Original masks / executable pools / inherited nonexecutions | 22 / 15 / 7 |
| Exact initial-state replays / completed bounded refinements | 42 / 42 |
| Execution or protected-movement blockers | 0 |
| Protected atom maximum motion | **0.0 Å** |
| Final common geometry passes / rejections | **41 / 1** |
| Changed selected candidate indices | **0/15** |
| Selected rounded heavy PDBs unchanged from their starts | **15/15** |

Every raw refined state retains CA S, and every THR retains CB R; no global
stereo findings or chirality reflection. One nonselected 2CI2 THR55 start (index 0)
reduces its score but introduces an OG1–I/52 O clash and is rejected by the
unchanged common gate. This illustrates why lower model energy cannot replace
the independent output audit.

Substantial generated-heavy motion occurs in only three nonselected starts:
1CTF THR76 index 1 (0.276725 Å), 1R69 VAL6 index 2 (0.152352 Å), and 2CI2 THR55
index 0 (0.778725 Å). The other starts move little under this declared bounded
configuration. Do not label this a no-op, a converged minimization, or a general
failure of every possible refinement method.

Selected rounded generated-atom RMSD medians, Å:

| Cohort | Executable / planned | Exact names | Generated-only VAL symmetry | Parent exact |
|---|---:|---:|---:|---:|
| Old eleven | 11 / 11 | 0.198491 | 0.198491 | 0.198491 |
| New-input eleven | 4 / 11 | 1.930166 | 1.463808 | 1.930166 |

All selected indices and rounded heavy outputs are unchanged, so this bounded
intervention provides **no selected-output accuracy rescue** on the retained
cohort. Native-control caveats and seven unsupported parent requests remain;
no substitution, threshold change or omission improves the denominator.

## Independent readback and checks

`scripts/audit_generated_chi_relax.py` does not invoke PyRosetta, minimize, rebuild
or evaluate scores. It binds each complete initial pose to its saved all-atom
mobility trace, independently re-derives allowed/protected ownership, validates
before/after identities and motion, checks unchanged score weights/term totals,
re-audits raw/rounded geometry and chirality, and replays strict complete-pool
selection before reading references. Raw failures and atom traces remain private.
No individual library means, weights, sigmas or tables are published.

Retained paths under `tmp/pdbfixer-research/`:
`generated-chi-relax-plan-01.json`, `generated-chi-relax-run-01/`, intermediate
audits and the final independent readback. The historical fixed-heavy experiments
and exact-state pair diagnosis remain unchanged and separately reproducible.

```bash
PYTHONPATH=src:. python scripts/audit_generated_chi_relax.py \
  --plan tmp/pdbfixer-research/generated-chi-relax-plan-01.json \
  --run tmp/pdbfixer-research/generated-chi-relax-run-01 \
  --output <existing-new-directory>/generated-chi-relax-evidence.json
```

Use fresh destination files. **738 focused tests passed in 27.68 s**; not the full
repository/platform suite. Production code/defaults, installed libraries,
dependencies, source fixtures and public capability remain unchanged; no commit
or push. Any further convergence/mechanism control must be separately frozen,
not tuned to these known native errors. Original tables, documented appropriate
license rights and explicit user approval are still required before production.

Subsequent [native-free stopping-policy controls](chi-stopping-controls.md) test
three predeclared numerical settings on every original feasible start. They
demonstrate sensitivity of coordinates/choices without evaluating native accuracy
or selecting a preferred profile. The bounded baseline result above is unchanged.
