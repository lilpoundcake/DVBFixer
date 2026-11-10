# Bounded actual-endpoint χ gradient guard: 39/42 meet the numerical criterion

Status: **private native-free numerical-control experiment**, not production,
native accuracy, general convergence or global-minimum qualification. Follows
[native derivative checkpoints](chi-checkpoints.md); the
[research-only licensing/production gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03)
remains unchanged.

## Frozen intervention

`research_chi_endpoint_guard.py plan` froze **all 42 original mode-C starts**,
all 15 pools and seven inherited nonexecutions. No native/error/gradient/geometry
filtering. The inner optimizer remains the same absolute-tolerance L-BFGS mode,
tolerance 1e−6, gmax 1e−4 **per degree**, maximum 200 iterations per call,
torsion-space mode, no neighbor list and unchanged objective/weights.

Run one inner cycle. On an independent pose and ScoreFunction clone, compute the
gradient of the **actually returned endpoint** with the same target-only
MinimizerMap/AtomTreeMultifunc API. If the maximum absolute component is at most
1e−4 per degree, stop; otherwise reenter the identical inner optimizer. Limit:
**three total cycles**, fixed before execution. No new initial angles, candidate
families, tolerance grid or best-profile selection. Three cycles is a diagnostic
resource bound, not a convergence guarantee or an accuracy-selected setting.

Each first cycle must reproduce the original mode-C **complete returned atom
trace and raw score exactly**. Only generated distal heavy atoms and target
sidechain H may move. Generated CB, observations, backbone H, context H/extras/
virtuals stay protected. No reflection, observation restoration or context repair.
Every phase/gradient is retained; independent raw/rounded geometry is checked
after execution and does not select which phases continue.

## Results

[Independent evidence](chi-endpoint-guard-evidence.json):

| Measure | Result |
|---|---:|
| Starts enrolled/audited; execution blockers | 42/42; 0 |
| Exact first-cycle mode-C replays | **42/42** |
| Starts requiring more than one cycle | 4 |
| Total inner calls / extra calls | 49 / 7 |
| Final actual-endpoint gradient criterion met | **39/42** |
| Three-cycle budget exhausted | **3/42** |
| Final common geometry passes | **40/42** |
| Both final gradient and common geometry pass | **37/42** |
| Protected movement / global stereo findings | **0.0 Å / 0** |

All 49 retained phases keep CA S and THR CB R. The common gate passes 47/49
phases. The same two final geometry rejections remain: 2CI2 THR55/start 0
(OG1–I/52 O clash) and 1CTF VAL56/start 2 (CG2–own backbone C clash). Meeting the
gradient criterion does not make either geometry acceptable.

### The four continued starts

Actual maximum gradient components, score units per degree:

| Case/start | Cycle 1 | Cycle 2 | Cycle 3 | Outcome |
|---|---:|---:|---:|---|
| 2CI2 THR22/2 | 0.001844293 | 0.000250206 | 0.000124301 | Budget exhausted |
| 1R69 THR27/1 | 0.000296562 | 0.000194645 | 0.000126400 | Budget exhausted |
| 1CTF VAL56/0 | 0.007424845 | 0.000002814 | Not run | Criterion met |
| 1CTF THR76/2 | 0.000951026 | 0.000697287 | 0.000494669 | Budget exhausted |

All four retain acceptable common geometry in each executed cycle. The guard
successfully distinguishes a small **previous** gradient from the returned
endpoint's gradient, and it corrects the VAL56 numerical stopping result within
the declared budget. But **three endpoints remain above the cutoff**. They are
not relabelled as converged, omitted from the denominator, or assigned more
iterations after seeing the result. THR22 is still above the cutoff even though
its gradient has decreased substantially.

## Interpretation and limits

This is a bounded test of an actual-endpoint stopping rule, not a validation of
the raw force-field minimum, stable candidate ranking, native reconstruction or
every analytic derivative. Raw scores still belong to the H-containing refined
pose, not a fresh rounded-heavy import. No candidate ranking, preferred profile,
new PDB publication or native RMSD evaluation is performed.

The actual-endpoint threshold is an explicitly measured numerical condition;
even its 39 passes are **not** a certificate of global minimum, general optimizer
convergence, or tolerance-independent coordinates. The prior finite-derivative
comparison caveat remains. All evidence is repeated numerical research on
already evaluated inputs, not an independent scientific cohort.

The derivative calculation uses separate clones and checks objective agreement
within the frozen 1e−8 absolute score bound; the mapped clone's numerical motion
bound is 1e−9 Å. The retained scientific endpoint is never reset after derivative
evaluation. The independent auditor verifies the saved phase coordinates,
ownership, continuation decisions, complete gradient inventory, all final score
weights/term totals, and the final result's binding to its last phase. It reads no
native reference, invokes no PyRosetta or minimization, and preserves every
execution failure, budget exhaustion and geometry finding.

## Reproduction

Private plan/run: `tmp/pdbfixer-research/chi-endpoint-guard-plan-01.json`,
`chi-endpoint-guard-run-01/`. Every `cycle-N.json` retains the complete atom trace,
raw score, gradient and decision. No raw library/table parameters are published.

```bash
PYTHONPATH=src:. python scripts/audit_chi_endpoint_guard.py \
  --plan tmp/pdbfixer-research/chi-endpoint-guard-plan-01.json \
  --run tmp/pdbfixer-research/chi-endpoint-guard-run-01 \
  --output <fresh-report.json>
```

**782 focused tests passed in 27.43 s**; not full repository/platform acceptance.
No production consumers/defaults/dependencies, fixtures, installed libraries or
public capability changed; no commit/push. Original tables/provenance/suitable
license rights and explicit user approval remain mandatory before production,
independently of the still-open scientific qualification gates.

Subsequent [pre-enrolled new-input comparison](prospective-chi-guard.md) executes
all declared arms without choosing a preferred numerical profile. All15 new guard
starts meet the gradient criterion, but five selected C/guard PDBs stay identical
and do not outperform the synthetic baseline on matched masks. Refusals,
shortfalls, unequal coverage and the stock-supported accuracy counterexample remain.
