# Native χ checkpoints: pre-step gradients differ from returned endpoint gradients

Status: **private native-free numerical-mechanism evidence**, not a production
fix, full function/iteration/return-branch trace or convergence certificate.
Follows [unit-aware endpoint derivatives](chi-derivatives.md); the
[research-only licensing and production gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03)
is unchanged.

## Frozen supported diagnostic route

The earlier Python virtual-callback feasibility controls produced no usable
trace and are retained separately. This experiment instead uses the supported
native **`AtomTreeMinimizer.deriv_check_result()`** API. It enrolls **all 42
original mode-C starts**, all 15 pools and seven inherited nonexecutions, with
no endpoint-gradient/native/geometry filtering.

The same target-only MoveMap, ScoreFunction, numerical options and pre-score are
used, following the pinned installed-revision `MinMover.cc` orchestration. Enable
only native derivative checking and disable diagnostic stdout in the first run.
The native checker uses its fixed five signed increments **0.001/0.002/0.004/
0.008/0.016°**; no new adaptive step grid or optimizer settings are selected.
All initial snapshots, before/after atom traces and complete final numerical
results must be compared byte-for-byte with mode C. Changed states would remain
retained and would not qualify as exact original-return evidence.

[Norm evidence](chi-checkpoints-evidence.json): **42/42 exact endpoint replays**,
zero execution blockers, **1,097 derivative-evaluation checkpoint records**, and
**0.0 Å protected movement**. This is repeated numerical evidence on existing
inputs, not 42 independent scientific samples. Parent geometry (40/42) and all
previous outputs/scores/selectors remain unchanged.

## Sparse native records are not zero gradients

`NumDerivCheckData` always retains the analytic **L2 norm**, but individual
`DOF_DataPoint` rows are populated only when the native analytic or numeric
derivative exceeds **0.001** in magnitude. Invalid/default rows are absent data,
not evidence of zero derivatives. Moreover, L2 norm does not equal the largest
component for a two-χ target: `norm/sqrt(n) ≤ max(abs(g)) ≤ norm`.

For THR22 the last norm was 0.0001033724 per degree: that alone could not establish
whether either component exceeded the 0.0001 cutoff. The first run therefore
did **not** claim to explain THR22 from its norm alone.

## Separately frozen verbose component amendment

The second run keeps all 42 starts and adds native verbose/stdout diagnostics.
All **16 exposed numerical option getters** are checked before/after constructing
the verbose options; they must match exactly. Only diagnostic flags differ.
Checked native-stdio flushing and labelled synthetic real-tracer BEGIN/END markers
bound every captured log. There are no Python overrides of C++ virtual callbacks.

[Component evidence](chi-checkpoint-components-evidence.json): again **42/42
exact endpoint replays**, all numerical checkpoint files identical to the quiet
run, **1,097 complete checkpoint vectors / 1,957 printed component values**.
Each printed vector's L2 norm is checked against the full-precision native record.
Text-format tolerance is separate from the scientific derivative or geometry gates;
cutoff classifications include an explicit print-precision uncertainty margin.
Both adapters restore; all logs and raw native data remain retained.

The quiet wrapper's `checkpoints.json` records its requested quiet diagnostic
flags; in this amendment the **actual native override flags** and numerical option
fingerprints are bound separately in `verbose-adapter.json`. Do not mistake the
outer request for the applied native stdout/verbose flags.

## Last computed gradient versus returned endpoint

Values below are score units **per degree**, not per radian:

| Case/start | Last computed gradient components | Returned endpoint maximum | Explicit cap |
|---|---|---:|---|
| 2CI2 THR22/2 | −9.55821122e−5, −3.93689471e−5 | 0.001844293 | No |
| 1CTF VAL56/0 | 5.66662531e−6 | 0.007424845 | No |
| 1R69 THR27/1 | 0.00029656222, −1.89785877e−5 | 0.000296561 | Yes |
| 1CTF THR76/2 | −0.000951025883, −7.65937529e−5 | 0.000951026 | Yes |

Thus both previously unexplained no-warning endpoints have **last computed
components below gmax**, despite a larger gradient at the exactly reproduced
returned endpoint. THR22 has 20 recorded derivative evaluations; VAL56 has 9.
The two explicit cap cases each have 201 derivative checkpoints, and their last
component maxima remain above the cutoff.

The pinned `Minimizer.cc` computes `Gmax` before the line-minimization call, then
tests that value after the call without necessarily recomputing derivatives at
the newly returned point (lines 674–744 in the retained revision). These source
and checkpoint observations support **pre-step-gradient stopping** for the two
no-warning cases; they are not evidence of missing gradients, a unit error in
Rosetta itself, or observed-coordinate movement.

This is still **not a directly recorded return-branch identifier or full function
call trace**. No claim of global minimization, tolerance-independent coordinates,
accurate native reconstruction or a generally validated optimizer follows.

## Reproduction and remaining scope

Private plans/runs: `tmp/pdbfixer-research/chi-checkpoints-plan-01.json`,
`chi-checkpoints-run-01/`, explicit `chi-checkpoints-plan-02.json` and
`chi-checkpoints-run-02/`. Native source/provenance is pinned in
`minmover-source-1a56185/` and the earlier AtomTree source directories. No library
tables or raw rotamer parameters are published.

Independent auditors read saved files only: `audit_chi_checkpoints.py` validates
all exact-versus-changed replays, native sparse/norm arithmetic, ownership and
settings; `audit_chi_checkpoint_components.py` binds verbose capture controls,
actual native flags/options, complete vectors and source/norm-report hashes.

```bash
PYTHONPATH=src:. python scripts/audit_chi_checkpoints.py \
  --plan tmp/pdbfixer-research/chi-checkpoints-plan-02.json \
  --run tmp/pdbfixer-research/chi-checkpoints-run-02 \
  --output <fresh-verbose-norm-report.json>
PYTHONPATH=src:. python scripts/audit_chi_checkpoint_components.py \
  --plan tmp/pdbfixer-research/chi-checkpoints-plan-02.json \
  --run tmp/pdbfixer-research/chi-checkpoints-run-02 \
  --norm-audit <fresh-verbose-norm-report.json> --output <fresh-components.json>
```

**774 focused tests passed in 26.71 s**; not full repository/platform acceptance.
No production defaults/consumers/dependencies, fixtures or installed library
changed; no commit/push. An actual-endpoint-gradient guard would require a separate
frozen research experiment, not a silent modification of this replay. Original
tables/provenance/license rights and explicit user approval remain mandatory
before production, independently of scientific qualification.

Subsequent [bounded actual-endpoint guard](chi-endpoint-guard.md) preserves the
same inner numerical settings and caps the experiment at three total calls per
start. It reports39/42 measured-gradient passes, three retained exhaustions and
37/42 combined geometry/gradient passes; no production or accuracy qualification.
