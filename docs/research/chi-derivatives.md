# Unit-aware analytic χ derivatives: four endpoints exceed mode-C gmax

Status: **private native-free numerical diagnosis**, not a production profile,
full derivative qualification or convergence certificate. Follows
[logging-only termination replay](chi-termination.md) and
[finite local probes](chi-stationarity.md). The
[research-only licensing/production gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03)
remains unchanged.

## The missing unit distinction

The pinned `MinimizerMap.cc` at installed revision `1a56185` multiplies atom-tree
PHI (radians) by **180/π** when forming minimizer variables, and divides by the
same factor when copying variables back. `atom_tree_minimize.cc` divides the
torsion derivative by that scale too. Thus the target-χ minimizer variables here
are **degrees** and their gradients are score units **per degree**.

The earlier finite probes correctly reported score units **per radian**, but those
values cannot be compared directly to mode-C `gmax=1e−4`. In per-radian units that
cutoff is about **0.00572958**, not 0.0001. Reporting a resolved finite descent at
the previous 1e−7 energy floor is not the same as exceeding this gradient cutoff.
No earlier raw result, finite-descent count or numerical setting is changed.

## Frozen experiment

`research_chi_derivatives.py plan` froze **all 42 saved mode-C endpoints**, all
15 pools and seven inherited nonexecutions, including both parent geometry
rejections. No native/gradient/motion filtering. Original full snapshots/scores
replay exactly; saved endpoints materialize by moving only saved mobile atoms,
with complete endpoint coordinates exact and the existing energy-recompute bound.

On a separate clone, set up the same target-only MoveMap, `MinimizerMap` and
`ScoreFunction.setup_for_minimizing`, then call **`AtomTreeMultifunc.dfunc`**.
No `Minimizer.run`, packing, optimization, parameter change or new selector.
Validate every node's target/CHI/PHI identity, unique χ mapping and 180/π scaling.
Both signed finite steps **0.01° and 0.005°** are applied to minimizer variables
from the endpoint vector; all raw atom traces and numerical results are retained.

Compare analytic gradients with both central differences and with the older
`Pose.set_chi` finite slopes after unit conversion. The predeclared diagnostic
bound is **1e−6 + 1e−3 × max(|analytic|, |numeric|)** per degree. It is a reporting
criterion, not a scientific geometry gate; all mismatches remain visible and
the bound is not fitted afterward.

Pinned public source/provenance is retained under
`tmp/pdbfixer-research/atomtree-source-1a56185/` and
`atomtree-derivative-source-1a56185/`. No rotamer tables were downloaded.

## Results

[Independent saved-only evidence](chi-derivatives-evidence.json):

| Measure | Result |
|---|---:|
| Endpoints enrolled/audited; execution blockers | 42/42; 0 |
| Analytic χ derivatives; signed probes | 62; 248 |
| Analytic-versus-numeric comparisons | 124 |
| Outside the predeclared reporting bound | **1/124** |
| Minimizer-variable versus external-χ numeric mismatches | **0/124** |
| Protected maximum motion | **0.0 Å** |
| Endpoints with analytic maximum above 1e−4 per degree | **4/42** |

The four above-cutoff endpoints:

| Case/start | Maximum | Explicit max200 warning |
|---|---:|---|
| 2CI2 THR22 / 2 | 0.001844293 | No |
| 1R69 THR27 / 1 | 0.000296561 | Yes |
| 1CTF VAL56 / 0 | 0.007424845 | No |
| 1CTF THR76 / 2 | 0.000951026 | Yes |

Units in this table are score units per degree. **15 of the 19 finite-descent
endpoints have analytic gradients at or below the declared cutoff**; their tiny
remaining descent is compatible with a nonzero stopping tolerance. The other
four remain above it: two explicitly reached the iteration cap, and the return
mechanism of **THR22/2 and VAL56/0 remains unresolved**. Endpoint gradients alone
do not prove which branch returned in the original L-BFGS execution.

The one derivative-comparison mismatch is **1R69 VAL24/start 0, χ1**, at 0.01°:
analytic 5.9069302e−5, numeric 5.7968805e−5 per degree; error 1.1004971e−6 versus
bound 1.0590693e−6. At the separately predeclared 0.005° step, error is 2.7511073e−7
and passes the same bound. This step sensitivity is consistent with finite-step
effects but does not establish their cause or justify dropping the coarse-step
failure. No new step grid or relaxed bound is introduced.

## Limits and reproduction

This is unit-aware consistency evidence on previously evaluated inputs. It does
not validate every energy term/chemistry, prove an exact stationary point or
global minimum, qualify probe geometry/chirality, or establish native accuracy.
Parent mode-C geometry (40/42), coordinates, scores and selectors stay unchanged.
No derivative probe is selected as a scientific output.

Private plan/run: `tmp/pdbfixer-research/chi-derivatives-plan-01.json`,
`chi-derivatives-run-01/`. The independent auditor invokes no PyRosetta or score
evaluation; it checks source/state/vector bindings, all probe atom traces,
ownership and units, computes both numerical comparisons, and reports all failures.
No raw endpoint angles, library parameters or tables enter the publication.

```bash
PYTHONPATH=src:. python scripts/audit_chi_derivatives.py \
  --plan tmp/pdbfixer-research/chi-derivatives-plan-01.json \
  --run tmp/pdbfixer-research/chi-derivatives-run-01 \
  --output <fresh-report.json>
```

**761 focused tests passed in 26.76 s**; not full repository/platform acceptance.
Production consumers/defaults/dependencies, fixtures and installed libraries are
unchanged; no commit/push. Original tables/provenance/appropriate license rights
and explicit user approval remain mandatory before production, independently of
scientific qualification. Exact return-path tracing remains a separate experiment.

Subsequent [supported native derivative checkpoints](chi-checkpoints.md) reproduce
all42 starts exactly and show that the two no-warning cases had last computed
components below the cutoff before returning endpoints with larger gradients.
This supports the pinned source's pre-step-gradient stopping mechanism without
claiming a full function/iteration or directly recorded return-branch trace.
