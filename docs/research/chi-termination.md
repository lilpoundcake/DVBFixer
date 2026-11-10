# Logging-only χ replay: two explicit iteration-limit stops

Status: **native-free diagnostic replay**, not a changed numerical profile or
convergence certificate. Continues the [finite local probes](chi-stationarity.md)
and [stopping controls](chi-stopping-controls.md). Research-only PyRosetta and
the [mandatory production/license gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03)
remain unchanged.

## Frozen logging-only intervention

`research_chi_termination.py plan` binds **all 42 original mode-C starts**, all
15 executable pools and seven inherited nonexecutions. No filtering by finite
descent, geometry, native error or score margin. Numerical settings, objective,
weights, library, atom ownership and initial states remain unchanged.

The sole Rosetta option change is `-unmute core.optimization.Minimizer`.
Require exact equality of original complete import snapshots, before/after atom
trace bytes **and every final result including numerical scores** against mode C.
This is a real minimizer replay, unlike the previous saved-coordinate finite-probe
materialization. It is repeated numerical evidence, not 42 new scientific samples.

## Initial capture failure and explicit amendment

The first run reproduced all 42 endpoints exactly, but native buffering defeated
per-apply file-descriptor capture: **two max-cycle warnings appeared after those
captures had ended**. Its zero per-case counters are not trustworthy evidence
that no warnings occurred. That complete attempt and its unscoped outer log remain
retained; no warning is assigned to a case from that attempt.

A synthetic diagnostic initially failed because Python bindings do not implement
`TracerProxy << str`; that failure is also retained. A second synthetic control
confirmed that a real tracer write plus native `fflush` reaches the active fd
capture. No optimizer-warning text is manufactured by this control.

A separately frozen amendment, `research_chi_termination_flush.py`, adds checked
native-stdio flushes before/after each fd boundary and **synthetic tracer begin/end
markers** in every capture. A missing marker or failed flush is an execution
blocker, never an absence-of-warning result. Both process-local adapters restore.
The markers are labelled synthetic and distinct from actual optimizer warnings.
No numerical setting or scientific acceptance threshold changes.

## Results

[Independent evidence](chi-termination-evidence.json):

| Measure | Result |
|---|---:|
| Controlled captures / exact mode-C endpoint repeats | **42 / 42** |
| Explicit `LBFGS MAX CYCLES 200 ... NOT CONVERGED` warnings | **2** |
| Explicit line-search-abort messages | **0** |
| Cases without either explicit message | **40** |

The max-cycle warnings belong to:

- **1R69 THR27, start index 1**;
- **1CTF THR76, start index 2**.

Both also have resolved finite descent in the earlier probe study. This is a
saved-only association over the complete common inventory, not a new enrollment
rule or proof of why their residual slopes have those values.

The two warnings provide explicit evidence that these replays reached the
iteration cap without the reported convergence test succeeding. They do **not**
explain the other **17/19 finite-descent endpoints**, nor do the 40 cases without
explicit messages certify normal convergence. The pinned L-BFGS implementation
contains return paths that need not emit either of these messages; no full
analytic-gradient or per-iteration trace has been collected here.

No extra iterations, tighter tolerances, new candidates or native RMSD reporting
were added. Parent mode-C geometry (40/42), all saved coordinates and all earlier
selectors/results remain unchanged. No production refiner is enabled.

## Reproduction and checks

Private retained artifacts under `tmp/pdbfixer-research/`:
`chi-termination-plan-01.json`, failed-buffer-attribution `chi-termination-run-01/`
and its outer log, synthetic diagnostics, explicit amended
`chi-termination-plan-02.json`, and controlled `chi-termination-run-02/`.
Individual `minimizer.log` files retain controls and genuine messages.

`audit_chi_termination.py` reads saved files only, independently verifies source
bindings, all four exact replay artifacts per case, native boundary controls,
log hashes, frozen runtime settings and explicit diagnostic text. It refuses
per-case attribution for the unamended buffered attempt. The finite-probe join
reads only its native-free report, not references or accuracy files.

```bash
PYTHONPATH=src:. python scripts/audit_chi_termination.py \
  --plan tmp/pdbfixer-research/chi-termination-plan-02.json \
  --run tmp/pdbfixer-research/chi-termination-run-02 \
  --output <fresh-report.json>
```

**756 focused tests passed in 26.37 s**; not full repository/platform acceptance.
No source fixtures, installed library, production consumers, defaults, dependencies
or public capability changed. No table parameters are published; no commit/push.
Original tables, documented suitable rights/provenance and explicit user approval
remain mandatory before production, independently of scientific qualification.

Subsequent [unit-aware analytic derivative checks](chi-derivatives.md) show that
15 of the 19 finite-descent endpoints already meet the degree-scaled gradient
cutoff. Four remain above it, including the two explicit cap events; return-path
diagnosis is still needed for the other two. This does not alter the logs above.
