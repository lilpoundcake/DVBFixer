# Target-χ stopping controls: bounded results are numerically sensitive

Status: **native-free numerical controls on already evaluated research inputs**,
not accuracy qualification, optimizer certification or a production profile.
Continues the [bounded χ-relaxation experiment](generated-chi-relax.md). No new
starting conformations, objective terms, weights, context repair or packing.
The [research-only/production approval gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03)
remains unchanged.

## Why another control, and what was frozen

The first bounded run moved many starts very little and explicitly did not claim
convergence. The upstream [Minimizer.cc at the reported PyRosetta revision](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/optimization/Minimizer.cc)
dispatches `lbfgs_armijo_nonmonotone` to fractional energy-change checking and
`lbfgs_armijo_nonmonotone_atol` to absolute checking. Its L-BFGS path also checks
a maximum-gradient cutoff. The installed default cutoff is verified as 1.0
before each control overrides it. Source bytes/provenance are pinned privately
under `tmp/pdbfixer-research/minimizer-source-1a56185/`; no table data are downloaded.

Before any new minimization, `research_chi_stopping_controls.py plan` froze this
ordered sequence on **all 42 original feasible starts in all 15 pools**:

| Mode | Energy criterion | Tolerance | Gradient cutoff | Maximum iterations |
|---|---|---:|---:|---:|
| Original baseline, unchanged | Fractional | 0.01 | 1 | 200 |
| A | Absolute | 0.01 | 1 | 200 |
| B | Absolute | 1e−6 | 1 | 200 |
| C | Absolute | 1e−6 | 1e−4 | 200 |

Baseline→A changes the criterion; A→B changes tolerance; B→C changes gradient
cutoff. This is a sequential sensitivity design, **not a full factorial** or a
search for the best native-RMSD setting. Backbone/jumps/other χ remain disabled;
only generated distal heavy atoms and target sidechain H may move. Generated CB,
all observations, context H/extras/virtuals remain protected. All seven original
unsupported masks remain inherited nonexecutions in each mode, not new samples.

Every start must repeat the original complete pose snapshot and numerical score
exactly. Runtime getter checks verify all six stopping/mobility settings on every
one of the 126 created minimizers; process-local adapters restore. Objective,
installed library hash, source/atom identities and score weights remain fixed.

## Results: motion and choices respond to stopping policy

[Independent native-free evidence](chi-stopping-controls-evidence.json):

| Mode | Starts audited | Common passes | Protected motion | Changed selected indices vs baseline | Changed rounded selected structures vs baseline |
|---|---:|---:|---:|---:|---:|
| A | 42 | 40 | 0.0 Å | 0/15 | 1/15 |
| B | 42 | 40 | 0.0 Å | 1/15 | 12/15 |
| C | 42 | 40 | 0.0 Å | 1/15 | 15/15 |

All 126 additional bounded runs complete without execution/mobility blockers.
All CA remain S; all THR CB remain R; global stereo findings remain empty. There
are **126 repeated starts, not 126 independent scientific samples**. Every mode
selects one independently geometry-accepted result for each of the same 15 pools.

Largest generated-heavy motion rises to 1.553417 Å in A and 2.550827 Å in B/C.
Raw final-score changes relative to the baseline range from −19.370642 to 0 in A,
−29.269935 to 0 in B, and −29.269935 to −2.99e−7 in C. Lower raw score is not
interpreted as better native accuracy. Score semantics still include raw refined
target sidechain H, not fresh rounded-heavy import scores.

B and C change the selected initial index only for **1PGB THR11**, from 0 to 2.
No native reference or accuracy file is read in this suite or its auditor; no
claim about that new choice's RMSD is made. No mode is selected as a preferred
research/production profile based on these outcomes.

## Rejections and remaining numerical uncertainty

All unchanged raw/rounded geometry rejections remain visible:

- A rejects 2CI2 THR22 index 2 (CG2–own backbone C clash) and THR55 index 0
  (OG1–I/52 O clash).
- B and C reject 1CTF VAL56 index 2 (CG2–own backbone C clash) and the same THR55
  index 0. These are not pruned to improve the reported coverage.

Between B and C, selected indices are stable, **but four rounded selected
structures still change** and the maximum per-start raw-score difference is
0.052735 Rosetta units. The latter occurs in a nonselected THR11 start. Therefore
the suite does not certify full convergence, a global minimum or tolerance-
independent coordinates. It establishes that the original bounded outcome was
stopping-policy-sensitive, not that a single numerical rule is physically correct.

Neither the higher-motion candidates nor their lower model scores are used to
adapt geometry cutoffs, weights, starting angles or the declared controls.
Further derivative/iteration/convergence controls require a separate protocol,
not adaptive tuning against already known native errors.

## Reproduction and scope

`audit_chi_stopping_controls.py` invokes no PyRosetta/minimization or score
evaluation. It binds all complete initial states, independently verifies atom
mobility and protected motion, raw/rounded common geometry and chirality,
unchanged score weights and final term totals, and complete-pool selection.
It reports all profiles, with no cross-profile selector and no native accuracy.
No library means, sigmas, raw weights or table rows enter the publication.

Private suite/plans: `tmp/pdbfixer-research/chi-stopping-plans-01/`; runs:
`chi-stopping-<profile-id>-run-01/`. The original relaxation run/negative result
and all previous evidence stay unchanged.

```bash
PYTHONPATH=src:. python scripts/audit_chi_stopping_controls.py \
  --plans tmp/pdbfixer-research/chi-stopping-plans-01 \
  --runs-prefix tmp/pdbfixer-research/chi-stopping- \
  --output <fresh-output.json>
```

**742 focused tests pass in 26.95 s**, not a full repository/platform acceptance
suite. No source fixture, installed library, production code/default/dependency
or public capability changed; no commit/push. Original tables, documented
appropriate license rights and explicit user approval must still be obtained
before any production transition, independently of scientific qualification.

Subsequent [finite local χ probes](chi-stationarity.md) test every saved mode-C
endpoint without optimization or native reporting. They find residual finite
descent at 19/42 endpoints, but do not identify the termination mechanism or
validate analytic derivatives. All stopping-control results above remain unchanged.
