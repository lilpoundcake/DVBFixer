# Mean-proposal context checks and residue-pair score attribution

Status: **finite research diagnosis, no new selector or production capability**.
Continues the [library-mean study](dunbrack-means.md) with all retained pools,
then both prior-versus-score disagreements. No native-error-based enrollment,
weight refit, candidate generation, minimization, packing or hydrogen optimization.
Installed PyRosetta remains research-only under the
[mandatory production gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03).

## All saved pools: context matches across arms

`scripts/audit_mean_comparability.py` binds all mean pools and corresponding parent
relative-rotation pools before readback. The
[comparability evidence](mean-context-comparability-evidence.json) checks **84
saved scored snapshots**, 42 mean and 42 parent states, covering all 15 executable
pools. All 22 original masks remain: seven inherited refusals/context blockers
are nonexecutions, not additional successful tests.

In **15/15 pools**, exact source PDB bytes, generated-atom ownership and target
identity match between arms. Across every feasible state in each pooled comparison,
observed/scorer-extra heavy coordinates, non-target hydrogen coordinates, virtual
coordinates, residue/atom type inventories and external/target graph descriptions
match exactly. Target generated heavy atoms and target H are intentionally excluded
from context equality. This closes the earlier cross-arm comparability caveat on
these states; it is not proof of correct context protonation or experimental truth.

## VAL equivalence reduces, but does not remove, the accuracy failure

Saved reporting already includes generated-only CG1/CG2 permutation for whole
VAL masks. No atom renaming, native-based selection or new score is performed.
For the four executable new-input cases, matched RMSD medians (Å):

| Mean-pool selector | Exact atom names | Generated-only VAL symmetry |
|---|---:|---:|
| Family prior | 1.210568 | 0.879635 |
| Fixed-heavy score | 1.930166 | 1.463808 |

For 1PGB VAL29 the erroneous score-selected result drops from 1.770814 to
**1.499963 Å** with methyl equivalence, still much worse than the prior-selected
**0.078348 Å**. VAL21 drops from 2.336914 to 1.427653 Å, also not a label-only
error. Old-eleven medians remain 0.198491 Å under either reporting convention.
These are retrospective metrics on the same four/eleven masks, not new validation.

## All two disagreements: exact repeated states before EnergyGraph accounting

`scripts/research_mean_pair_energies.py` enrolled **all nonempty prior/score
disagreements**, 1PGB VAL29 and 1R69 THR26, without reading native coordinates.
Before interpreting each state's residue-pair graph, the unchanged fixed-heavy
scorer repeats its original import/score: all four complete compressed pose
snapshots and numerical scores/terms must match exactly. Graph readout then must
leave every heavy, H and virtual atom coordinate unchanged.

The original two attempts stopped on JSON container equality at the 1R69 terminal
OXT inventory: the live API returned a tuple while the saved JSON contained a list.
Both failed runs and the instrumented diagnostic remain retained. A separately
frozen v2 driver, `research_mean_pair_replay.py`, normalizes JSON container types
only; it adds no score/coordinate/snapshot tolerance and changes no numerical value.
All four exact-state checks pass under that explicit transport correction.

The [independent pair evidence](mean-pair-energies-evidence.json) reads saved graph
data without importing PyRosetta or evaluating scores. All edge sums for `fa_atr`,
`fa_rep` and `fa_sol`, plus `fa_dun` onebody sums, close against their original
weighted totals. The predeclared bookkeeping bound is
`1e-5 + 1e-7 * sum(abs(terms))`; actual pair-delta discrepancies are about 1e−13.
This bound concerns cached-energy accounting, not exact original-state replay or
scientific geometry acceptance.

**For the three audited pair terms, no non-target residue pair changes in either
case.** Each `fa_dun` difference
belongs solely to the target residue. Thus these recorded differences are not
caused by a changed distant context, mismatched atom types or untracked context-H
motion between these exact states.

## Where the model preferences arise

All deltas below are score-selected minus prior-selected, in weighted Rosetta
score units. They explain model arithmetic, not physical/native causality.

### 1PGB VAL29: extra attraction/solvation outweigh repulsion and prior

Total selected-minus-prior score is −0.404808. Pair totals remain the earlier
`fa_atr` −1.597526, `fa_sol` −0.867703 and `fa_rep` +0.790042, while target
`fa_dun` worsens by +1.012951.

Largest attraction gains involve **PHE30 (−0.597736), ALA26 (−0.592192) and
THR18 (−0.482900)**. The same neighbors account for major solvation gains and
repulsion increases. The bad native-accuracy choice is favored by this balance,
not by an absent empirical-prior term or a methyl-name-only discrepancy.

### 1R69 THR26: context repulsion relief overrides the family prior usefully

Total score difference is −9.906042; pair `fa_rep` changes by −9.448179, dominated
by **VAL24 (−5.065842), PHE44 (−3.457978) and ARG43 (−1.470374)**. The target
`fa_dun` worsens by +1.023393. Here the context correction favors the more accurate
mean proposal (0.240256 Å versus the prior's 2.022584 Å).

These counterexamples still oppose a blanket prior-upweighting recipe. They do
not establish which atom-level/H interactions or context/model errors explain
the native differences; residue-pair attribution cannot answer that alone.

Focused suite: **731 passed in 26.12 s**, including synthetic context-ownership,
complete-disagreement enrollment, accounting and exact JSON transport controls.
This is not a full repository or cross-platform acceptance run.

## Reproduction and limits

Private frozen artifacts under `tmp/pdbfixer-research/`:
`mean-comparability-plan-01.json`, `mean-pair-energies-plan-01.json`, failed
`mean-pair-energies-run-01/` and `-02/`, explicit
`mean-pair-energies-plan-02.json`, completed `mean-pair-energies-run-03/`.
Instrumented and single-state diagnostic attempts remain separate from the four
reported final repeats. No raw library/table parameters are published.

```bash
PYTHONPATH=src:. python scripts/audit_mean_comparability.py audit \
  --plan tmp/pdbfixer-research/mean-comparability-plan-01.json \
  --output <fresh-comparability.json>
PYTHONPATH=src:. python scripts/audit_mean_pair_energies.py \
  --plan tmp/pdbfixer-research/mean-pair-energies-plan-02.json \
  --run tmp/pdbfixer-research/mean-pair-energies-run-03 \
  --output <fresh-pairs.json>
```

No source fixture, scientific threshold, installed library, production consumer,
dependency or public capability changed; no commit/push. Context completion,
generated-only relaxation and independent accuracy qualification remain separate
unexecuted gates, not conclusions of this saved-state diagnosis.
