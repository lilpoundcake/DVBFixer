# Legacy starting-state controls: H variability and OpenCL repeat variability are distinct

Status: **bounded native-free mechanism controls**, not a production change or a
complete causal explanation of the original THR2 Cβ inversion. Continues the
[nonexact phase replay](legacy-stereo-trace.md).

## Protocol frozen before execution

`research_legacy_start_controls.py` enrolled all **six original prepared inputs**
and their `.dat` files with two declared labels, seeds **1 and 2**. Python's
random state is saved, seeded only around native `Modeller.addHydrogens`, then
restored. No other optimizer settings, forces, observations or platforms change.
The unchanged production CLI still uses `minimize --no-solvent --atom-naming
standard`, default physics and automatic platform selection.

The installed OpenMM `addHydrogens` source explicitly draws `random.random()` for
new H positions, then minimizes them in an internal Context. The source hash is
bound in the plan. A Python seed controls those draws, **not every possible source
of numerical variability** in the internal H-placement minimization.

At every original native minimization boundary, retain the full atom map/bonds,
System, Integrator and State XML, same-platform checkpoint and actual platform
properties. No native reference or accuracy evaluation is used.

The same initial protocol also declares a separate **16-call phase-1 replay**:

- two previously identified cases: **1CRN THR2** (original endpoint failure) and
  **1ENH THR27** (comparison control), selected before these controls execute;
- both seed labels, two repeats per mode;
- captured-platform checkpoint restart versus Reference XML-state restart;
- identical captured system and native tolerance/iteration budget, no additional
  seeds, retry-until-failure, parameter search or preferred-mode selection.

Reference is an explicit research numerical control, **not a production platform
change**. These are previously studied inputs, not independent scientific samples.

## Capture: 12 completed workflows, but only eight H interventions

[Capture evidence](legacy-start-capture-evidence.json): **12/12 workflows**,
**24 native phase calls / 48 boundary snapshots**. Each boundary also has XML and
checkpoint evidence. All actual platforms were **OpenCL, Apple M3 Pro, single
precision**. Hooks and Python random states restored; original prepared/dat files
remain intact.

The cases follow different existing production branches:

| Inputs | H calls across both labels | Interpretation |
|---|---:|---|
| Two 1ENH masks | 0 | Keep-H fast path: seed labels do not intervene; these are numerical-repeat controls |
| Four 1CRN masks | 8 | H is stripped/rebuilt through the existing variant-aware path |

Thus it would be incorrect to call all 12 runs independent seeded H experiments.
For each 1ENH mask, both complete phase-1 starting coordinate arrays are identical,
including H. For each 1CRN mask, starting heavy coordinates match exactly between
the labels, but H coordinates differ: exact-name H RMS differences about
**1.068–1.084 Å**, maxima **2.022–2.048 Å**.

These are **exact-name metrics, not symmetry-minimized physical H displacement**;
chemically equivalent methyl H permutations contribute. For example, THR2's
unique hydroxyl H `HG1` differs by approximately **0.134544 Å**, whereas named
methyl H differences are much larger. Do not interpret every named-H difference
as an equally large physical rearrangement.

All THR2 post-phase boundaries in this capture have CB R. THR6 remains CB S, as it
already was in the prepared source. The earlier THR1 near-planar Cα endpoint does
not recur in these two labelled runs. No sampled boundary has a negative protein
Cα volume; this is not a guarantee about every point inside a native call.

### Explicit audit correction: tiny heavy-coordinate changes are not zero

The first saved-only auditor incorrectly required bit-exact heavy coordinates
across `addHydrogens` and failed. Actual numerical displacement is about
**1.283×10⁻⁶ Å** in the OpenCL H-placement path. Original auditor source and failure
log remain retained. A hash-bound reporting amendment removes that incorrect
assumption and reports `exact: false` and the measured displacement verbatim.

No scientific run, seed, geometry threshold, optimizer setting or coordinate was
changed. This is **not an observation-preservation pass obtained by widening a
gate**: the report makes no such pass claim. Heavy coordinates remain identical
*between the two seed-labelled starts*, a different comparison.

## Saved complete-start replay: 16/16 executed

[Replay evidence](legacy-start-replay-evidence.json): every replay's initial
heavy **and H coordinates exactly match** the captured phase-1 starting array.
Same-platform runs load the original checkpoint with recorded OpenCL properties;
Reference runs load the identical serialized system/integrator/state. All use
native tolerance **10 kJ mol⁻¹ nm⁻¹**, maximum **1,000 iterations**.

| Repeated starting state | OpenCL maximum heavy endpoint difference, Å | Reference difference |
|---|---:|---:|
| THR2 / seed 1 | 0.005183 | 0 |
| THR2 / seed 2 | 0.023293 | 0 |
| THR27 / label 1 | 0.030307 | 0 |
| THR27 / label 2 | 0.005301 | 0 |

All four declared OpenCL repeat pairs differ, despite identical full starting
coordinates/checkpoints. All four Reference pairs repeat their endpoints exactly
for heavy atoms and H. The two THR27 labels share the same coordinate start, so
these are **four declared pairs, not four independent starting structures**.
Cross-platform endpoints also differ; their maximum heavy differences range
**0.001981–0.079692 Å** for the declared first-repeat comparisons.

This establishes numerical repeat variability on the tested captured OpenCL route
separately from changed H starting coordinates. It does not identify every kernel
cause, prove general Reference determinism, or make Reference a superior scientific
profile. No mode is selected by native accuracy or silently adopted in production.

**All sixteen replay endpoints keep target CB R.** Therefore these controls do
not reproduce or fully explain the original THR2 R→S event. The original wrong-CB
output remains immutable and valid evidence; a negative replay does not erase it.

## What the captured force inventory establishes

The serialized systems contain the existing Cα signed-volume guard. Its relevant
particle tuple is **N–CA–C–CB**, centred at **CA**, not CB. No CB-centred periodic
improper is present for the tested target residues. The listed bond/angle/proper
torsion/nonbonded terms are **not absent**, and may still influence stereochemical
energetics; this inventory is not a force-isolation experiment.

The narrow conclusion is that the current explicit Cα safety invariant does **not
also implement a signed-volume Cβ invariant**. Together with the original THR2
prepared β volume of only +0.198390 Å³, this identifies an unprotected, poorly
shaped starting centre worth investigating. It is not yet proof of a unique
causal path or a qualified production fix.

## Retention and continuation

Private prefix `tmp/pdbfixer-research/legacy-start-controls-` contains the initial
plan, 12-run capture, separately bound replay plan and 16-run replay. Independent
audits are `legacy-start-capture-audit-02.json` and
`legacy-start-replay-audit-02.json`; original auditor/failure and reporting
amendment are retained separately. No original source, failed output or historical
denominator was replaced. No raw scientific runtime data are promoted to public
completion capability.

Independent saved-only reproduction:

```bash
PYTHONPATH=src:. python scripts/audit_legacy_start_controls.py capture \
  --plan tmp/pdbfixer-research/legacy-start-controls-plan-01.json \
  --run tmp/pdbfixer-research/legacy-start-controls-capture-01 \
  --output <fresh-capture-evidence.json>
PYTHONPATH=src:. python scripts/audit_legacy_start_controls.py replay \
  --plan tmp/pdbfixer-research/legacy-start-controls-replay-plan-01.json \
  --run tmp/pdbfixer-research/legacy-start-controls-replay-01 \
  --output <fresh-replay-evidence.json>
```

Next work must keep start-state variability, native minimization repeat variability
and missing Cβ protection distinct. Do not extend a seed/platform grid until a
preferred failure appears, restore observations, swap distinct THR branches, or
claim a repair without an independent final stereo audit. No production/default/
dependency/library change, commit or push was made.

Verification: **810 focused tests passed in 28.00 s**; targeted Ruff and the
17-file completion-package mypy check pass. Both independent saved-only evidence
readbacks are byte-identical. This is not full repository/platform acceptance.

A subsequent [private Cβ wall/stress experiment](beta-wall.md) checks an explicit
target-THR intervention without new seeds or production changes. A finite wall
preserves R in the declared stress endpoints, but only two of six stress-plus-wall
endpoints pass geometry; the other four are refused, not presented as repairs.
