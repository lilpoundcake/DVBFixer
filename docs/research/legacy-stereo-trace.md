# Legacy phase trace: endpoints do not replay, so the THR2 crossing is not localized

Status: **negative logging-only mechanism probe**, not an exact replay, causal
attribution, repair or production qualification. Follows the
[full legacy comparator](legacy-comparator.md).

## Frozen observer experiment

Before execution, `research_legacy_stereo_trace.py` enrolled **all six retained
prepared inputs** with their exact `.dat` files: four THR cases and two VAL
controls. Each ran once through unchanged production minimize with
`--no-solvent --atom-naming standard`, default physics/platform and no new seeds,
starting angles, tolerances, forces or stereochemical repair.

A supported Python `Simulation.minimizeEnergy` wrapper records all actual Context
atom coordinates before and after each original native call, forwarding the
original arguments and return untouched. This is not a Python C++ virtual callback
or full iteration trace. Inputs, production source hashes and observer code are
bound before execution; failures and changed endpoints must remain retained.

The independent auditor reads saved snapshots only, uses direct determinants and
RDKit target CIP, verifies every snapshot hash and complete phase inventory, and
compares final real outputs to the original comparator. No native reference,
accuracy metric or reconstruction library is read during this probe.

## Actual results, not the hoped-for replay

[Evidence](legacy-stereo-trace-evidence.json): **6/6 executions completed**, exactly
two native minimization calls each: **12 calls / 24 full Context snapshots**.
Observers restored; original prepared inputs and `.dat` files remain unchanged.

**Zero of six final PDBs repeat the original comparator byte-for-byte.** Maximum
heavy-coordinate differences from the corresponding original endpoints are:

| Case | Maximum endpoint difference, Å | Observed outcome |
|---|---:|---|
| 1ENH THR6 | 0.057236 | CB remains S throughout sampled boundaries |
| 1ENH THR27 | 0.047170 | CB remains R |
| 1CRN THR1 | 0.903382 | CB remains R; Cα becomes nearly planar |
| 1CRN THR2 | 2.130130 | CB remains R, original R→S endpoint failure not replayed |
| 1CRN VAL8 | 0.278880 | Target Cα remains S |
| 1CRN VAL15 | 0.133525 | Target Cα remains S |

Before the first native call the heavy-atom identities match the prepared inputs;
Context versus input-PDB differences are about **1.28–3.34×10⁻⁶ Å**. The retained
full H-containing contexts remain available for further controlled work.
This probe does not distinguish effects of H reconstruction, platform numerical
behavior, or observer-related numerical perturbation. In particular, it is **not
evidence that the observer preserved the original optimization trajectory**.
No repeat-until-failure or post-hoc forced Reference replay was performed.

### THR2: no localized crossing in this run

Its prepared β signed volume is **+0.198390 Å³**: positive, but small. The phase-1
after-state is **+2.432270 Å³**, phase-2 after-state **+2.440301 Å³**, rounded final
**+2.437493 Å³**, all with target CB R. The original full-comparator output with
CB S remains unchanged and independently verified. These are different retained
outcomes, **not** a correction or disproof of the original failure.

The observer therefore **does not locate the original THR2 R→S transition**.
Endpoint snapshots also cannot exclude a transient crossing inside a native call.

### THR1: a distinct near-planar Cα result

Prepared Cα volume is **+0.403414 Å³**, then **+0.080152 Å³** after phase 1 and
**+0.080089 Å³** after phase 2. Rounded final volume is **+0.082859 Å³**. RDKit
assigns CB R but **does not assign target Cα CIP** at those post-phase/final states.
The production Cα wall uses a fraction of the initial signed volume; this outcome
is consistent with its relative-volume design, not evidence of an absolute
tetrahedral-quality guarantee. No force-isolation experiment establishes causality
here. The volume stays positive: this is **not a negative-volume D-Cα finding**.

No sampled Context boundary has a negative protein Cα volume. That is a boundary
measurement only, not a whole-trajectory no-inversion certificate. No scientific
geometry threshold is relaxed to call the near-planar result acceptable.

## Limits and continuation

The useful outcome is the explicit **failure of exact endpoint replay**, plus
retained phase evidence and a new near-planar endpoint. Neither a preferred
optimizer setting nor an observation-preserving replacement is qualified.
All original comparator results, failures and denominators remain unchanged.

Future mechanism work needs separately declared control of the starting complete
state/system and numerical platform, with explicit H-placement/context evidence.
Do not infer the cause from this one nonexact replay, silently change production
platforms, or extend an adaptive search until a desired failure appears.

Private plans/runs are under `tmp/pdbfixer-research/`:
`legacy-stereo-trace-plan-01.json`, `legacy-stereo-trace-run-01/` and
`legacy-stereo-trace-audit-01.json`. Independent reproduction:

```bash
PYTHONPATH=src:. python scripts/audit_legacy_stereo_trace.py \
  --plan tmp/pdbfixer-research/legacy-stereo-trace-plan-01.json \
  --run tmp/pdbfixer-research/legacy-stereo-trace-run-01 \
  --output <fresh-evidence.json>
```

No production/default/dependency/library/public-capability changes or commit/push.
This native-free probe does not use PyRosetta tables; the separate research-only
library/production approval restrictions remain unchanged.

Verification: **803 focused tests passed in 26.26 s**. Fresh saved-only evidence
readback is byte-identical. This is not full repository/platform acceptance.

Subsequent [complete starting-state controls](legacy-start-controls.md) separate
H-start differences from OpenCL repeat variability using saved systems/checkpoints
and bounded Reference controls. They still do not reproduce the original THR2
CB inversion or establish its complete causal mechanism.
