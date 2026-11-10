# Rigid methyl rotation lowers initial energy but does not prevent inward HA

Status: **negative prospective initial-condition mechanism control**, following
the [exact pair attribution](ha-nonbonded-pairs.md). No production change or repair
qualification. In particular, this is not a new explanation of original THR2.

## Frozen matched control

[Complete evidence](methyl-start-control-evidence.json) retains **48/48 executed
cells**, 24 original H starts × `unchanged` / `energy_selected_methyl`, with no
retries, timeouts or budget nonexecutions. These are the existing four masks,
two H seeds and three source-geometry arms on **one protein**, not 24 independent
proteins. All existing heavy atoms remain immobile in the original H-only System.

Only generated THR24 **HG21/HG22/HG23** may change initially, by a common rigid
rotation about the fixed CG2-to-CB axis. Every other initial coordinate remains
exact, including GLN16 HA. Original C–H radii, axial coordinates, H–C–H internal
geometry and all pair distances among the three H are preserved. Source input
ownership excludes observed H. No new H, force, parameter, exclusion, seed or
whole-system minimization is introduced.

Both arms evaluate the fixed **24-angle grid 0–345° in 15° steps**, using original
whole-System *initial potential energy*. The selected arm uses its minimum;
energies within 10⁻⁸ kJ/mol are tied by shortest rotation then signed angle.
No final geometry or native error is used to choose the angle. All grid energies
are retained and independently recomputed with the native Reference force field.
The original **1.0 kJ/mol/nm / maximum 50-iteration / Reference** H minimizer stays
unchanged. Callback reporting does not request an early exit.

Plan: `tmp/pdbfixer-research/methyl-start-plan-01.json`, SHA256
`d52e8c56ca650e04a212737c9f9110169f15063c52de3d0a86d7a62643560461`.
Private grid coordinates and paths: `methyl-start-run-01/` beside the plan.
Budgets were 900 seconds per cell, 3,600 seconds total. Private paths and raw
coordinates/gradients do not transfer through Git.

## Independent acceptance and negative result

All **24 unchanged controls exactly reproduce their original H endpoints**.
The auditor independently verifies signed axis angles, preserved radii/internal
geometry, all nonmethyl starting coordinates, original serialized System and all
grid energies/selection. Existing heavy coordinates are exact at every callback.
Every mobile-H physical force matches the negative native callback gradient:
maximum component difference **0.0 kJ/mol/nm** across 2,400 callbacks, not only HA.
All 2,496 initial/callback/final nodes and final H-bond deviation summaries remain
published. There is no claim of complete H geometry or whole-system acceptance.

| GLN16 HA result | Unchanged | Selected methyl rotation |
|---|---:|---:|
| Cells with an inward callback | 2/24 | 2/24 |
| Cells with an all-acute callback | 1/24 | **2/24** |
| Final all-acute HA | 1/24 | **1/24** |

The extra transient all-acute event is retained, not hidden behind the unchanged
final count. The control does **not** prevent the known final inward HA.

## Known full-generated THR24 cases

Both H seeds independently select **330° (−30°)** from their initial energy grids.

| Measurement | Seed 1 unchanged → selected | Seed 2 unchanged → selected |
|---|---:|---:|
| Initial whole-System energy, kJ/mol | 145674.319 → 145561.857 | 148491.074 → 148229.376 |
| Initial HA nonbonded outward force, kJ/mol/nm | −56105.474 → −46482.938 | −142933.015 → −118877.690 |
| Initial HA–HG23 distance | 1.356899 → 1.377462 Å | 1.310981 → 1.327682 Å |
| First callback CA–HA radius | 2.506358 → 1.983431 Å | 5.373800 → 0.528361 Å |
| Final outward cosine | −0.994893 → **−0.995005** | +0.912712 → +0.923532 |

For seed 1 the initial energy and inward force decrease, but the HA still crosses
at the first callback and ends inward with angles **64.215/72.834/78.130°**
(N/CB/C). Seed 2's selected path acquires a transient all-acute event, although
its final HA is not all-acute. Lower initial global energy is not demonstrated to
be an acceptance predictor or a cure for this known failure.

## Why this does not test tetrahedral methyl initialization

The installed native initializer places multiple H near the same parent-directed
vector before its internal relaxation. Static inspection of these original
starts gives H–CG2–H angles **6.13/18.55/18.88°** in the seed-1 focal case and
**4.22/2.97/7.10°** in seed 2, not approximately tetrahedral angles. A common rigid
rotation deliberately preserves that clustered geometry. This experiment tests
its environment-dependent orientation, **not** replacement with a canonical
tetrahedral H arrangement. It cannot be reinterpreted as rejection or acceptance
of that separate construction.

The next separate prospective control will compare unchanged, fixed-phase
tetrahedral and initial-energy-selected tetrahedral methyl starts, with the same
original C–H radii and no other initial coordinate changes. It must freeze all
angles, phase/tie rules, budgets and matched cases before execution. Do not tune
the existing grid, pick successful seeds, silently modify the force field or
promote a neighboring-HA rewrite into production.

Validation: **17 focused tests** cover independent geometry, selection, identity
and refusal controls; the full fast suite passes **2,880 tests**, three optional
Torch/Biotite skips and 33 slow tests deselected. Targeted Ruff passes; independent
publications regenerate byte-identically. This is not broad scientific/platform
qualification or a resolution of original THR2. Empirical-table production and
license gates remain unchanged.

Code: `scripts/research_methyl_start_control.py`,
`scripts/audit_methyl_start_control.py`.
