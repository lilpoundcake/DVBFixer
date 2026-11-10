# Tetrahedral methyl initialization avoids the retained inward H-stage endpoint

Status: **case-specific prospective H-stage control**, not complete H geometry,
whole-system acceptance, a production repair or original THR2 causality. Follows
the [negative rigid-rotation control](methyl-start-control.md). Production remains
unchanged.

## Frozen three-arm experiment

[Complete evidence](tetra-methyl-control-evidence.json) retains **72/72 completed
cells**: all 24 original H starts × unchanged / fixed-phase tetrahedral /
initial-energy-selected tetrahedral. There are no retries, timeouts, blockers or
budget nonexecutions. The cohort remains four masks, two original H seeds and
three source-geometry arms of **one protein**, not 24 independent proteins.

Only generated THR24 HG21/HG22/HG23 may change initially. Preserve each original
C–H radius; set the three directions to canonical sp3 angular geometry about the
fixed CG2-to-CB axis. A deterministic CA-projection frame defines azimuth; the H
labels receive phases 0/120/240°. The fixed-phase arm uses azimuth zero. The
selected arm uses the original whole-H-System initial energy on exactly 24 phases
(0–345°, 15° spacing), with the same 10⁻⁸ kJ/mol tie rule as the preceding control.
No final/native accuracy is scored. Undefined frames or observed H are refused.

All nonmethyl starting coordinates remain exact, including GLN16 HA; all heavy
coordinates remain exact throughout H minimization. No force, parameter, mass,
exclusion, seed, H regeneration or whole-system stage is changed. The original
Reference / tolerance 1.0 kJ/mol/nm / maximum 50-iteration H minimizer is reused.
Native callback reporting only copies data and returns `False`.

The plan is private at `tmp/pdbfixer-research/tetra-methyl-plan-01.json`, SHA256
`97d88724b5f984f87512d500357ff2b8c4096c62590c5d4b3198340e0cfd2102`.
Budgets were 900 seconds per cell and 5,400 seconds total. Raw coordinates,
candidate H grids and gradients remain in `tetra-methyl-run-01/` beside it;
those paths do not transfer through Git.

## Independent acceptance and full denominator

All **24 unchanged H endpoints reproduce exactly**. The independent auditor
checks original radii, canonical H–C–H/H–C–CB angles, right-handed azimuth, untouched
atoms, all grid energies/selection and original System serialization. Canonical
H–C–H angles are 109.471221°; the original clustered starts range from about
1.7° to 27.4°, and their complete per-case angle inventories are retained.

At **all 3,600 callbacks**, physical forces on every mobile H equal the negative
native objective gradient, maximum component error **0.0 kJ/mol/nm**. All 3,744
initial/callback/final nodes, including final H-bond deviations, remain published.

| GLN16 HA descriptor | Unchanged | Tetrahedral zero | Tetrahedral selected |
|---|---:|---:|---:|
| Any inward callback | 2/24 | 1/24 | **0/24** |
| Any all-acute callback | 1/24 | **0/24** | **0/24** |
| Final all-acute | 1/24 | **0/24** | **0/24** |

The zero-phase seed-1 focal trajectory still transiently crosses inward; do not
report both tetrahedral arms as having no crossings. The selected arm has none
at recorded callbacks, not a claim about every internal line-search evaluation.

## Known full-generated THR24 / seed-1 case

The energy-selected canonical phase is **30°**, chosen before minimization.

| Measurement | Unchanged | Tetrahedral zero | Tetrahedral selected |
|---|---:|---:|---:|
| Initial energy, kJ/mol | 145674.319 | 143782.351 | 143721.550 |
| Initial HA nonbonded outward force, kJ/mol/nm | −56105.474 | −19459.103 | −13734.044 |
| First callback CA–HA radius | 2.506358 Å | 0.835825 Å | **0.470526 Å** |
| Final outward cosine | −0.994893 | +0.919153 | +0.913471 |
| Final CA–HA radius | 1.076222 Å | 0.989758 Å | **0.971753 Å** |

Changing only the initial three-H angular construction, while keeping all original
forces and other coordinates, prevents the known final all-acute HA in this
H-stage case. The preceding rigid rotation did not. This supports a specific
initial-angular-geometry mechanism, not a general theorem about every methyl or
a demonstration that a lower initial energy alone predicts success.

## Important remaining H strain

Avoiding the inward/all-acute endpoint is **not complete H geometry acceptance**.
The selected focal endpoint has HA angles **103.584/89.102/129.787°** (N/CB/C),
not approximately 109.5°, and its CA–HA bond is about **0.118 Å shorter** than the
force-field reference. Its first reported step also compresses this bond severely.
Seed 2 likewise retains compressed CA–HA geometry. The remaining fixed heavy
environment still imposes strain; zero recorded inward callbacks is a much
narrower result than a qualified initial-H method or a well-converged minimum.

No final heavy Cα moves during this study, so its later flattening cannot be
declared prevented on this evidence alone. The existing post-H repair also
changes some H between this boundary and full-system minimization. The separately
prepared downstream bridge reproduces all **24 original whole-system starts
exactly in native nm coordinates**, using the unchanged production repair once;
downstream whole-system sampling/evaluation is a distinct control, not inferred
from this H-stage result.

Validation through this H-stage checkpoint: **16 new focused tests** and targeted
Ruff pass; full fast suite **2,896 passed**, three optional Torch/Biotite skips and
33 slow tests deselected. Independent publications regenerate byte-identically.
The bridge has separate focused/preflight acceptance. Original THR2 causality,
remaining target THR/VAL failures, broader chemistry/platform qualification and
empirical-table production/license approval remain open.

Code: `scripts/research_tetra_methyl_control.py`,
`scripts/audit_tetra_methyl_control.py`.
