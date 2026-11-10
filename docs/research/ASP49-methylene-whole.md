# Whole-phase CH₂ control: ASP49 rescue transfers a defect to ALA151 locally

Status: **one-protein native-phase mechanism control**, not production repair or
general qualification. Follows the [prospective CH₂ H-stage experiment](ASP49-methylene-control.md).
All earlier adverse results, source refusals and production safeguards remain.

## All 18 continuations, six exact matched controls

[Complete evidence](ASP49-methylene-whole-evidence.json) retains **18/18 completed
cells and all 54 initial/phase-1/phase-2 states**. These are six existing 1TIM
contexts (three methyl source arms × two seeds), each with unchanged / local
LEU48 CH₂ / all-carbon CH₂ H output. There are no new seeds, retries, timeouts,
execution blockers, budget nonexecutions or favorable-case exclusions.

Use the unchanged native H-label repair once, followed by the captured object-index
mapping into full labels. All six unchanged starts **and both endpoints reproduce
their corresponding prior source-arm whole trajectories exactly in nm**, all
atoms. Each phase's serialized System agrees with that matched baseline. Original
forces, anchors, State/integrator and the native tenfold phase-2 positional-weight
reduction remain unchanged: **Reference / no solvent / tolerance 10 kJ/mol/nm /
maximum 1,000 iterations per phase**.

These are **two native phases before production's emergency recovery**, not the
complete production pipeline. No emergency reflection, observation restoration,
candidate bundle or PDB is emitted. No production/default/threshold changed.
Budgets are 1,800 seconds per cell and 10,800 seconds total.

Plan SHA256:
`1d9871f16b8b53aade2f9415141675833707f3de857a9865bf81f88ae6e8a238`.
Private plan/run/audits are `ASP49-methylene-whole-*` under
`tmp/pdbfixer-research/`; those runtime paths do not transfer through Git.

## Aggregate counts hide the transferred defect

All raw and PDB-rounded signed Cα/THR-Cβ volumes, both focal sites, H descriptors
and original-heavy motion remain reported. The 0.5 Å³ shape descriptor is
unchanged and **not** a complete geometry/chemistry/steric acceptance gate.

| Final descriptor across six contexts | Unchanged | Local LEU48 CH₂ | All-carbon CH₂ |
|---|---:|---:|---:|
| Negative Cα | 0 | 0 | 0 |
| Cα at/below 0.5 Å³ | **1** | **1** | **0** |
| THR-Cβ at/below 0.5 Å³ | 0 | 0 | 0 |
| All-acute HA | 0 | 0 | 0 |

The unchanged flag is the previously retained ASP49/seed-1/selected-source
endpoint at **0.498283 Å³**; its unmodified-methyl-source counterpart is also
substantially flattened at **0.500677 Å³**, just above the threshold. Preserve
both continuous values rather than overstating their binary split.

The local arm rescues ASP49 in every context, but its single global flag moves to
**ALA151 / seed 2 / graph-tetrahedral-zero source**. Equal aggregate flag counts
therefore do not mean equal structures or absence of a new adverse outcome.

## ASP49 rescue is sustained downstream

For seed 1 / unchanged methyl source, all arms start at identical heavy volume:

| ASP49 Cα signed volume | Unchanged | Local CH₂ | All-carbon CH₂ |
|---|---:|---:|---:|
| Whole-system start | 2.575429 Å³ | 2.575429 Å³ | 2.575429 Å³ |
| Phase 1 | 1.065215 Å³ | 2.588732 Å³ | 2.586063 Å³ |
| Phase 2 | **0.500677 Å³** | **2.612701 Å³** | **2.610218 Å³** |

All local-arm final ASP49 volumes lie approximately **2.609–2.613 Å³**; all-carbon
volumes lie **2.608–2.612 Å³**. Both retain near-tetrahedral focal HA angles instead
of the native inward/distorted arrangement. This supports the specific controlled
initial-condition mechanism without directly rewriting ASP49 HA or its heavy frame.

## New ALA151 flattening survives both native phases

For seed 2 / graph-tetrahedral-zero source, original ALA151 heavy geometry is
identical across arms:

| ALA151 Cα signed volume | Unchanged | Local CH₂ | All-carbon CH₂ |
|---|---:|---:|---:|
| Whole-system start | 2.315224 Å³ | 2.315224 Å³ | 2.315224 Å³ |
| Phase 1 | 2.369564 Å³ | **0.460936 Å³** | 2.440499 Å³ |
| Phase 2 | **2.466134 Å³** | **0.445384 Å³** | **2.507008 Å³** |

The local arm's final ALA151 HA remains inward (outward cosine **−0.982127**), with
N/C/CB-to-HA angles **84.318/91.124/85.969°** and CA–HA **1.088678 Å**. It is no
longer formally all-acute because one angle exceeds 90°, yet the observed-heavy
Cα is strongly flattened. Zero final all-acute count and zero D-Cα therefore miss
this shape defect. The substantial change from 2.466 to 0.445 Å³ is not merely a
borderline threshold classification.

Only two LEU48 H were changed initially, approximately 27 Å away. Their direct
initial pair interaction with ALA151 HA is absent beyond the native H-System
cutoff, as documented in the H-stage note. This does not isolate indirect
environment coupling from global optimizer history and does not establish a
nonlocal physical force or unique pair causality.

## Narrow all-carbon positive result is not broad qualification

All-carbon CH₂ has no half-Å³ findings in these six contexts; its minimum final
Cα volume across all centres is **2.162403 Å³**. Its maximum final H-angle deviation
is **11.915°**, versus **28.629°** unchanged and **25.182°** local. Median per-cell
maximum H-bond deviations are **0.03211/0.03214/0.03235 Å**, respectively. These
are descriptive readbacks, not new favorable acceptance thresholds.

The broad arm changes **620 generated H** at 310 carbon-CH₂ parents, not two.
This is one protein with correlated source arms, no native-accuracy comparison,
not validation on glycans/ligands/PTMs or other force fields/platforms. Earlier
1SHG fixed-phase adverse outcomes remain unchanged. Do not promote a broad
initializer or call this complete-complex/convergence qualification.

Original heavy observations are not frozen during whole relaxation. For the new
ALA151-failure context, **1,870 original heavy atoms** move with RMS/max
**0.285218/1.881413 Å** unchanged, **0.287578/1.921888 Å** local and
**0.284049/1.910602 Å** all-carbon. No observations were restored afterward.

## Validation and next boundary

Five focused native-baseline/parameter-drift/focal-readback/saved-only IO tests
pass. Full fast suite through initial downstream tests: **3,004 passed**, three
optional Torch/Biotite skips and 33 slow deselected; the later end-to-end saved-only
audit test also passes separately. Ruff, agent-map and generated CLI/GUI checks
pass; independent evidence publications regenerate byte-identically.

Next: freeze a source-diverse carbon-multi-H factorial control using the **existing
complete native starts and known adverse reconstruction contexts**, not a new
success-selected protein/seed search. Keep carbon-CH₂ and methyl contributions
separate and include their fixed-zero combination, all unchanged endpoints and
source/geometry refusals. Label reused positive cases and intentionally poor-heavy
stress contexts explicitly; none establishes general repair qualification.
Original unforced THR2 Cβ causality and empirical-table production provenance/
license/approval remain unresolved. Production/public completion remains disabled.

Code: `scripts/research_methylene_whole_control.py`,
`scripts/audit_methylene_whole_control.py`.
