# Reconstruction follow-up: natural gaps and the 1HE8 anchor failure

Date: 2026-11-12. User direction: small artificially hidden reference panels are
insufficient; prioritize real incomplete structures and make difficult cases such
as 1HE8 work. This records the next scientific scope, not completed new sampling.

## What success means

The task is usable gap completion, **not recovery of an unknown crystal pose**.
Primary endpoints are complete intended sequence/heavy atoms, correct junctions,
no severe clashes or D-Cα, exact nonmobile observations, bounded explicitly mobile
observations, compatible simultaneous gaps and retained topology warnings.
Reference RMSD is optional descriptive evidence for artificial-mask controls only:
it is not an acceptance gate, donor-selection objective or reason to reject a
geometrically acceptable reconstruction. Historical scoring/output is unchanged.

The small panels and independent 1PGA/2CI2 study establish transport and bounded
workflow evidence; they do not constitute a representative production benchmark.

## Existing real sources, not more artificial masks

Retained experimental downloads already include 9CQD, 1N8Z, 1JPS, 1AHW, 1G9M,
1GC1, 1BRS, 1I2M, 1KTZ and 1HE8. Inventory the original biological assemblies,
deposited sequence/observations, actual missing spans, partner context, partial
occupancy, missing observed heavy atoms, disulfides and heterogens before choosing
new backend cells. Earlier cc94 runs on this cohort are historical evidence,
**not** new RFD3/Boltz/Protenix execution.

The retained assembly-1 preparation inventory has substantial natural-gap tasks:

| Entry | Target residues, including gaps | Chains | Natural missing spans (residue counts) |
|---|---:|---:|---|
| 1HE8 | 1131 | 2 | 31,4,21,9,92,7,6,13,9,24 |
| 1AHW | 1294 | 6 | 3,8,8,3,8,8 |
| 1N8Z | 1041 | 3 | 9,3,4,10 |
| 1G9M | 949 | 4 | 4,12,4 |
| 1GC1 | 948 | 4 | 11,13,4 |
| 1JPS | 658 | 3 | 1,4,8,4,7,8 |
| 1I2M | 618 | 2 | 7,5,39,4,6,4 |
| 9CQD | 530 | 3 | 8,14,3,3,3,2,40 |
| 1KTZ | 468 | 4 | 12,18,10,6,12,18,10,6 |

These are inventory counts from the retained preparation summary, not fresh
admission or successful reconstruction claims. In particular, 1N8Z's initial
preparation refused severe remaining geometry defects; retain the separate
follow-up lineage rather than call that initial refusal a pass. Small 1BRS
(199 target residues, two two-residue gaps) can remain an interface/control case,
not the main evidence for large-structure readiness.

Keep original-source hashes and refusals. No trimming to a tiny protein merely to
make an inference shape fit, no occupancy rewriting and no missingness invented
by reference masking. An explicitly scoped protein-only projection remains a
projection, not completion of omitted chemistry. Compare the same deposited
missing-region task across backends; where whole-context inference cannot fit,
declare source-only local conditioning separately and assess the complete assembly.

Freeze the admitted cases, actual sequence placements, generated masks, source
preparation, backend/context arms and serial pressure-v2 call budget before the
next samples. Keep six development seeds, failures and native/restored/repair
stages; do not silently expand to a full factorial or reuse validation seeds for
tuning. Reuse existing donors for the first 1HE8 relaxation investigation.

## What the anchor experiments actually did

`scripts/refine_patchwork_all_gaps.py` constructs **one simultaneous optimization**:
all generated atoms in all gaps, plus selected observed anchors, are mobile.
Other observed heavy atoms are massless/fixed but retain their interactions in
the OpenMM force-field system. Hydrogen placement is temporary; output ownership
and source-record audits are independent. Anchors have positional restraints and
explicit movement checks relative to the original source, not ratcheted references.
Bonded/soft-steric preconditioning and full-force-field minimization are separate
recorded phases. This is local energy minimization, not global conformational
sampling or a demonstrated untangling algorithm.

In the retained full 1HE8 patchwork experiment:

- Ten gaps, **1712 generated atoms +148 atoms in19 observed anchor residues**.
- Initial bond/clash/amide findings were removed, but a Ramachandran finding
  remained (also represented by an aggregate generated-region finding).
- Three observed residues exceeded the historical movement gate:
  A321 backbone RMS0.5291/Cα0.5003 Å; A353 Cα0.5020 Å; A895 backbone RMS0.5447,
  Cα0.5312 and heavy maximum1.2612 Å. Historical limits were0.5/0.5/1.25 Å.
- Therefore the rejected result does **not** mean that anchor relaxation had no
  effect. It improved the proposal but did not pass all recorded gates.

The separate richer-context two-anchor control moved **all generated atoms plus
A895/A902**, not only those two residues. It retained an amide at A901–A902 with
ω≈−144.6° and failed its separately declared tighter anchor envelope. Its limits
must not be confused with the earlier full-patchwork envelope. No RMSD-to-crystal
gate caused either rejection.

Evidence paths under the retained `tmp/` tree:

- `1HE8-patchwork-all-generated-anchors-v1/attempt-01/`:
  `protocol.json`, `movements.json`, `summary.json`, all trial structures.
- `mac-reconstruction-2026-11-11/1he8-two-anchor-engineering-feasibility/`:
  `declaration.json`, `summary.json`, independent audits and final structure.

## Next diagnostic sequence

1. Re-audit exact movable atoms, peptide/disulfide topology and fixed neighboring
   atoms for each retained 1HE8 attempt; distinguish geometry and movement failures.
2. For A901–A902 and the remaining Ramachandran center, inspect force components,
   gradients, convergence and constraint/positional-wall activity. Energy reduction
   alone does not demonstrate removal of a diagnosed junction defect.
3. Compare direct anchors against a predeclared sequence-expanded window and
   spatial neighboring-residue shell while keeping **all generated regions mobile
   together**. Existing `spatial_relaxation.py` is implementation, not proof that
   this additional controlled 1HE8 matrix has executed successfully.
4. Retain initial/final and useful intermediate geometry/topology diagnostics.
   A clean endpoint cannot certify an unsaved trajectory. Do not solve a ring
   piercing by permitting atoms to pass through the protein during optimization.
5. If a starting donor needs a different conformation, compare another compatible
   donor or separately declared local reconstruction. More minimizer iterations
   and looser anchor limits are not assumed solutions.

Do not silently relax scientific gates to pass a favorite case. If a movement or
geometry policy is changed, declare it explicitly, justify it against the actual
completion objective and independently report old/new outcomes. New whole-1HE8
donor approval remains a separate manual decision; no public capability or P3 is
enabled by this plan.

## Publication and navigation

Code/tests/reports and checksum-bound research archives are published on a new
branch from `main`; local environments, weights, browser profiles and unrelated
force-field directories are not part of that publication. See the
[comparison report](../research/mac-reconstruction-comparison-2026-11-11.md),
[original handoff](mac-reconstruction-portfolios-handoff.md) and
[evidence archive instructions](../research/evidence/mac-reconstruction-2026-11-12/README.md).
