# ASP49 inward H paths: delayed methylene-H repulsion, not the original methyl impulse

Status: **retrospective saved-path force attribution**, not a controlled pair
intervention or production repair. Follows the
[source-diverse whole-phase result](graph-methyl-whole.md). Original unforced
1CRN THR2 Cβ causality remains separate and unresolved.

## Frozen complete matched readback

[Independent evidence](graph-methyl-ASP49-attribution-evidence.json) covers **all
three 1TIM arms and both existing seeds**: six contexts, **312** saved initial/
callback/final nodes and **1,172,184** HA–other-particle pair observations. Selection
of ASP49 is retrospective to its already-known H/whole outcomes, not blinded
discovery. No unsuccessful seed or matched positive context is omitted.

Fresh Reference Contexts only read forces at unchanged saved coordinates. No new
H, minimization, integration, force ablation, exclusion/parameter edit or native
accuracy scoring occurs. Measurement groups are changed only on a cloned System;
restoring group metadata must restore its exact original serialized force field.
Analytic pairs reuse the tested native Lorentz–Berthelot LJ/reaction-field Coulomb
formulas with exact native exceptions. All full three-dimensional pair sums close
to unmodified native NonbondedForce, maximum component error **1.048e−9
kJ/mol/nm**. Independent summed force groups also match the original native
callback gradient, maximum component error **5.821e−11 kJ/mol/nm**.

Every node retains force-group components, methyl-H/other-H/heavy and same/other
residue totals, and exactly the three most-negative outward-axis pairs. All pair
terms remain private with checksums. The outward axis is opposite the sum of the
three unit CA–N/CA–C/CA–CB bonds in the unchanged observed-heavy frame. A negative
outward cosine is an inward hemisphere descriptor, not the same as all three
HA angles being acute.

Plan SHA256:
`b7273ab667e61c394fb9717516d41ea7d3ca6b5f6adfc13c95e3f594e858c1e2`.
Private plan/run/audits are `graph-methyl-ASP49-attribution-*` below
`tmp/pdbfixer-research/`; these raw paths do not transfer through Git.

## Initial local force does not reproduce the 1SHG first-step mechanism

For seed 1, ASP49 HA initially points outward, cosine **+0.991295**, CA–HA **1.0 Å**
in every arm. Its initial nonbonded outward force is only
**+16.424 / +17.077 / +18.744 kJ/mol/nm** (unchanged/zero/selected), while bond
stretching contributes **+2538.318** and angles **+119.911** in all three. Heavy
and nonmethyl-H initial forces are identical across arms; the small initial
nonbonded differences are entirely the changed methyl-H arrangement.

Thus the much lower selected initial global energy does not produce a large
direct initial inward impulse on this HA. The prior
[1SHG first-step methyl contact](ha-nonbonded-pairs.md) must not be copied as the
explanation of this new protein's trajectory.

## Seed-1 crossing is delayed and dominated by LEU48 methylene H

| Saved trajectory | First inward callback | Inward callbacks | Final outward cosine |
|---|---|---:|---:|
| Seed 1 unchanged | 11 | 39 | −0.959424 |
| Seed 1 graph zero | none | 0 | +0.999687 |
| Seed 1 selected | 12 | 38 | −0.958368 |
| Seed 2 unchanged | none | 0 | +0.999677 |
| Seed 2 graph zero | none | 0 | +0.999685 |
| Seed 2 selected | **26** | **5** | +0.984010 |

Callbacks are zero-based and do not expose every internal line-search evaluation.
In seed-1 unchanged, immediately before the first recorded crossing:

| Node | HA outward cosine | CA–HA | LEU48 HB2 distance | HB2 outward force |
|---|---:|---:|---:|---:|
| Callback 9 | +0.998310 | 1.184603 Å | 1.285698 Å | −79598.175 |
| Callback 10 | +0.929693 | **0.768799 Å** | **1.240141 Å** | **−145750.713** |
| Callback 11 | −0.674820 | 1.511720 Å | — | — |

Force units are kJ/mol/nm. At callback 10 the HB2 projection is almost entirely
LJ (**−145715.849**), while all methyl-H partners together contribute only
**−29.343**. The leading pair after crossing changes to LEU48 C;
the complete evolving list is retained, rather than assuming one pair dominates
every later node.

The actual native topology identifies **LEU48 HB2 as a methylene hydrogen**:
parent CB is bonded to CA, CG, HB2 and HB3. It is **not in the graph-methyl
intervention**. In selected seed 1, callback 11 likewise has HA–HB2
**1.146865 Å**, inward pair projection **−416822.384** (LJ −416780.277), followed
by the first recorded inward state at callback 12. The graph-zero trajectory
avoids this crossing even though it changes only methyl H initially; that is
evidence of coupled optimizer/environment sensitivity, not direct initial HB2
rewriting or an isolated proof of pair causality.

## Seed-2 selected transient: a different path remains adverse

The selected seed-2 path crosses inward at callbacks 26–30 and returns outward
by its endpoint. Before this crossing, **ILE86 HD13**, a genuine graph-methyl H,
has distance **1.489215 Å** at callback 25 and inward projection **−5218.624**
(LJ −5214.777). At callback 26 the largest inward pair is instead **ASP49 CG**,
while total nonbonded outward force has already become positive. This transient
must not be erased by its outward final state or the earlier all-acute count.

A quasi-Newton/global line-search path is not determined by one scalar pair
projection at a single callback. These are exact force/path observations, not
weight-only causality, an ablation, or a convergence certificate at the fixed
50-iteration limit.

## Next controlled boundary

A separate prospective initial-condition control can test canonical geometry of
the **two generated LEU48 CB hydrogens**, leaving ASP49 HA, all heavy atoms, other
H and every force unchanged. Include all three existing source arms and both
seeds, matched unchanged endpoints and a separately declared all-carbon-methylene
arm. Freeze the heavy-plane formula, deterministic H-slot assignment, original
radii, refusal rules and fixed Reference/1.0/50 budget before execution. Do not
choose the slot or seed by final geometry, increase iterations, remove the pair
or silently broaden the production initializer.

Validation: six focused frame/native-pair/audit tests and targeted Ruff pass.
The complete fast suite through this attribution checkpoint passes
**2,978 tests**, three optional Torch/Biotite skips and 33 slow deselected.
Independent publications regenerate
byte-identically. No empirical tables, raw coordinates or force-field tables are
redistributed, and production/default/threshold/license gates are unchanged.

Code: `scripts/attribute_graph_methyl_ha.py`,
`scripts/audit_graph_methyl_attribution.py`.
