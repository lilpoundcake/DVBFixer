# Source geometry at fixed χ1: both THR rebuild arms pass this bounded final audit

Status: **private source-preparation intervention**, not a production fix, broad
qualification, native-accuracy improvement or complete explanation of the original
THR2 inversion. Follows the [finite Cβ wall experiment](beta-wall.md).

## Why this comparison differs from another wall-strength search

The wall stress retained correct Cβ handedness but left four angle failures. This
experiment instead changes the **geometry of atoms originally generated during
preparation**, without selecting another rotamer or adding a new stereochemical
force. The scientific question is whether upstream geometry preparation changes
the subsequent pipeline outcome under a controlled numerical procedure.

Before execution, `research_start_geometry.py` froze all six original legacy
prepared structures and `.dat` files, seeds **1 and 2**, and three arms:

| Arm | THR source change |
|---|---|
| `raw` | No heavy-coordinate change |
| `beta_only` | Reconstruct OG1/CG2 with canonical illustrative bond lengths/tetrahedral angles at the **original CB** |
| `full_generated` | Reconstruct CB from observed N/CA/C, then reconstruct OG1/CG2 |

Both rebuild arms retain the original measured **N–CA–CB–OG1 χ1**, before PDB
rounding. They do not consult native coordinates, enumerate alternative angles,
rank candidates, swap chemically distinct branches or reflect observed atoms.
The geometry functions are the existing illustrative research constructions, not
a newly calibrated residue/backbone-dependent geometry library.

All three target heavy atoms must be recorded as **added atoms in the original
`.dat` provenance**, even for the distal-only arm. Thus already wrong THR6 is
eligible here only because its affected atoms were generated, not authoritative
observations. This is a different task from preserving an originally observed
wrong-handed sidechain. VAL rebuild arms are explicit identity controls, not an
attempt to assign a stereogenic Cβ to VAL.

## Matched procedure and immutable source observations

Every cell runs the unchanged production CLI with explicit research controls:

```bash
dvbfixer minimize prepared.pdb --dat prepared.dat --no-solvent \
  --atom-naming standard --rebuild-h --platform Reference -o minimized.pdb
```

The isolated wrapper also selects **Reference inside `Modeller.addHydrogens`**;
the CLI platform flag alone does not control that internal Context. Python's
random seed is set only around H reconstruction and then restored. Every arm,
including `raw`, forces the same strip-and-readd-H procedure. This intentionally
differs from earlier default keep-H/OpenCL workflows and is not their exact replay.

Only declared generated-heavy coordinate columns are rewritten on private input
copies. Every other record, including observed heavy coordinates, names, serials,
occupancy, insertion codes and metadata, remains byte-identical. Old H records
are unchanged in those private intermediate copies, then explicitly stripped and
rebuilt **before** the physical minimization. No user input or `.dat` is altered.

Unrounded χ1 preservation is independently checked within 1e−8 degrees. The
largest measured χ change after ordinary PDB rounding is **0.039884°**, retained
as a rounding measurement, not silently reset or used to choose a new angle.

All planned cells run even if their initial geometry audit fails. That is an
explicit **mechanism experiment**, not a production policy admitting bad inputs
or publishing their partial products.

## Results: 36/36 workflows, 72 native phases

[Independent evidence](start-geometry-evidence.json):

| Arm | Completed | Final heavy-geometry passes | THR CB R | THR final geometry passes |
|---|---:|---:|---:|---:|
| Raw | 12/12 | **10/12** | 6/8 | 6/8 |
| Distal-only (`beta_only`) | 12/12 | **12/12** | 8/8 | 8/8 |
| Full generated sidechain | 12/12 | **12/12** | 8/8 | 8/8 |

There are zero execution blockers and no reported emergency Cα-repair messages.
All eight VAL identity contrasts have **exactly matching final heavy coordinates**
between `raw` and the corresponding no-op arm. No H-geometry qualification follows
from the heavy-atom identity controls.

The two raw failures are **THR6 under both seeds**, retaining its already wrong
CB S. Both geometry-reconstruction arms finish with CB R and accepted target
geometry for THR6, without changing source observations or selecting a different
χ1. This supports an upstream generated-geometry intervention on these inputs,
rather than merely increasing a soft signed-volume penalty.

**The original unforced THR2 R→S event still does not occur in this raw control.**
Its ordinary Reference/H-controlled runs already pass. Therefore the experiment
does not identify a unique cause of that original event or show that the proposed
rebuild prevents its exact trajectory.

### Initial failures are not erased by final success

Distal-only reconstruction leaves the original CB and its poor Cα geometry in
place. Its initial audit still flags near-planar Cα/angle findings and, in some
cases, target clashes. Full reconstruction improves these initial shapes, but
THR1 still has an initial OG1–A/37 O clash at the retained source χ1. All such
findings remain in the evidence, even though the corresponding final results pass.

This is not a qualification of either input construction as a complete clash-free
generator. Nor does final success imply that every intermediate optimization
state had acceptable geometry.

## Causal scope: same procedure is not identical H or identical potential

The intervention precedes H reconstruction. Changed heavy geometry can change
the subsequently placed H coordinates, even with the same algorithm/seed/platform.
Those H differences are **mediators**, not coordinates held identical across arms.

Distal-only reconstruction leaves the input CB fixed. Full reconstruction also
changes the generated-CB positional-restraint anchor and the input-derived Cα
guard reference volume. The original settings/force-field policy remain the same,
but this is **not a fixed-potential isolation of just one Cartesian coordinate or
one Cβ angle**. All actual systems and phase states are retained for further
controlled disaggregation. No undocumented causal claim is inferred from 12/12.

Final minimization legitimately moves original observations under the original
restraint policy. Their displacements are measured against the authoritative
**masked original input**, excluding the generated target atoms, and reported
separately from all-prepared-heavy motion. No original coordinate is restored
after it has moved; source-rewrite preservation is not final-output preservation.

## Independent audit and remaining gates

`audit_start_geometry.py` verifies complete source/arm/seed inventories, input
record changes, `.dat` ownership, retained χ1/rounding evidence, H-platform/seed
restoration, phase snapshots, final atom identity and raw final PDB bytes. It reads
no native reference coordinates and invokes no minimization or library queries.

The validator uses actual-frame heavy coordinates and captured external bonds.
Existing measured-χ target bond/angle/clash checks, global Cα/THR-Cβ volume checks
and RDKit target CIP are kept separate from source construction. Raw and rounded
phase coordinates and the actual final PDB are audited. This does not qualify H
geometry, every observed-residue angle, unrelated chemistry, or native accuracy.

The sample remains six deliberately selected masks from **two proteins**, with
two H seeds rather than twelve independent biological examples. ILE/PTMs, partial
observed sidechains, other source classes and new-input qualification remain open.
No arm is promoted to a production default or chosen using held-out RMSD.

Private artifacts under `tmp/pdbfixer-research/`: `start-geometry-plan-01.json`,
`start-geometry-run-01/`, original `start-geometry-audit-01.json` and independent
readback. Original inputs/logs/failures and all initial geometry findings remain.

```bash
PYTHONPATH=src:. python scripts/audit_start_geometry.py \
  --plan tmp/pdbfixer-research/start-geometry-plan-01.json \
  --run tmp/pdbfixer-research/start-geometry-run-01 \
  --output <fresh-evidence.json>
```

**820 focused tests passed in 28.41 s**; targeted Ruff and the 17-file completion-
package mypy check pass. Saved-only evidence regeneration is byte-identical.
This is not the full repository/platform suite. No production/default/dependency/
library/public capability changed; no commit or push.

Subsequent [new-source qualification and saved-force disaggregation](prospective-geometry.md)
show that the old distal-only force lists actually match their controls, but the
old12/12 geometry result does not generalize. New admitted inputs give2/8,3/8,5/8
gate passes and no native-accuracy superiority; all H-scope refusals remain visible.
