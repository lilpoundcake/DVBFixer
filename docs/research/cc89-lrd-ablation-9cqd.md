# 9CQD: paired local refinement diffusion ablation

Date: 2026-10-23. Research-only experiment; not a qualified backend or a
reproduction of PATCHR's published Boltz/Protenix benchmark.

## Question and frozen matrix

Can local neural re-denoising repair the supplied 9CQD gap junctions better
than the current generated-only OpenMM refinement? Compare:

1. Initial sampling without versus with per-step fixed-atom reinjection.
2. LRD noise/mobility on generated atoms only versus including nearby original
   observed residues.
3. OpenMM-only, LRD-only and LRD→OpenMM processing.

There are **10 distinct final candidates**, not 12: the two OpenMM-only controls
do not depend on an LRD mobility setting. Both initial seed-7 candidates were
regenerated locally in the same pinned environment. No candidate or seed was
selected using reference coordinates.

See the [experiment plan](../plans/cc89-local-refinement-ablation.md) and
[PATCHR comparison](patchr-gap-reconstruction-comparison.md). Joint generation
of nearby gaps is **deferred at user request** until a suitable test structure
is available; it is not part of this experiment.

## Inputs and implementation

- Original PDB SHA-256:
  `6bcc9638c3d299ea12151c3894b7e6c3482ebd2fdfe1d8e54dcc1f6f6316a281`.
- Original archived request SHA-256:
  `713d5fbc971b46bdc69cbf8a4215e216d8747b1cedf39b7d6f8e464c33e32246`.
- Execution request changes only the profile label to the private research
  identity; sequence, masks, contexts and target placement are unchanged.
- Targets M/W/X: 27/214/221 residues; missing W26–28 and X140–141.
- Protpardelle-1c cc89 checkpoint:
  `dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483`.
- Pinned PyTorch 2.6.0 MPS, float32, fallback disabled; existing Apple portability
  and per-step callback patches in an isolated engine copy.
- Initial generation: 500 steps per gap, seed 7, step scale 1.2, s_churn 200,
  source-only independent local contexts, physical sequence positions and
  neutral denoiser chain phase. No sequence design/MPNN.

Research code: [runner](../../scripts/run_cc89_lrd_ablation.py),
[pure protocol helpers](../../scripts/cc89_lrd_protocol.py),
[unit tests](../../tests/test_cc89_lrd_protocol.py).

The LRD stage calls the pretrained cc89 neural denoiser, not OpenMM and not a
new full-noise structure generation. It starts from the frozen raw candidate
with masked Gaussian noise at sigma 1.5 Å, then performs 25 Euler denoising
updates using a rho-7 schedule, positive minimum sigma 0.0064 Å, terminal sigma
zero, step scale 1.2 and no additional churn. Inactive atoms are hard-projected
to their reference coordinates after every update. This reset is intrinsic to
locality in both LRD arms; the **reinjection experimental factor refers to the
initial 500-step generation**, not another LRD reset toggle.

Noise is seeded by 7 + gap index and shared at common atoms across arms.
All four LRD calls use the same original-observation centroid for a given gap,
independent of mobility and initial-sampling arm. No global transform or final
original-coordinate overwrite is used to hide observed movement.

The window uses target-sequence positions, not residue-number arithmetic:

| LRD policy | W window | X window | Total mobile heavy atoms |
|---|---|---|---:|
| Generated-only | W26–28 | X140–141 | 30 |
| Observed-mobile | W23–31 | X137–144 | 100 (30 generated + 70 observed) |

Absent source sidechains are not invented as context. Other generated gaps
are not used as conditioning. Both gap views are denoised independently and
only their disjoint mobile updates are merged.

OpenMM-only and LRD→OpenMM use the unchanged v8 scientific refiner on CPU,
one thread, no solvent, 500 iterations, eight restarts, 0.25 Å perturbation,
seed 7. It moves generated atoms only. In an observed-mobile combination,
observed atoms remain frozen at their LRD positions during OpenMM: allowing
OpenMM to additionally move those atoms would be another experimental factor,
not an assumption made silently here.

## Independent checks and evidence

Every candidate is checked against the **original** source and request with
the normal validator. Source-displacement failures are not removed from strict
acceptance. Geometry excluding those failures is reported separately and is
never called production acceptance.

For observed-mobile arms, diagnostics also cover the entire released window
and one residue beyond each edge, including new outer junctions and clashes.
This avoids relying only on the narrower normal generated-region neighborhood.
Inactive PDB atom records must remain byte-identical, including identity,
coordinates, occupancy and B factors. D-Cα and severe overlap checks remain on.

Final evidence directory:
`test/errors_diffusion/9cqd-lrd-factorial-seed7-fixed-frame/` (ignored).
It contains original inputs, both initial candidates and per-gap traces,
LRD traces/masks, code snapshots and hashes, all processed candidates and
`comparison.json`. No candidate is auto-published by this script.

An earlier exploratory run under
`test/errors_diffusion/9cqd-lrd-factorial-seed7/` used the inactive-mask centroid,
which changed the input frame between mobility arms. It is preserved as a
confounded pilot, **not the final comparison**. Its two initial candidates were
copied unchanged into the corrected-frame run, and only LRD and postprocessing
were repeated. No sigma/window/seed was tuned against the results.

Sampling evidence: each initial gap has zero callbacks in the off arm and
500 callbacks in the on arm; on-arm post-projection error is zero. Every LRD
gap has 25 denoiser updates, 50 first-layer RoPE calls, and zero inactive-state
projection error. Shared-frame equality was verified across all four LRD arms.

## Results

**0/10 candidates pass strict validation. None passes geometry even if source
displacement is reported separately rather than counted as a geometry failure.**
All ten candidates have zero detectable D-Cα. No gate was relaxed.

Counts below are the normal independent validator's generated/junction
non-planar amides and severe overlaps; expanded-window diagnostics additionally
check moved observed residues and their outer boundaries.

| Initial reinjection | Processing | LRD mobility | Non-planar amides | Severe overlaps | Max observed displacement (Å) | Strict pass |
|---|---|---|---:|---:|---:|---|
| Off | OpenMM-only | — | 1 | 0 | 0 | No |
| Off | LRD-only | Generated-only | 1 | 4 | 0 | No |
| Off | LRD→OpenMM | Generated-only | 2 | 0 | 0 | No |
| Off | LRD-only | Observed-mobile | 0 | 2 | 2.9385 | No |
| Off | LRD→OpenMM | Observed-mobile | 0 | 1 | 2.9385 | No |
| On | OpenMM-only | — | 2 | 0 | 0 | No |
| On | LRD-only | Generated-only | 1 | 4 | 0 | No |
| On | LRD→OpenMM | Generated-only | 2 | 0 | 0 | No |
| On | LRD-only | Observed-mobile | 0 | 2 | 2.9330 | No |
| On | LRD→OpenMM | Observed-mobile | 0 | 1 | 2.9330 | No |

### Without moving original atoms

OpenMM-only reproduces the earlier failures: W28→W29 omega −132.2902° without
reinjection; W25→W26 146.2973° and W28→W29 −116.0406° with reinjection.

Generated-only LRD does not resolve W28→W29 (about −129.8° / −129.7°), leaves a
junction C–N outside 1.20–1.45 Å and introduces four severe overlaps.
Following it with OpenMM removes those overlaps and fixes C–N connectivity,
but leaves two non-planar amides and a gross Ramachandran outlier at W25:

- Reinjection off: W27→W28 omega 144.359°, W28→W29 −142.880°.
- Reinjection on: W27→W28 omega 147.960°, W28→W29 −142.969°.

Observed-heavy coordinates remain exact in these arms.

### Releasing original boundary residues

Observed-mobile LRD removes the non-planar-amide finding in the generated
neighborhood, but moves 70 original heavy atoms (maximum approximately 2.94 Å),
introduces a gross Ramachandran outlier and leaves bad C–N connectivity plus
two severe overlaps. It is not a successful geometry repair hidden only behind
the fixed-coordinate gate.

OpenMM following this LRD fixes the generated-region C–N distances and removes
one clash. However, W25 and W29 remain Ramachandran outliers, and the
W24/GLY:O–W25/SER:O overlap remains approximately 1.01 Å (O–O distance 2.03 Å).
Both atoms are original observations released by LRD and frozen at those new
positions by the generated-only OpenMM stage. The result still fails geometry
even independently of its source displacement. A subsequent minimizer that
also releases observed atoms is a different, untested policy.

### Interpretation and limits

For this one structure/seed/schedule, allowing source-side motion trades an
amide-planarity problem for source drift and other geometry defects.
Generated-only LRD and extra initial reinjection do not provide an accepted
reconstruction. This does not establish that PATCHR itself fails this example:
the experiment uses a cc89 low-noise loop with declared settings, not the
PATCHR Boltz/Protenix LRD implementation. No hyperparameter sweep or reference
accuracy claim was made.

The full fast Python lane passes **1447 tests**, with two optional-dependency
skips. Nine new protocol tests cover boundary identity/positions, noise schedule,
static projection and mobility-independent source centering. Ruff and agent-doc
checks pass. Production samplers and validation thresholds are unchanged;
the runner remains an explicitly research-only script.

## Extension: minimize generated atoms and immediate anchor residues

At the user's request, six additional OpenMM-containing variants release the
four immediate anchors **W25/W29 and X139/X142**, in addition to W26–28 and
X140–141. The enlarged minimization mask contains 54 existing heavy atoms:
30 generated and 24 observed anchor atoms. It does not release the whole ±3
LRD window. No atom identity is reclassified as generated in the original
request or in independent validation.

The exact same frozen raw/LRD candidates are reused; only the minimization
mask changes. Seed, force field, CPU/thread policy, 500 iterations, eight
restarts and 0.25 Å perturbation remain unchanged. PDB records outside the
explicit minimization mask are checked byte-for-byte against that minimizer's
input. For an observed-mobile LRD input, this preserves earlier LRD movement
outside the immediate anchors; it does not hide it by restoring source atoms.

Results: **6/6 pass the tested geometry checks**, including expanded-window
and outer-junction diagnostics. All have zero non-planar amides, severe overlaps
and detectable D-Cα; the prior gross Ramachandran outliers are also absent.
**0/6 pass strict acceptance**, solely because original-coordinate drift exceeds
the unchanged preservation thresholds.

| Initial reinjection | Input to generated+anchor OpenMM | Geometry pass | Max anchor movement (Å) | Max other observed movement (Å) | Strict pass |
|---|---|---|---:|---:|---|
| Off | Raw (no LRD) | Yes | 1.6514 | 0 | No: source motion |
| Off | Generated-only LRD | Yes | 1.6849 | 0 | No: source motion |
| Off | Observed-mobile LRD | Yes | 3.5715 | 2.9120 | No: source motion |
| On | Raw (no LRD) | Yes | 1.9311 | 0 | No: source motion |
| On | Generated-only LRD | Yes | 1.7141 | 0 | No: source motion |
| On | Observed-mobile LRD | Yes | 2.9425 | 2.9124 | No: source motion |

For the raw/no-reinjection input, maximum anchor-backbone movement is 1.2614 Å
at W29:N; the maximum heavy-atom movement is W29:CG2 (1.6514 Å). Thus these
are real anchor adjustments, not merely hydrogen optimization or numeric drift.

For this structure/seed, permitting the immediate anchor residues to relax is
enough to remove the observed geometry defects without neural LRD or initial
reinjection. This is not evidence that moving experimental anchors is always
appropriate or that the resulting loop is experimentally accurate. A future
production mode would need an explicit user-approved mobility policy and its
own acceptance contract; the existing strict contract is not weakened here.

New evidence: `processed-with-anchors/`, `comparison-with-anchors.json`, frozen
script/protocol copies, settings and input bindings in
`anchor-minimization-manifest.json`, under the same fixed-frame evidence directory.
All six PDBs were additionally copied unchanged to
`test/errors_diffusion/9cqd-lrd-visual-review/` with `__anchors-mobile.pdb`
suffixes. That folder now contains 18 structures: two raw, ten original
processed outputs and six anchor-mobile additions. Previous files are intact.

Execution:

```bash
OPENMM_CPU_THREADS=1 python scripts/run_cc89_lrd_ablation.py evaluate \
  --minimize-anchors \
  --workspace test/errors_diffusion/9cqd-lrd-factorial-seed7-fixed-frame
```

The final fast lane now passes **1449 tests**, with two optional-dependency
skips. Eleven protocol tests include explicit anchor-mask and insertion-code
coverage. Ruff and agent documentation checks pass.

## Anchor backbone RMSD and MODELLER comparator

RMSD uses original source coordinates, **N/CA/C/O without any rigid alignment**.
Each residue contributes four atoms; the total is over all sixteen atoms, not
an arithmetic mean of four residue RMSDs. Units: Å.

| Reconstruction / LRD input | Initial reinjection | W25 | W29 | X139 | X142 | All anchors |
|---|---|---:|---:|---:|---:|---:|
| Diffusion raw → generated+anchor OpenMM | Off | 0.451 | 0.888 | 0.243 | 0.150 | 0.518 |
| Diffusion generated-only LRD → generated+anchor OpenMM | Off | 0.595 | 0.865 | 0.265 | 0.168 | 0.548 |
| Diffusion observed-mobile LRD → generated+anchor OpenMM | Off | 1.513 | 2.173 | 0.783 | 0.954 | 1.461 |
| Diffusion raw → generated+anchor OpenMM | On | 0.410 | 0.752 | 0.258 | 0.188 | 0.457 |
| Diffusion generated-only LRD → generated+anchor OpenMM | On | 0.413 | 0.952 | 0.237 | 0.180 | 0.540 |
| Diffusion observed-mobile LRD → generated+anchor OpenMM | On | 1.449 | 2.086 | 0.801 | 0.936 | 1.411 |
| MODELLER, no subsequent OpenMM | Not applicable | 2.300 | 2.523 | 5.975 | 5.496 | 4.404 |
| MODELLER → generated+anchor OpenMM | Not applicable | 2.603 | 2.565 | 5.851 | 5.261 | 4.338 |

The native MODELLER comparison uses version 10.8, one initial model and two
fast loop refinements, default `--pin-input`, default MODELLER seed −8123.
It is not described as a seed-7 MODELLER run. Selection is the lower native
molpdf loop candidate (47.73271), not an original-coordinate/RMSD selector.
The additional OpenMM stage uses the same seed-7 generated+anchor protocol as
the diffusion outputs. Target sequences M/W/X 27/214/221 and observed-to-target
placements are locked to the diffusion request through the existing
`scripts/run_locked_modeller_comparator.py`; no raw full-length FASTA tails are
reintroduced.

MODELLER's initial automodel optimization moves the observed structure globally;
`--pin-input` restricts only subsequent loop refinement, not that initial stage.
Since this table does not superpose coordinates, global/frame or chain movement
also contributes to MODELLER's numbers. They measure literal source-coordinate
preservation, not experimental loop accuracy or isolated local deformation.
The native MODELLER path can additionally complete source sidechains and infer
CONECT records: it is not the same generated-only identity/link contract as
diffusion. Applying that strict contract flags extra outside-generated atoms,
source drift and additional links; these are recorded rather than suppressed.

MODELLER evidence is retained under
`test/errors_diffusion/9cqd-modeller-comparator/`: locked targets, source,
work/model log, settings, both outputs, native `.dat`, independent audit and
full-precision `anchor-backbone-rmsd.{json,tsv}`. Metadata is reconciled only in
the independent audit view, never by silently altering the saved MODELLER PDBs.
Both PDBs and the eight-row RMSD table were added to the visual-review folder,
which now contains 20 PDBs.

### Updated metrics after W-chain Kabsch alignment

At the user's request, each of these eight candidates is rigidly aligned to
the original source using **844 common observed W-chain backbone atoms
(N/CA/C/O)**. The observed W25/W29 anchors are included in that chain-wide fit;
generated residues have no source counterpart and are excluded. The same proper
rotation and translation is applied to the entire complex. Chain X is **not**
aligned separately, so its relative position with respect to W remains part of
the comparison. Anchor RMSD is computed in double precision before PDB rounding.

| Reconstruction / LRD input | Initial reinjection | W25 | W29 | X139 | X142 | All anchors |
|---|---|---:|---:|---:|---:|---:|
| Diffusion raw → generated+anchor OpenMM | Off | 0.451 | 0.876 | 0.239 | 0.147 | 0.512 |
| Diffusion generated-only LRD → generated+anchor OpenMM | Off | 0.594 | 0.860 | 0.262 | 0.166 | 0.545 |
| Diffusion observed-mobile LRD → generated+anchor OpenMM | Off | 1.506 | 2.149 | 0.763 | 0.946 | 1.446 |
| Diffusion raw → generated+anchor OpenMM | On | 0.411 | 0.748 | 0.257 | 0.186 | 0.455 |
| Diffusion generated-only LRD → generated+anchor OpenMM | On | 0.416 | 0.944 | 0.234 | 0.177 | 0.536 |
| Diffusion observed-mobile LRD → generated+anchor OpenMM | On | 1.440 | 2.061 | 0.779 | 0.928 | 1.396 |
| MODELLER, no subsequent OpenMM | Not applicable | 0.696 | 1.093 | 0.566 | 0.421 | 0.738 |
| MODELLER → generated+anchor OpenMM | Not applicable | 1.422 | 1.649 | 0.423 | 0.582 | 1.147 |

Units: Å. The total remains RMSD over sixteen anchor-backbone atoms, not the
mean of four per-residue RMSDs. The W-chain fit RMSD is 0.323 Å for MODELLER
and 0.345 Å for MODELLER→OpenMM. Much of MODELLER's earlier unaligned 4.404 Å
anchor RMSD is therefore rigid-frame motion; it should not be interpreted as
local anchor deformation. Alignment does not change internal geometry or
retroactively satisfy the production exact-source-coordinate contract.

Eight aligned full-complex PDBs, rigid transforms and full-precision
`anchor-backbone-rmsd-kabsch-W.{json,tsv}` are saved in
`test/errors_diffusion/9cqd-lrd-visual-review/kabsch-W/`. The original twenty
top-level PDBs and the earlier unaligned table remain unchanged.

### Clarification: MODELLER without subsequent OpenMM

Clarified 2026-10-27: the **MODELLER, no subsequent OpenMM** row is the normal
loop-reconstruction output, not the MODELLER→OpenMM variant. MODELLER's own
initial optimization and fast loop refinement were performed; "without OpenMM"
does **not** mean an entirely unoptimized structure.

After the W-chain Kabsch alignment defined above, this variant has anchor
backbone RMSD (Å): W25 **0.696**, W29 **1.093**, X139 **0.566**, X142 **0.421**;
combined over all sixteen N/CA/C/O atoms: **0.738 Å**. These are the existing
measurements, not a new calculation or model run.

Corresponding inspection files:

- Before alignment:
  `test/errors_diffusion/9cqd-lrd-visual-review/modeller__raw.pdb`.
- After W-chain alignment:
  `test/errors_diffusion/9cqd-lrd-visual-review/kabsch-W/modeller__raw.pdb`.

The separate `modeller__openmm__anchors-mobile.pdb` files include the additional
generated+anchor OpenMM minimization and correspond to the **1.147 Å** aligned
combined RMSD row. Do not conflate these outputs.
