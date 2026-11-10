# Nine-χ close-contact controls: diversity increased, selected outputs unchanged

Status: **separate private sampling study**, not a production geometry-policy
change, empirical rotamer library, public completion capability or packer
qualification. Continues the [backbone-selected close-contact controls](rosetta-contact-controls.md).

**Main result:** seven matrices now have two varying source-feasible axes and show
full-score nonadditivity. Nevertheless, dense scoring selects **exactly the same
12 conformations as the coarse three-χ subset**. No accuracy or coverage gain over
that subset is established. Conditional factorized disagreements do not change
this conclusion and require the reference-pose caveat below.

## Enrollment before new sampling

`scripts/research_rosetta_dense.py enroll` froze every case from the previous
contact plan: all four pairs and all three mask classes, with the zero eligible
1CRN-pair outcome inherited unchanged. No reconstruction-success filtering,
replacement masks, new structure download or new sidechain-truth selection.

For paired missing branches, the explicitly ordered χ grid is:

```text
-60, 60, 180, -80, -40, 40, 80, 160, -160 degrees
```

It contains the original three candidates first, followed by ±20° offsets
(180°+20° wraps to −160°). At most nine proposals per paired target and **81 joint
combinations per case**. This illustrative grid is not an empirical probability
model or a calibrated conformational prior. The anchored mask still has exactly
one observed-partner algebraic closure: rotating it would move fixed observations.

The same `l_cb` and paired-branch internal-coordinate construction are used, with
no minimization or coordinate reflection. Only explicitly missing atoms are built.
The first three paired proposals/audits and the anchored proposal must reproduce
the frozen coarse results exactly before any new score is computed.

### Grid admission is separate from geometry and production

The original generator and its three-staggered audit remain unchanged and still
reject off-grid proposals. This new research adapter checks exact membership in
the nine-χ grid, then calls the existing **measured-χ comparator geometry**. Raw and
rounded bond/angle, steric, named-branch and chirality thresholds are unchanged.
Anchored validation remains the original anchored audit. Arbitrary torsion
metadata, unsupported empirical authority, observation movement and stereo
inversion are not accepted.

Private materialization independently rechecks source/joint geometry, exact source
records, rounded additions, global Cα/THR-Cβ and target RDKit CIP. Its standard
`DatRecord` sidecar owns only generated heavy atoms; the description explicitly
labels the new research regime. It is **not** the production or old-grid
materializer. Existing public readers/capabilities are not changed.

The unchanged score runner's compatibility field `old_source_accepted` here carries
the **new private grid's raw-joint eligibility**; it does not assert that these
off-grid candidates pass the old three-staggered generator policy. The enrollment
and score-plan schemas identify the separate nine-χ experiment.

## Exact scoring and independent readback

The unchanged `research_rosetta_joint.py score` still uses the same full-pose
`ref2015`, Rosetta revision `1a56185c2592611dec4c9c75ddc9468cd2227c1f`, seed 7,
fresh-H imports, two independently checked target views, fixed-context signatures
and complete-feasible-pool selection. Rejected raw combinations are scored only
as retained diagnostics; they cannot win or receive a materialized sidecar.
No packing, H optimization, relaxation or truth-based selection occurs.

`scripts/audit_rosetta_dense.py` replays saved proposals, source bindings, individual
audits/winners and the complete Cartesian search **without coordinate generation**.
It checks raw/rounded materialization, stereo, both physical/topology imports,
typed virtual separation, bounded terminal OXT, score-term sums/weights and every
selected private PDB/`.dat`. All **38 original coarse/native poses** reproduce
identical PDB bytes, full imported snapshots and exact scores across experiments.

[Portable independent evidence](rosetta-dense-controls-evidence.json) retains every
matrix and truth-free diagnostic, bound to the frozen inputs/implementations and
raw score evidence. Local snapshots remain required for full independent replay.

| Measure | Count |
|---|---:|
| Frozen / admitted cases | 12 / 12 |
| Complete source-feasible combinations | 141 |
| Raw joint rejections / additional rounded rejections | 12 / 0 |
| Private eligible outputs / generated-only sidecars | 129 / 129 |
| Native control poses, never eligible | 12 |
| Scored full poses / exact target-view imports | 153 / 306 |
| Selected private PDB/`.dat` pairs | 12 |
| Two-varying-axis / single-option-axis matrices | 7 / 5 |
| Changed coarse-score winners | **0** |

All 153 scored poses pass exact physical-heavy/topology and stereo readback.
Original observations remain unchanged. No rejected combination has a sidecar.
Native measured-geometry **9/12** and the three correlated THR89 source-contact
refusals remain the earlier separate control, not reclassified by dense sampling.

## What nonadditivity now measures

The enrolled diagnostic is unchanged:

```text
Delta(i,j) = E(i,j) - E(i,j0) - E(i0,j) + E(i0,j0)
```

The reference is the first source-feasible index pair, not native or best RMSD.
The complete source-feasible Cartesian score table, including diagnostic joint
rejections, is required. Full-score Rosetta units are not kcal/mol or an isolated
atom-pair interaction. Target H can differ with target heavy geometry; context
H/type/topology remains fixed under the original scoring contract.

| Pair / mask | Source-feasible dimensions | Max absolute Δ | Reference jointly eligible? |
|---|---:|---:|---|
| THR89/VAL92 distal | 3×5 | 0.42017 | Yes |
| VAL2/THR40 distal | 9×3 | 0.00552 | Yes |
| VAL2/THR40 whole | 5×3 | 0.00786 | Yes |
| THR43/THR51 distal | 8×4 | 48.69339 | **No** |
| THR43/THR51 whole | 7×2 | 5.44118 | **No** |
| THR51/THR69 distal | 4×3 | 0.06571 | Yes |
| THR51/THR69 whole | 2×4 | 0.02776 | Yes |

The other five matrices have one-option axes and remain uninformative, not
evidence of no physical interactions. These are correlated masks on four pairs
in one previously studied protein, not seven independent system validations.

## Conditional factorized baseline and its limitation

Before execution, `scripts/rosetta_dense_diagnostics.py` fixed the surrogate:

```text
E_factorized(i,j) = E(i,j0) + E(i0,j) - E(i0,j0)
```

It reports both independent axis minima and a surrogate minimum restricted to
the **same complete joint-geometry-feasible set** as full scoring. This is a
conditional reference-axis surrogate, not an independently fitted interaction
graph or a physical decomposition into independent residue energies.

Full score and the feasible factorized surrogate disagree only on THR43/THR51:

| Mask | Full-score reduction versus surrogate choice | Surrogate / full-score pooled RMSD |
|---|---:|---:|
| Distal/distal | 7.103 Rosetta units | 1.484 / 0.140 Å |
| Whole/whole | 4.632 Rosetta units | 0.315 / 0.209 Å |

**Both use a geometry-rejected reference pose.** The large reference-grid
nonadditivity and ranking differences must not be presented as a robust practical
packing improvement. No alternate reference was chosen after seeing the results.
Independent conditional minima and independent clearance winners both remain
jointly feasible in all twelve cases: there is no demonstrated winner rescue.

Most importantly, every full-score winner is already present and selected in the
coarse grid. Matched coarse and dense median pooled exact-name RMSD are both
**0.202370 Å**; there are zero changed output conformations, not merely equal
medians. The density expansion exposes interactions and rejected competitors,
but supplies **no added selected accuracy or coverage** here. Native accuracy
remains reporting-only; no template arm was executed and no superiority follows.

## Reproduction and next research gate

Retained roots: `tmp/pdbfixer-research/dense-enrollment-01/`, `dense-plan-01/`,
`dense-score-01/`, and saved audit revisions. Research-selected structures live in
`dense-score-01/<case>/selected-candidate.{pdb,dat}`, not public bundles/GUI outputs.
Exclusive new destinations are required:

```bash
# Scientific interpreter; frozen parent contact evidence is required.
PYTHONPATH=src:. python scripts/research_rosetta_dense.py enroll --output <enrollment>
PYTHONPATH=src:. python scripts/research_rosetta_dense.py prepare --enrollment <enrollment> --output <plan>
# Existing user-authorized licensed PyRosetta interpreter:
PYTHONPATH=. python scripts/research_rosetta_joint.py score --plan <plan> --output <scores>
PYTHONPATH=src:. python scripts/audit_rosetta_dense.py --plan <plan> --scores <scores> --coarse-scores <original-contact-scores> --output <audit.json>
```

**626 focused tests passed in 26.84 s**; targeted Ruff passed. Synthetic tests
confirm exact coarse containment, old-policy refusal of new angles, frozen
observations, THR inversion rejection, unchanged anchored closure, complete matrix
requirements, factorized limitations, generated-only ownership and saved-candidate
replay without invoking a generator. Not a full repository suite or broad chemical
qualification.

Next gate: specify a descriptive reference-robust analysis on the **already saved matrices**,
retaining every reference choice and stating where a feasible rectangle does not
exist; do not resample, select a favorable reference or call this new independent
prospective qualification. New empirical sampling,
context/protonation calibration and independent-family/full-legacy comparison
remain separate. These results do not justify an approximate interaction graph,
production routing or replacing PDBFixer consumers.

Final saved-candidate readback reproduces byte-identically. Package mypy,
generated CLI/GUI checks, agent maps, old evidence and all 49 fixture hashes also
pass. Final artifact checks are retained in
`tmp/pdbfixer-research/validate-dense-final.py`.
