# WP2K: graph-fragment torsion grid and retained-context steric ranking

Date: 2026-11-20; baseline `19025c8`, persistent
`feature/chemical-reconstruction`. Freeze before new conformer search/controls.

## Scope and documentation plan

Private, nonempirical heavy-conformer baseline on the SAME twenty observed
source/state declarations and whole-sidechain masks as WP2J. No new input source,
empirical rotamer probabilities, native-H model, production integration or licensed
data. PyRosetta research permission is not production permission; this route uses
none of its data. Native-reference coordinates cannot guide generation/ranking.

Documentation: this protocol; a pre-search baseline-signature ledger derived from
retained WP2J outputs; validation JSON; active handoff and task map. Commit and pin
the ledger before search, without changing its signatures after outcomes. Pending
links: none. User's continued implementation request authorizes this bounded slice.

## Immutable baseline and declarations

Rebuild the SAME source requests using unchanged WP2J source/state selection,
metadata preparation, full retained complex, opaque fractional context and CB-cone/
chain constructor. Preserve ASH/GLH/HIP as explicit experimental states, PRO's
explicit leaving H, GLY not applicable and missing eligible interior MET. Neither
the source nor generated masks may be substituted. Bind the baseline's rounded
generated-coordinate signature (identities/elements/coordinates only, no absolute
path metadata) to the pre-search ledger. No native coordinates in that ledger.

Start only from a complete baseline that passed raw AND rounded existing geometry.
Do not use rotations to conceal failed source/chemical/base geometry. Preserve
every original WP2J score and refusal; known-case repeats are not independent folds.

## Graph-defined axes and rigid downstream fragments

Compile axes from the exact projected target graph, not per-residue chi tables.
Use the existing sanitized complete chemical witness (virtual ALA chemistry only)
to inspect RDKit bond conjugation/hybridization. No witness coordinate becomes an
observation or output. All requested target heavy atoms must be present as fixed
observations or generated mask, no measured target H/D, temporary/mobile ownership
or source movement. The caller supplies observed root identities; the cohort uses
N/CA/C/O. Keep the generic graph operation separate from this protein adapter.

Eligible axis: heavy single nonaromatic, nonconjugated bond; endpoints must be
SP2/SP3 with at least one SP3. Deleting that bond must split the target graph; all
roots remain upstream. The downstream heavy component, excluding the axis endpoint
itself, must be nonempty and wholly generated. No moving observed fragment or
generated external link. Cyclic axes, both-SP2/conjugated bonds and empty terminal
heavy motions are excluded with explicit reasons. PRO ring puckering is not sampled.
Orient axes away from roots, order by upstream root distance then exact identities.
Choose heavy torsion anchors lexically among appropriate upstream/downstream
neighbours using nondegenerate source authority frames, not experimental holdouts.

Grid = Cartesian product of offsets `(0, -120, +120)` degrees relative to the
unchanged baseline torsions, parent axes before children. Maximum 256 combinations;
refuse rather than truncate a larger grid. All-zero is index 0 and preserves the
baseline addition coordinates exactly. No probability/prior, optimizer, random seed,
continuous phase search, reference fitting or success-dependent grid refinement.

Use RDKit's existing Get/SetDihedralDeg owner to rotate complete downstream
fragments, preserving branches, stereo and cyclic substructures rigidly. Populate
the scratch witness conformer with actual target observed/generated heavy points;
only requested generated heavy atoms can be returned. After each nonzero operation
require every target observation bit-identical; reject, NEVER restore, any moved
observation. Nonmoving generated points/axis endpoints may differ by at most 1e-10 A
numerically. No H from the scratch chemical witness is a generated/qualified H.

## Separate construction and independent audit

Add a distinct digest-bound search plan/strategy and dedicated entry point. The
ordinary single-proposal constructor must reject the search strategy rather than
silently ignoring its grid. Retain baseline, all attempts/offsets, raw/rounded
coordinates, geometry audits, trace audits, scores, failures and selection.

Reuse the existing candidate auditor without loosening WP2J heavy geometry:
0.10 A / 12 degrees, separation 0.7 A, heavy planes 0.05 A and stereo minimum 0.05.
Independent trace measurement checks ALL nondegenerate local heavy four-atom
paths involving generated atoms: sampled central bonds must have their declared
offset from baseline; all unsampled central bonds remain unchanged. Raw trace
tolerance 1e-6 degrees, rounded 0.25 degrees. These are provenance/transport checks,
not wider bond/angle/plane acceptance bounds. Auditing must not call the rotator.

## Frozen retained-context ranking, not an energy or extra acceptance gate

Use only existing `diagnose.geometry.VDW_RADII_ANGSTROM` radii, with its file/content
pin. Score every unique nonbonded pair with at least one generated heavy atom,
including other generated atoms and ALL retained observations. Exclude 1-2/1-3
paths from authoritative target/resolved graphs plus explicit covalent links;
do not invent unknown context bonds. Fractional opaque atoms remain included
without occupancy downweighting. Unknown element radii refuse ranking, never drop
an atom. Maximum 1,000,000 eligible pairs; no distance/context crop to fit budget.

For clearance `d - radius_i - radius_j`, report sum of squared negative clearances
and worst overlap (A^2/A). This is a steric-overlap heuristic, NOT LJ/electrostatic/
H-bond/solvation energy or calibrated posterior. There is no new VDW pass threshold.
Only attempts whose raw AND rounded geometry AND trace audits pass are selectable.
Choose smallest ROUNDED sum of squared overlaps, then smallest rounded worst
overlap, then enumeration index. Preserve raw/rounded scores independently. Native
reference metrics must never be consulted to break ties or rescue an attempt.

## Comparison, tests and remaining gates

Twenty observed entries, not twenty native successful executions: GLY N/A and
ineligible interior MET remain. Retain all candidates even if quality degrades.
Only AFTER selection extract the same exact-name native heavy holdouts and score
baseline/selected/all complete candidates using the existing scorer. Report paired
baseline/selected RMSDs separately from hard geometry and steric scores, no required
improvement count, post-hoc threshold or symmetry relabelling. Do not reclassify
WP2J's poor native TRP or any older failure. No claim that this grid is Dunbrack.

Test graph/ring/amide/observed-anchor exclusions; full multi-chi masks and rigid
branches/rings; proper rigid equivariance; exact fixed observations; complete
enumeration/budget refusal; no-reference construction/ranking interfaces; trace
corruption; deterministic ranking and pair accounting; unknown-radii/size failures;
baseline identity and old default rejection/replays. Independently rehydrate saved
plans/attempts, remeasure geometry/all torsions/sterics/selection/holdouts without
construction. Replay old controls byte-for-byte as regressions only. Run focused,
split slow boundaries, Ruff, strict selected/source-bound typing, fixture/maps/
generated checks and full fast tests; keep legacy typing and WP2F timeout visible.

Heavy/H integration, native-H/name quality, empirical backbone-conditioned priors,
continuous orientation/refinement, broader complex/linked/metal chemistry and
public/default/production/dependency/license gates remain OPEN regardless of scores.
