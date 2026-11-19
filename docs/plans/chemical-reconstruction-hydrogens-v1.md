# WP2F: complete selected-residue H with fixed heavy observations

Date: 2026-11-19. Prospective protocol, baseline `a7fde14`; continue only on
`feature/chemical-reconstruction`. Freeze before new full-H control execution.

## Scope and reuse

Private complete-H proposal for one explicitly resolved internal peptide residue
in a retained source. Reuse WP2C exact inventory naming, WP2D pinned amide metrics,
WP2E source preparation/opaque-context policy and independent source audit, shared
ownership/admission, cone/chain/planar constructors and measurement-only audit.
The old observed adapter's H refusal and all original request/plan defaults remain
unchanged. This is a separate research route, not public partial-H or whole-source
protonation. Complete heavy inventory and no measured target H/D are required;
do not strip existing H, rebuild heavy atoms, choose a different residue/state or
silently generate only a convenient H subset.

Documentation plan: this protocol, a hydrogen validation JSON, active handoff and
agent task map. Pending links: none. No new source downloads, PyRosetta/Dunbrack
access, optimization, pH/state inference, publication or dependency replacement.

## Declared group/geometry policy

Every target H in the exact projected graph is generated once, every source atom
remains observed, and every other chain/heterogen/water is retained. Bind source,
preparation, inventory, complete H mask, authority and plan; no source-coordinate
payload or deleted experimental H reference enters the constructor.

Use these graph-admitted groups only, with observed-heavy construction anchors:

1. Backbone N-H: existing real incoming peptide C plus N/CA, pinned ff19SB
   N/H/XC types and WP2D metrics. Free-CCD N-H coordinates are not used.
2. Single-bonded carbon with three heavy neighbours/one H, or two heavy neighbours/
   two H: existing two-angle cone for each H, using the first two heavy neighbours
   in exact lexical atom-name order. Signed frame orientation comes from explicitly
   named CCD ideal H geometry, not an L/S shortcut. Keep BOTH fixed heavy stereo
   constraints and complete four-ligand CIP constraints. Extra H-name frames are
   labelled geometric conventions when the centre has no chemical stereo label.
   Independent sibling H-H angles, full CIP and clashes still must pass.
3. Single-bonded C with one heavy neighbour/three H, or +1 N with the same
   inventory: existing chain placement with one shared heavy frame, chosen from
   topology in lexical name order. Preserve source CCD phase/name correspondence;
   do not scan phases, optimize rotamers or relabel equivalent H after failure.
4. O with one single-bonded heavy neighbour/one H: the same deterministic heavy-
   chain frame, retaining CCD hydroxyl phase. No hydrogen-bond optimization claim.
5. Aromatic C with two heavy neighbours/one H: existing planar exterior placement,
   explicit scalar length/angle metrics from pinned CCD ideal geometry and an
   equal-priority geometric compromise, NOT force-field stiffness calibration.

Unlisted hybridization/inventory, absent/degenerate heavy frames or measured target
H refuses. Choose frames from topology and ideal authority only; never switch to
another frame based on observed success, generated clash or experimental reference.
All anchors are fixed observations. The graph's local explicit H count, formal
charges and atom/bond identities remain authoritative and unchanged.

## Independent checks and immutable policy

Introduce a separately digest-tagged generated-metric plan for different audit
scopes, preserving old plan digests. Fixed observed-heavy geometry keeps the WP2E
0.10 A / 12-degree bounds. Every bond/angle involving a generated atom must use
0.02 A / 5-degree bounds; both scopes use 0.7 A minimum separation and the same
0.02 A unused anchor bound. H-involving planes use 0.005 A; fixed heavy planes
remain at 0.05 A. Stereo normalized minimum 0.05. No observed atom moves and no
bound is selected after results. PDB rounding applies only to generated H.

Measure raw and rounded candidates independently, including H-H group geometry,
context clashes, every declared scalar and stereo/plane constraint. Also run
RDKit 3D CIP on the complete existing three-unit chemical witness with target
coordinates replaced by observed heavy + generated H. This is still CONDITIONAL
target CIP: virtual ALA partners are not actual source neighbours, and no full-
polymer chemical qualification follows. Keep geometry and CIP outcomes separate.
Retain every raw/rounded construction, geometry and CIP failure without reflection,
source repair, alternate phase, fallback engine or success-based mask reduction.

## Twelve prospective entries

Synthetic controls: five complete-H views of ALA/THR/ILE/LYS/PHE with projected CCD
heavy atoms observed and the WP2D synthetic incoming C at 121.96 degrees, plus
one ALA incoming-angle-100 control. All other target H must now be absent/generated,
not copied as observations as in WP2D. Hypotheses: five nominal ideal consistency
passes; the distorted case retains an amide-angle rejection. These are same-
authority synthetic consistency controls, not experimental H holdouts.

Observed controls (no atom withholding; generate ALL canonical target H):

| Source | Selection |
|---|---|
| `protein_context/1UBQ.pdb` | first file-order interior chain-A ALA |
| `protein_context/1UBQ.pdb` | first file-order interior chain-A ILE |
| `protein_context/1UBQ.pdb` | first file-order interior chain-A LYS |
| `protein_context/1UBQ.pdb` | exact chain-A PHE45 |
| `regressions/1FR2.pdb` | first file-order interior chain-A THR |
| `protein_context/1A3N.pdb` | first file-order interior chain-A ALA |

Selection uses identity/order only, never completeness or geometric success. Keep
admission/completeness/frame refusals counted. 1A3N is expected to retain the altloc
refusal. No observed-H success count is predicted. Sources/masks overlap known
research inputs; five targets from 1UBQ are NOT five independent protein folds.

Check frozen protocol, CCD/FF/naming and observed source pins before exclusive
output reservation. Save original/prepared bytes and full request/plan/results,
source audit, inventory/group ledger and raw/rounded CIP evidence. No candidate
PDB/.dat/public bundle. Re-read originals unchanged. Independently read back saved
requests/candidates and remeasure them without proposal construction.

## Regression and interpretation

Test complete H ownership, group/frame choices, same-parent H-H geometry, strict
generated versus fixed-heavy audit scope, proper rigid transforms, wrong/missing
authority/type/metric, fixed D-C-alpha refusal, measured-H/partial-heavy/partial-
occupancy refusal, source preservation, retained-context clashes, CIP corruption
and exact old defaults/refusals. Run focused heavy/H/source tests, Ruff, scoped
strict mypy and full fast suite; retain baseline typing failures separately.

This does not establish measured-H transport, arbitrary partial H, all20 chemistry,
native proton coordinates/rotamer accuracy, global hydrogen-bond quality, actual
neighbour/full-polymer chemistry, broad retained-complex heavy/H qualification or
production/public/default/dependency/data-license permission.
