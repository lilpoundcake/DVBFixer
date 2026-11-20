# WP2P: prospective independent-source PRO mask comparison

Date: 2026-11-21; baseline `a1c1094`; persistent
`feature/chemical-reconstruction`. Freeze this protocol and source/mask manifest
BEFORE new admission, proposal, projection, PDBFixer or held-out scoring controls.

## Scope, authorization and documentation plan

User approved beginning the proposed independent comparison after WP2O. This is
a bounded PRO missing-heavy study on three source systems NOT used by WP2M/O,
not all-residue or broad complex qualification. The sources have appeared in older
research; they are independent of the two PRO development cases, NOT unseen data
or a blinded benchmark. Do not claim external validation or select successful sites.

Documentation: this protocol, adjacent frozen source/mask manifest, validation
record, active handoff/task discovery and a local quality note. Pending links: none
until the validation record is created. No version bump beyond 0.10.2, main merge,
public API/CLI/GUI/Homology/DatRecord/materialization, default/dependency change or
production/data-license approval. No new licensed data, rotamers or diffusion.

## Enrollment: nine cells, no outcome-dependent substitution

Use unchanged tracked raw 1TIM (enzyme dimer), 1CRN (disulfide-rich crambin), and
2CI2 (chymotrypsin inhibitor). Paths/SHA256 are frozen in
`chemical-reconstruction-comparison-cohort-v1.json`. They are distinct source
systems, not nine independent proteins. Preserve full retained context, waters,
heterogens, chain case, insertion codes and source annotations. No protein-only
crop, repair, altloc selection, occupancy rewrite or substitution to pass admission.

Select the FIRST file-order ATOM chain, then FIRST interior PRO with the existing
selector: 1TIM A/2, 1CRN A/5, 2CI2 I/25 (metadata-only discovery, no constructors or
geometry/native scores inspected). Independently run masks `CD`, `CG/CD`, and
`CB/CG/CD` on that SAME residue: 3 sources x 3 masks = NINE cells. Shared-system
and nested-mask dependence must be explicit. Never fall through to another PRO.
Absent/unsupported source or mask stays refused in the nine-cell denominator.
Natural additional missing target atoms must refuse the exact mask, not enlarge it.

Use existing explicit PRO authority and source preparation/probe. Before science,
independently audit original/prepared byte transport, exact mask and literal
incoming label/annotations, and the existing state-invariant type proof. Keep
HIS charges/states unassigned. Any raw source altloc/partial-occupancy limitation
remains a retained source refusal. Altloc-selected inputs require a later separate
declaration and are NOT silently added here.

## Three heavy arms and asymmetric auxiliary H evidence

1. **Original proposal:** one existing deterministic graph-internal proposal with
   the original heavy request/plan; retain partial/raw/rounded candidates even
   when its declared geometry rejects. No grid, template rerank or phase search.
2. **Projection -> H:** only after the existing baseline admission passes, one
   existing WP2O v2 real planar-port projection, same start/objective/cap/settings.
   If projection passes, exactly one existing complete-H proposal at the rounded
   heavy checkpoint. H cannot move/reselect heavy. Retain every refusal/failure.
3. **Seeded PDBFixer:** an explicitly separate raw heavy baseline, not a fabricated
   graph-internal/torsion/projection Candidate. Call the EXISTING
   rebuild_missing_atoms_with_retry helper (internal seeds 1..5), once per cell,
   on OpenMM Reference. These are NOT five independent proposals, nor production
   platform selection. Do not call addMissingAtoms directly, reflect, minimize,
   protonate, add H, strip context or repair to obtain acceptance.

PDBFixer receives ONLY a source-bound masked prepared PDB, never the withheld
reference. Filter only the exact masked atom records and annotation records
referencing them in a separately retained adapter artifact; record removed lines.
Keep all other atom/metadata bytes. Disable missing-residue insertion; scope
missingAtoms to the exact target names, disable unrequested terminal additions.
No other source incomplete residue is completed. Reject ambiguous topology keys,
unexpected/missing atoms or known-atom motion >1e-4 A. Measure and retain temporary
topology drift; never hide it by splicing an apparently fixed result. Run the
existing global L-chirality assertion; no D repair. Extract only exact additions
into its own DTO with PDBFixer/OpenMM/resource/platform provenance and rounding.
Retain raw additions and chirality failure if post-build global L validation fails.
Missing scientific dependencies/resources are environmental UNAVAILABLE, not
scientific rejection; unexpected program errors abort with existing evidence.

The old single-chain diffusion template baseline is NOT reused as an input gate:
its protein-only/single-chain admission is incompatible with retained complexes.
Reuse its existing seeded helper and exact-identity principles, not its crop or
diffusion request DTO. No full legacy Modeller/chirality-safe comparison claim.

H is auxiliary evidence ONLY for projection->H. Neither other arm constructs H;
do not compare their missing H as failures or claim H/native naming superiority.
No experimental H reference is present. Conditional ALA-witness CIP is not actual
neighbour/full-polymer chemistry. The primary three-way endpoint is heavy atoms.

## Common measurements, never constructor relabelling

Use a separate measurement-only additions interface factored from existing
candidate audits/held-out scoring. All three arms' raw/rounded coordinates are
measured under the SAME original+v2-port request, frozen heavy AuditPolicy and
constraint set, regardless of which constraints the constructor used. Do not
relabel PDBFixer as a graph strategy or feed failed external candidates to the
projector/H acceptance path. Keep each arm's own admission/audit separately.

Report complete exact-name fixed-frame heavy RMSD, maximum displacement and atom
count; no fitting, symmetry renaming, permutation, reference selection or native
coordinates in construction. Extract masked references ONLY after all available
three-arm proposal/H phases for that cell return. Score rejected but complete raw
and rounded proposals; incomplete outputs have NO partial score.

Common geometry includes bond/angle/explicit port/plane/stereo constraints and
generated-context separation with authoritative 1-2/1-3 exclusions. Report the
unique eligible generated-involving severe-overlap PAIRS below 0.7 A as a separate
measurement, with minimum eligible separation. This is the existing very-severe
overlap guard, NOT a vdW/MolProbity clashscore or force-field quality claim. Opaque
context's unknown internal bonds remain unassessed; do not fabricate exclusions.
Round to PDB 0.001 A and audit both phases. All observations fixed; no new output
writer or accepted public structure. Shared measurement extraction must preserve
every historical audit, digest, finding and replay byte.

Policies unchanged: heavy bonds 0.10 A/angles 12 degrees; H 0.02 A/5 degrees;
separation 0.7 A; heavy/H planes 0.05/0.005 A; stereo 0.05; projection movement
cap 2 A, one solver start, existing solver settings. No tuning on this cohort.
Native-reference geometry may be reported after proposals as a separate diagnostic
under the same bounds, not an admission filter or permission to relax thresholds.

## Accounting, readback and stopping rules

No required pass/win rate. Report per-cell raw/rounded results for every available
arm, source/proof/baseline/projection/H/template refusal and environmental failure.
Keep denominator nine and three source systems. Aggregate paired RMSD deltas only
on explicitly named complete pairs; list excluded cells/reasons. Also report
geometry-pass and RMSD-worsening counts without treating geometry as native quality.
Do not pool raw/rounded phases, nested masks, seeds, replays or repeats as independent
trials. No statistical/general superiority claim or completion of scientific gates.

Tests use tracked/synthetic source/mask fixtures and mocked rebuilds: masked-input
nonleakage, exact context bytes/identity, unavailable resource classification,
wrong/extra/missing/ambiguous atoms, fixed motion, global D refusal, no H/fallback,
distinct strategy provenance, common measurements, strict DTO readback and refusal
accounting. Independent saved readback reopens pinned sources/authorities and masked
artifacts, remeasures all additions/RMSD/geometry/overlaps/rounding/H without any
constructor, PDBFixer or optimizer. Tampering must fail. Retain private evidence
below `tmp/chemical-reconstruction-comparison-v1/`, not system /tmp.

Run focused/boundary/full-fast, Ruff, scoped/source-bound typing, maps/generated/
version/fixture checks and historical byte-identical replays. Preserve ambient
814 typing diagnostics, missing GUI dev dependencies, WP2F timeout, WP2L/v1 HIS/
altloc refusals and both WP2O native worsenings. If an environmental blocker prevents
an arm, keep it visible and do not silently use another prefix/platform/source.
Any later expansion/repair/model/threshold change needs a separately frozen protocol.
