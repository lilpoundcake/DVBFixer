# Experimental structure coordinate completion: current handoff

Updated: 2026-11-17. This is the current continuation entry point for the retained
Mac reconstruction work. It supersedes conflicting next-action instructions in
the older handoffs, not their evidence, frozen requests or recorded failures.
This documentation update does not launch or authorize additional experiments.

**Current priority:** the agreed finite anchor-mobility calibration in section 3
supersedes the historical autonomous next-case continuation. Do not resume
unrelated reconstruction batches. The frozen **100-cell grid was launched on
2026-11-15** after four-case full-context preflight, adapter verification and a
green fast suite. All 100 cells and independent readback are now accounted for;
the first saved-panel trade-off map is complete. The original protein-adapter
gate flags have a scope mismatch and are not a scientific success-rate table;
use the separately bound corrected readback below. **The subsequent failure
investigation found polar-H collapse, nonstationary returned states and a steric
exclusion bug. The map is descriptive evidence of that executed setup, not an
unconfounded calibration of anchor-policy effectiveness.** Do not start a new
cohort or choose movement limits until those numerical/diagnostic issues have
separately declared controls.
3C7D/1HE8 remain blocked, not executed or replaced. Use run/readback receipts,
not process success or this status line, to establish individual outcomes.

## 1. Purpose and source data

We work with **experimentally determined structures deposited in public structural
databases**, principally the Protein Data Bank. The task is to complete absent
coordinates for residues whose identities and order are already known from the
deposited structure's sequence records and associated experimental construct
information. We are not designing new amino-acid sequences or new proteins.

The distinction is between a known sequence and incomplete experimental
coordinates: an experimentally studied construct can contain residues that were
not resolved in the deposited coordinate model. We estimate coordinates for those
residues while preserving the observed scaffold and relevant structure context.

Before completion, bind the accession, original file digest, assembly, chain and
insertion-code identities, construct sequence, sequence-to-coordinate mapping and
exact missing-residue mask. A database-reference difference, engineered deletion
or omitted domain is not automatically an unresolved coordinate gap. If the
deposited construct does not establish the requested residue identities, report
the ambiguity rather than inventing sequence or borrowing an unverified canonical
sequence.

**Coordinate completion is not sequence generation.** ML backends do generate
candidate atomic coordinates; these are computational proposals for the fixed
known sequence, not new experimental observations. Do not describe predicted
missing coordinates as experimentally measured or claim that no coordinate
generation occurs. Existing fields such as `generated_atoms` retain their spelling
and mean computationally completed atom ownership.

## 2. Tools and boundaries

The workflow uses existing structural-biology tools and publicly documented
modeling methods, including open-source tools where applicable. Research
comparators include RFD3, Protenix, Boltz, Protpardelle and MODELLER. Record the
actual tool, revision, checkpoint and license for each executed arm; not every
backend is applicable to every structure.

Open source code, downloadable weights and permission to redistribute them are
different properties. MODELLER has separate licensing requirements; checkpoint
terms must also be checked independently. Do not label the entire stack freely
licensed or treat a public structure accession as evidence of model accuracy.

Current scope:

- Complete known-sequence missing residues in retained experimental structures,
  focusing on loops/local omissions in a sufficiently observed scaffold.
- Preserve observed coordinates and source context; keep computational additions
  distinguishable from observations through atom-level provenance.
- Check structural consistency, geometry, source preservation and reproducibility.
- Do not introduce sequence design, mutations, novel sequence combinations,
  binding/activity optimization or biological-function engineering into this task.
- Domain-scale omissions, 7N8X and the four already deferred focused submodels
  remain deferred. Do not silently crop them into eligible tasks.
- Do not expand to synthesis, expression or biological testing. None is part of
  this coordinate-restoration handoff.

## 3. Current direction and settled state

### Test-set work is data curation, not a modeling run

Keep these operations distinct: reading source metadata, defining a test manifest,
downloading a deposited file, and executing coordinate completion. A request to
inspect or organize the test set does not authorize the subsequent operations.
For the retained list, reuse local metadata and source files by default.

The inspected collection helpers (`scripts/collect_real_gap_candidates.py`,
`scripts/acquire_real_gap_sources.py` and
`scripts/collect_additional_long_gap_candidates.py`) use experimental method,
observed/missing residue counts, gap spans and sequence-cluster grouping, with
deterministic selection and a pinned historical accession list. Their inspected
selection logic does not rank biological activity, pathogenicity or binding
performance. This code inspection does not establish suitability of every
individual accession or explain a particular assistant refusal.

For further test-set curation, document structural coverage, duplicate/family
grouping, construct-sequence evidence, data quality and reproducible split rules.
Preserve accession titles, organism/function annotations and source provenance;
do not conceal identities or relabel biological objects to avoid review. Do not
expand curation into selection for enhanced biological properties. If a specific
entry or operation needs separate assessment, state the concrete issue rather
than treating all public experimental structures as either automatically safe
or categorically prohibited.

**Use the existing retained list.** Source recruitment and model-size calibration
remain paused. The old 200-entry ambition, length quotas and large job matrices
are historical planning targets, not current execution requirements.

The user has removed a separate case-admission review phase as a prerequisite to
continuing restoration on that list. The 177 pending review labels are provenance,
not a new mandatory queue. This does not remove per-task sequence/identity,
context, supported-chemistry, numerical, resource or output-validation checks, nor
does it reinstate excluded/deferred structures.

At handoff creation, previously declared calculations were settled; this document
did not itself queue another restoration job. The user subsequently authorized
autonomous implementation; the additional retained outcomes are recorded below.
Input preparation accounted for 233 source tasks:
194 ready, 31 refused and eight time-blocked. These are preparation outcomes,
not 194 accepted restorations or a completed independent evaluation cohort.
No final evaluation split is frozen.

### Latest retained example: 1KTZ

- The private four-chain candidate contains 468 residues and **710 completed
  heavy atoms** within its frozen mask. Those additions are modeled, not measured.
- Independent readback confirms 2,986 original protein observation records and
  326 original water records unchanged, global L-Cα geometry and the reviewed
  missing THR stereocentre. The separately confirmed native target has 21 S–S
  pairs; the original 18 explicit source records remain preserved.
- Consistent full diagnostics report no new geometry findings in the current
  export, eight unchanged source warnings and **two real missing OXT errors at
  receptor B136/D136**. Those OXT atoms are outside the frozen 710-atom mask;
  the errors are not baseline findings and must remain visible.
- Receptor terminal chemistry is deferred to future user-facing preparation,
  not a blocker for the original restoration task or continuation on the retained
  list. Do not ask again about the already confirmed TGF-β3 redox target.
- Do not silently add OXT, caps, hydrogens, charges or protonation choices.
  The candidate is not chemically complete or MD-ready. Its result does not
  establish missing-coordinate accuracy, broad backend qualification or topology/
  trajectory correctness.

Evidence: [final candidate and provenance](../../tmp/mac-real-gap-continuation-2026-11-14/1KTZ-restored-heavy-export-attempt02/README.md),
[independent readback](../../tmp/mac-real-gap-continuation-2026-11-14/1KTZ-restored-heavy-export-attempt02/readback-receipt.json),
[consistent validation reproduction](../../tmp/mac-real-gap-continuation-2026-11-14/validation-consistency-review-attempt01/report.json).
These local artifacts are not guaranteed to arrive through Git; report missing
evidence rather than recreating it or assuming a pass.

### Autonomous full-context continuation: four additional cases

Historical authorization: the user explicitly requested continued implementation without confirmation at
each routine step. Select bounded next tasks from the retained list, declare
their inputs/backend/seeds/budget, execute, independently audit, and continue to
other eligible tasks after a case-specific refusal. Do not stop merely because
an old handoff did not name the next accession. This authorization does not
remove the chemistry/identity/context/resource checks or permit silent scope
changes, new recruitment, broad job matrices or public capability changes.

Four additional single-protein full-context MODELLER controls executed on Mac:

| Entry | Retained result |
|---|---|
|1CC6|Native fixation passes for all 3,374 prepared-protein/FAD/PHB/water atoms. Raw and generated-only derivatives remain rejected: four angle warnings after the original CG continuation, five new geometry findings after the separate standard-stereo continuation.|
|9M3K|Native fixation passes for all 2,670 prepared-protein/water atoms. A serial-alias materialization failure was retained and recovered from saved native bytes without another model. The separate standard-stereo derivative retains one new boundary-angle warning and is rejected.|
|1X0X|One seed-104 native proposal passes the heavy-coordinate contract without extra refinement: 23 completed atoms, all 2,917 prepared-protein/NAD/sulfate/water atoms fixed, global L and zero new full-geometry findings.|
|5HJI|One seed-104 native proposal passes the same contract without extra refinement: 18 completed atoms, all 3,116 prepared-protein/adenosine/water atoms fixed, global L and zero new full-geometry findings. Its 24 previously completed fixed-conditioning atoms are separate from those 18 new gap atoms.|

The four cases are development-exposed, not an independent evaluation cohort.
Source preparation and fixed conditioning remain distinct from the new gap atoms;
the canonical `.dat` sidecars own only the declared gap-completion atoms.
The successful candidates retain 190/253 source warnings respectively, with no
CLI errors; neither result establishes coordinate accuracy or MD readiness.
No threshold was relaxed, observed anchor released, new source downloaded or
public capability enabled. Failed raw stages and failed derivatives remain failed.

The private adapter carries all retained context into the native model, selects
only the exact generated atoms and skips documented full-model `to_iupac()` atom
swaps. Native method-call evidence and direct pre-materialization coordinate
audits are retained. This is not a production MODELLER default change. An offline
`to_iupac()` control on the older 6Z38 initialization produced 94 changed named
atoms rather than its final 90; it does not establish the exact historical call
sequence or retroactively pass 6Z38. Source serial aliases are resolved without
changing coordinates or any other observed record field; topology-reload metadata
comes from the saved native proposal, never by normalizing experimental occupancy.

Evidence:
- [Four-case proposal declarations and retained stages, first pair](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-modeller-controls-attempt01/plan.json),
  [second pair](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-modeller-controls-attempt02/plan.json).
- [Independent negative-case readback](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-stereo-independent-readback-attempt01/receipt.json).
- [Independent successful-case readback](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-next-independent-readback-attempt01/receipt.json).
- [Private review exports and provenance](../../tmp/mac-real-gap-continuation-2026-11-14/restored-heavy-review-attempt01/index.json).
- [Consolidated outcomes and development-exposure overlay](../../tmp/mac-real-gap-continuation-2026-11-14/autonomous-full-context-continuation-summary.json).
- [Software verification](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-modeller-controls-attempt02/software-verification-final.json):
  2,866 fast tests passed, three optional Torch/Biotite skips and 33 slow tests
  deselected; final focused tests, Ruff and agent-documentation checks pass.

Review exports are deliberately not public diffusion bundles or manual donor
approvals. Public CLI/API/Homology/GUI integration and broader qualification
remain separate work. Local evidence must be transferred explicitly with a
continuation; missing evidence is not a pass.

### Current agreed stage: finite anchor-mobility calibration

The immediate task is **calibrating observed-anchor mobility and collecting
evidence from unsuccessful restorations**. Pipeline modifications, a required
reconstruction-success percentage and broader qualification are deferred until
those data are available. The first deliverable is a map of trade-offs, not a
preselected winning threshold or a production-policy change.

Start with retained proposals and controls, then assess dependencies on a separate
prospectively selected cohort from the existing inventory. The later cohort has
no approved membership, size or split and is not automatically a held-out
validation set. The old 200-entry/80–40–80 targets remain superseded.

#### Saved-proposal panel

Use 4–6 structures and exactly one retained base proposal per structure. The
provisional expanded panel is:

| Entry | Internal gaps, residues | Context and outstanding input checks |
|---|---|---|
|3C7D|10|NAD and pyruvate: verify actual force-field participation, not just retained output records.|
|6Z38|9|Context includes 3,115 ANISOU records. Bind exact proposal/native transport lineage; the historical 90 changed named observations remain a failure.|
|1NON|4 × 18|Four target chains, 358 waters; verify exact joint proposal and full context.|
|2BE1|2 × 8, 2 × 10, 2 × 20|Target chains A/B each have 339 residues and D has eight; preserve D and all 26 waters.|
|1HE8|Historical ten-gap stress case; exact internal/terminal breakdown not yet bound|Verify exact retained proposal, source mapping and whole-entry metal/ligand context. Earlier protein-only experiments are not whole-context qualification.|
|1KTZ|2 × 18|Four protein chains, 326 waters, already confirmed 21-pair disulfide target; bind exact retained 710-added-atom candidate and lineage. Do not re-ask the redox decision.|

This is **not an executable file manifest** or a guarantee of six admitted cases.
Do not substitute a protein-only projection, delete context or silently replace
a blocked structure. Record demonstrated blockers and coverage limitations;
membership changes require an explicit amendment. Bind original/prepared source,
context, request, base proposal, generator, seed and prior continuations before
execution. One proposal per case does not establish generator-independent effects.

Main scientific cases must have meaningful internal gaps. 1X0X/5HJI are short
N-terminal technical controls, not main successes; 1CC6/9M3K also lack internal
gaps. 6ZH0 has a one-residue internal gap, but its remaining A1–A2 defect is at
the N-terminal reconstruction junction, not evidence about internal-loop mobility.
Preserve these outcomes without inflating the main panel.

Each cell starts from the **same base proposal**, never another cell's optimized
derivative. All declared generated atoms, including terminal additions already
in the proposal, remain jointly mobile. Report internal and terminal junctions
separately. Only internal-gap flanks define the observed mobile windows; there
is no spatial-window arm initially.

#### Agreed grid and movement contract

| Factor | Agreed values |
|---|---|
|Observed windows|Zero motion; one observed flanking residue on each side of each internal gap; two consecutive observed residues outward on each side. Union overlaps once; never bridge another gap or chain boundary.|
|Anchor atoms|Observed heavy atoms in declared residues, uniformly weighted; no B-factor/ANISOU/occupancy weighting. All observed heavy atoms outside the window remain fixed.|
|Envelopes, Å|Per residue: (N/CA/C/O backbone RMSD, Cα displacement, maximum heavy-atom displacement) = (0.10, 0.10, 0.50), (0.25, 0.25, 0.75), (0.50, 0.50, 1.25), (0.75, 0.75, 1.50).|
|Positional strength|1,000 / 10,000 / 100,000 kJ/(mol·nm²), with energy `0.5*k*displacement²`. Historical common-refiner base: 10,000.|
|Envelope mechanism|Existing one-sided penalty walls, fixed coefficient 1,000,000 kJ/(mol·nm²), separate from positional strength; strict independent final movement audit.|
|Reference frame|Original source observations, no fitting or ratcheting between phases. Disclose earlier source preparation separately.|
|Count|Two mobile windows × four envelopes × three strengths + one zero-motion control = **25 per proposal**, 100–150 for 4–6 admitted cases.|
|Successful zero control|Still execute all 25 cells to measure unnecessary-motion costs; research only, not production escalation behavior.|

Envelopes participate in optimization, not just post-hoc classification. The
walls are soft penalties, **not mathematically hard constraints**. Independently
check per-residue movement both before and after PDB serialization. Exceeding a
limit is a failure: no clipping, rounding away excursions or widening thresholds.
These are experimental engineering controls, not physical uncertainty estimates.

Preserve occupancy and ANISOU. Where an anchor moves, disclose its unchanged
source tensor as retained metadata, not a newly refined experimental tensor.
No tensor transformation or fabrication is authorized.

#### Agreed physical and numerical setup

Use one common refiner across proposals; retain proposal-generator identity as
a separate factor instead of tuning each backend independently.

| Setting | Agreed setup |
|---|---|
|Engine|OpenMM CPU, one thread; record actual versions/properties, no silent fallback.|
|Base force field|`amber19/protein.ff19SB.xml` + `amber19/tip3p.xml`; bind XML resources and all additional context parameters. These files alone do not parameterize arbitrary ligands/metals.|
|Environment|Nonperiodic `NoCutoff`, no added solvent box or implicit solvent; preserve original waters and relevant heterogens.|
|Bond treatment|No bond-length constraints, `rigidWater=False`; preserve declared covalent/disulfide graph.|
|Temporary chemistry|One disclosed chemical model and identical initial temporary hydrogens per case across cells. Existing helper uses pH 7.0; record actual variants, charges and H inventory. Do not infer ligand microspecies or metal coordination from that heuristic.|
|Preparation|Up to 100 bonded-geometry + 75 first soft-steric + 75 second soft-steric iterations, following the versioned existing conditioning schedule.|
|Main optimization|Up to **1,000 full-force-field iterations**; up to 1,250 across four phases. No restarts or automatic extensions.|
|Tolerance|Existing 10 kJ/(mol·nm), identical across cells. Normal API return is not proof of convergence.|
|Other forces|Declare/freeze existing chirality, trans-amide, peptide-angle terms and soft-core schedule; no tuning per cell or backend.|
|Initialization|One declared preparation seed and identical pre-optimization temporary particle state per case, reused across all cells; minimization, not MD sampling.|
|Resources|Sequential execution, 30 minutes per cell including preparation/all phases; retain existing RSS 8 GiB, host-pressure and 10 GiB evidence-cap safeguards. Check available storage before launch.|

Temporary H atoms are auxiliary model particles, not observed coordinates or
publicly owned gap atoms. Record temporary heavy-atom completions and their
force-field role; never silently extend output ownership or call them observed.
In particular, 1KTZ B136/D136 OXT chemistry stays deferred and visible. Preserve
all source records outside authorized coordinate changes in materialization.

#### Execution gates (protocol; executed readiness updates below)

1. Bind exact files/digests for sources, contexts, requests, proposals, native
   transport audits, generated masks and anchor lists. Keep prepared conditioning
   distinct from original observations. Freeze a finite machine-readable manifest.
2. Correct the private window selector: historical `sequence2` in
   `scripts/run_1he8_traced_mobility.py` takes a ±2 axis neighborhood around each
   direct anchor, wider than the agreed two flanking residues. Do not reuse it
   unchanged for this arm.
3. Before research optimization, preflight actual full-context topology and
   parameters: ligands, ions, waters, connections, particle identities, charges,
   variants and fixed/mobile masses. Output record preservation alone does not
   establish force-field participation. Verify usable PDBFixer resources, not
   merely successful import.
4. Reuse valid existing chemistry where supported. A new ligand-charge, metal or
   chemical-state model needs a separately agreed specification; never choose it
   based on reconstruction success. Record a blocker without deleting context.
   Supported cases need not wait for every other case, but each must pass its own
   frozen preflight before its grid starts.
5. Test the private adapter: exact windows, identical initial temporary states,
   unchanged references, phase budgets, wall participation, all-generated joint
   mobility, native fixation, context preservation and independent serialized
   readback. Preserve production defaults.
6. Freeze executed-code and resource hashes and provide a ready/blocked ledger
   before launch. This document does not authorize running unverified inputs or
   unsettled chemistry. At protocol freeze, no study variant had executed.

Source inspection covered `scripts/run_1he8_traced_mobility.py`,
`src/dvbfixer/model/diffusion/boundary_refinement.py`,
`src/dvbfixer/model/diffusion/geometry_conditioning.py` and
`src/dvbfixer/model/diffusion/anchor_policy.py`. Existing 1KTZ native-context
evidence provides a route to inspect, not a pass for the new anchor adapter.
Source inspection did not establish whole-panel parameterization readiness.

#### Required reporting and first-stage closure

Retain for every cell:

- Exact input/code bindings, ownership and anchor identities; full temporary
  system/topology, initial particle state, chemical choices and force definitions.
- Actual phase calls, iteration observations where available, elapsed time,
  resource stops, energies, final forces and active force groups. Do not invent
  an API-unexposed termination reason. Different restraint Hamiltonians do not
  provide directly comparable energy-based quality scores.
- Native fixation before source copying; unchanged outside-window records;
  per-residue backbone RMSD, Cα and maximum heavy displacement before and after
  serialization, in the source frame.
- Whole-complex/junction geometry, sidechain stereochemistry, global L, declared
  links, Ramachandran/chi12, steric findings including relevant context/waters,
  and changes to source contacts/secondary structure. Bind analysis methods and
  settings before execution; do not invent acceptance thresholds from outcomes.
- Unchanged source findings, introduced defects, terminal-chemistry limitations,
  generated motion and observed-anchor motion as separate fields. Retain raw and
  final failures; no accepted-only reporting or public publication.
- Distinct outcome flags for geometry failure, envelope violation, numerical
  failure, resource exhaustion, input/chemistry blocker and implementation error.
  Multiple flags may coexist. A timeout or code bug is not a physical mobility
  failure.

Report paired effects relative to the zero control and same base proposal:
junction improvements versus anchor motion, lost/new contacts, secondary-structure
changes and introduced defects. Show all cells and per-case/cross-case trade-offs.
Twenty-five derivatives are not 25 independent structures; repeated chains are
not independent biological samples. Geometric acceptance is not coordinate
accuracy, proof of an unthreaded path or MD readiness.

**First-stage completion means a complete map of outcomes**, not a minimum success
rate: every declared variant has completed or has a documented specific blocker;
independent reports and failures are retained; numerical/resource stops are
classified; no unresolved implementation errors remain; and the paired summary
states coverage limitations. A blocked variant is not an executed variant.
Finding no suitable regime is a valid negative result, not authorization for
undeclared retries, wider limits or changed geometry thresholds.

Review the map with the user before selecting a policy, modifying the pipeline or
starting the separately declared prospective cohort. No final validation cohort
or required reconstruction-success fraction is agreed.

#### Executed preparation: 2026-11-14

The user authorized starting the agreed work. Private tooling is now in
`scripts/anchor_calibration_preflight.py`, with an independent saved-particle audit
in `scripts/audit_anchor_calibration_preflight.py` and focused regressions in
`tests/test_anchor_calibration_preflight.py`. Production code/defaults are unchanged.

The first declaration exposed historical prepared-source differences; it did not
execute any worker. The second declaration explicitly replays those already-bound
preparations, without performing a new source repair: 52 previously moved protein
atoms in 6Z38 and 31 in 2BE1, with original coordinates retained separately.
Per-atom original/prepared coordinates and fixed-conditioning additions are in
each case's `historical-preparation.json`. These changes are not anchor motion
produced by this study. Native observation-fixation failures in the selected
source-restored MODELLER derivatives remain disclosed and are not retroactively
passed.

Five guarded force-field preflight workers completed; **zero calibration variants
executed**. These workers may perform the existing temporary PDBFixer/hydrogen
preparation, but stop before the first research minimization call. Process success
is not force-field readiness or scientific acceptance.

| Entry | Historical attempt02 preflight outcome |
|---|---|
|6Z38|Temporary system constructed; independent saved-coordinate/mass/nonbonded readback passed. 6,443 particles, 3,082 fixed heavy atoms in the two-residue preflight window. One auxiliary OXT at A432 is present only in the temporary FF model, not output ownership.|
|1NON|Temporary system constructed; independent particle readback passed. 12,482 particles, 5,250 fixed heavy atoms in the two-residue preflight window. Very large finite starting energy remains visible; this is not a geometry/convergence pass.|
|3C7D|OpenMM template matching stopped at terminal HIS (temporary residue index 403), missing an external C connection. NAD/pyruvate participation has not yet been established; this is the first encountered blocker, not proof they are parameterized.|
|2BE1|Template matching stopped at terminal VAL (temporary index 709), missing an external C connection. No automatic cap or additional output atom was introduced.|
|1KTZ|Template matching stopped at terminal ASP (temporary index 294), missing an external C connection. Deferred B136/D136 chemistry remains unchanged. The new mobile-anchor model has not passed its 21-pair graph check.|
|1HE8|Input assembly retains four LINK records and GNP/Mg context. The reused unlinked-context merger refuses it; no FF worker executed. This is an adapter/context limitation, not a failed reconstruction or demonstrated physical incompatibility.|

At that stage, only **two temporary systems** passed particle readback. No case was
declared grid-execution ready: shared-initial-state reuse across all cells,
cell-specific auxiliary forces, full movement/geometry/contact reporting and the
independent grid adapter remain to be implemented and verified. Do not promote
these two preflight snapshots to completed 25-cell studies.

Before changing terminal chemistry, ligand charges or the metal model, agree an
explicit temporary-force-field specification separately from output completion.
Do not bypass the failures with omitted waters/ligands, unreported terminal atoms,
`ignoreExternalBonds`, relaxed geometry or hidden topology changes. Continue
supported implementation work without treating blocked cases as executed seeds.

**Subsequent explicit user decision:** in the temporary FF copy only, model the
uncapped target-chain C termini as free carboxyl termini under the agreed pH-7
setup and allow the auxiliary OXT atoms needed for that model. Record exact
residue identities/additions. This does not extend source/output/`.dat` ownership
or settle final user-facing terminal chemistry, including 1KTZ B136/D136. Existing
terminal C covalent bonds must fail closed rather than be removed to fit this
model. Ligand microspecies/charges and GNP/Mg remain separate unresolved models.

Preparatory attempt03 repeated the same preflight with an immutable `protocol.md`
copy and pre-bound XML resources, preserving attempt02 and its readback. This
avoids making future verification depend on a living handoff file that is updated
after execution. No research variant was executed in either attempt.

**Attempt04, after the explicit temporary-terminal decision:** independent
particle and base-force-field readback now passes for **6Z38, 1NON, 2BE1 and
1KTZ**. 2BE1 has one auxiliary OXT at D26; 1KTZ has auxiliary OXT at B136/D136
and its actual bond graph matches the exact approved 21 S–S pairs. None of those
auxiliary atoms extends output ownership. 3C7D now reaches the distinct missing
NAD-template blocker; its NAD/pyruvate model is not established. 1HE8 remains
blocked at the linked GNP/Mg adapter/model gate.

`scripts/anchor_calibration_grid.py` implements cell-specific force/mass setup
from the captured **unfrozen base System plus one initial NPZ per case**, without
rerunning hydrogen placement for each cell. It retains all generated atoms jointly
mobile, writes phase traces, compares native and serialized anchor movements,
checks native fixation and outside-window physical records, records full
diagnostics and emits generated-only `DatRecord` sidecars. Implementation errors
stop execution rather than masquerading as scientific refusals.

Descriptive environment metrics are declared before optimization: 4.5 Å observed
heavy-atom residue proximity, excluding covalently connected residue pairs, and
MDAnalysis C3 assignments on continuous peptide segments of at least six residues,
using auxiliary FF backbone H (missing-H donors disabled). The latter is an
**intrasegment comparison**, not full interchain-sheet DSSP qualification; its
reference is the same initial completed proposal, not experimental missing atoms.
Neither metric introduces an acceptance threshold or selects a policy.

The serial executor now accepts an optional **shorter** per-job timeout. The grid
charges the measured shared-preparation time against each variant's 1,800-second
ceiling, so reusing preparation does not silently grant extra optimization time.
Historical/default jobs still have 1,800 seconds; extensions are refused.

Four supported cases permit a finite **100-cell** first series, with 3C7D/1HE8
retained in the blocked ledger, not replaced or counted as executed. Final
software verification and the frozen execution declaration precede launch; no
research variant had executed at that preparation checkpoint.

**Execution launched 2026-11-15:**
`tmp/anchor-calibration-2026-11-14/grid-attempt01/plan.json` freezes 25 cells each
for 6Z38, 1NON, 2BE1 and 1KTZ (100 total), in that order. The background serial
executor uses the existing host/RSS/disk safeguards and shortened per-cell time
ceilings accounting for shared preparation. No alternate proposal, resampling,
restart, threshold change or production-default change is authorized.

Software gate: the first fast-suite invocation had 2,896 passes, three optional
skips and six missing-CLI `PATH` failures. With the environment bin directory
added only for the child invocation, the final complete fast suite passed
**2,903 tests, three optional skips, zero failures/errors**. The initial failed
report remains retained, not relabeled. Final frozen-adapter verification also
reconstructed all four initial states, checked finite forces and agreement with
the original common-refiner force states without executing minimization.

The live run retains per-phase states/reporters, exact native/serialized mobility
and `.dat` readback, full diagnostics and descriptive environmental metrics.
`summary.json` is written on exit, including partial/resource/code-error exits;
unlaunched cells remain unlaunched. Implementation/invariant errors stop the
executor and require an explicit retained correction, not a hidden retry. Do not
claim the 100 cells completed from the existence of the plan or launch receipt.

**Execution completed (before independent reanalysis):** all 100 guarded jobs
returned with retained reports and no resource stop. Ninety-four reached the
final reporting stage; six stopped at a global Cα-chirality gate (one 1NON and
five 2BE1), with rejected phase states retained. No repair, repeat optimization
or additional sampling was launched.

Post-run inspection found an **evaluation-scope implementation error**: the
protein-only `assessment_adapter(request, source)` was given full-context PDBs.
This produced `outside-generated-identity-mismatch` for waters,
`pdb-metadata-consistency` for the full ANISOU-bearing materialization, and
`retained-explicit-link-unexpected` for the separately approved additional 1KTZ
S–S annotations. Do not interpret the original all-failed mask flags as 0/100
scientific restorations. Do not erase those original reports either.

`scripts/audit_anchor_calibration_grid.py` now performs a separately declared,
read-only audit with **zero optimization calls**:

- Bind retained outputs and preserve the original reports/coordinates.
- Independently audit the full inventory, native fixed coordinates, physical
  outside-window records, ANISOU/occupancy/annotations, actual serialized mobile
  coordinates, generated-only `.dat`, global L and the exact approved 21-pair
  1KTZ link target. Independently reproduce rejected native D-Cα states.
- Build a protein-request evaluation view with every fixed/generated/anchor
  coordinate equal to the actual final full PDB. Source-record restoration,
  fitting and hidden anchor-motion erasure are forbidden. The frozen request is
  unchanged; context is assessed separately rather than silently omitted.
- Recompute protein-mask geometry and native/serialized movement under unchanged
  limits, distinguishing internal-gap and terminal findings.
- Run complete CLI diagnose on the prepared baseline and each retained final PDB,
  with waters and the unchanged 0.6/0.9 Å steric thresholds. The worker's earlier
  chemistry/Ramachandran/chi12/steric subset is not a substitute for full CLI
  postflight; in particular, retain deferred 1KTZ B136/D136 OXT errors explicitly.

The original grid runner must not be reused for another cohort without fixing
its assessment scope. Its private current implementation is now corrected to
use a coordinate-identical protein evaluation view and explicitly calls the
chemistry/Ramachandran/chi12/sterics result a subset, not full CLI postflight.
Future declarations use revision `anchor-calibration-grid-v2-scope-corrected`.
**That revision has not run new optimization.** The executed v1 code archive,
inputs, all original failures and result JSON remain unchanged. No production
pipeline/default was modified. Focused post-readback tests: 48 passed; Ruff and
agent documentation checks pass. The earlier complete 2,903-pass suite is a
separate pre-execution result, not a new full-suite claim for the readback edits.

#### Completed first saved-panel map

Independent full-context/request-view readback completed without resampling or
optimization. All 94 returned final candidates preserve native fixed observations,
outside-window records/ANISOU/occupancy and generated-only ownership. The six
rejected phase states independently reproduce D-Cα; a complementary audit also
confirms their native outside-window observations remained fixed. No rejected D
state was reflected, exported as accepted or relabeled a pass.

| Entry | Declared/reported cells | Returned final candidates | Protein-mask geometry + motion passes | Independently reproduced D-Cα refusals |
|---|---:|---:|---:|---:|
|6Z38|25|25|0|0|
|1NON|25|24|0|1|
|2BE1|25|20|0|5|
|1KTZ|25|25|18|0|

Interpretation boundaries:

- **6Z38:** all cells have zero local findings touching internal generated
  residues; the local failures are dominated by ten N-terminal clash findings.
  Do not present this as 25 failures to repair the internal nine-residue junction.
  Full CLI also retains a terminal OXT error and other introduced findings.
- **1NON/2BE1:** returned candidates retain real internal-gap geometry/steric
  defects. Local internal-residue finding counts range 24–40 / 15–20 across the
  returned cells; none clears the restoration mask and movement gates. One saved
  proposal per structure does not establish that every alternative proposal or
  backend would fail.
- **1KTZ:** the zero-motion control already passes the declared protein mask.
  Seventeen of 24 mobile cells also pass mask+motion; seven exceed the strict
  movement envelopes. Observed proximity-contact loss increases by up to two
  relative to zero, with no measured intrasegment C3 changes. This is evidence
  about unnecessary-release cost on this case, not a universal anchor threshold.
- **No cell has zero introduced full-CLI findings.** For 1KTZ the two deferred
  B136/D136 OXT errors remain, plus two CLI `chain_break` warnings explicitly
  describing ASP136 → HOH201 boundaries. Those protein-to-water warnings are not
  evidence of unclosed peptide junctions; retain the raw diagnostic and track the
  diagnostic limitation rather than silently removing warnings or changing
  production code. Mask success is not whole-entry chemical completeness.

The finite first-stage deliverable is complete: **a 100-row paired trade-off map,
failure evidence, two blocked cases and explicit coverage/metric limitations**.
It is not a selected/validated optimal policy, successful recovery of all panel
members, coordinate-accuracy validation or MD qualification. Twenty-five cells
from one retained proposal are not independent biological samples. No pipeline
change, extra seed, wider envelope, new source or prospective cohort has been
launched on the basis of these results.

**Next decision belongs with the user:** review this map and choose the separately
declared prospective cohort and/or bounded pipeline experiments for the retained
internal-gap failures. Freeze that next scope and analysis contract before new
execution; do not revive the old large cohort or replace blocked cases silently.

#### Failure investigation: 2026-11-17, read-only findings

The user requested a refusal investigation, **not new reconstructions or loosened
thresholds**. `scripts/analyze_anchor_failures.py` scanned all 100 saved endpoints
and 11 witness histories (initial state plus four recorded phases). It creates no
OpenMM Context, performs no minimization, repairs no coordinates, and preserves
all original reports. Four focused regression tests pass; production code and
thresholds remain unchanged.

| Case | Concrete saved-state evidence | Interpretation |
|---|---|---|
|6Z38|The initial proposal has seven severe generated-heavy contacts. The fixed-anchor endpoint retains ten; generated-heavy RMS force is 13.34 kJ/(mol·nm), while the massless-zeroed all-component diagnostic is 4.48. THR A1 Cβ orientation is +23.90° initially and +58.85° after phase 3, then opposite canonical R after the full-FF phase (−22.35°).|A bad terminal proposal remains in a geometrically bad endpoint, and the refiner introduces a sidechain stereocentre defect. The Cα wall does not protect THR Cβ. This is not an internal-anchor-window failure.|
|1NON|Initial actual-topology control finds 82 severe generated-heavy contacts; the closest is generated THR C76:CG2 to observed HOH C254:O at 0.269 Å. Fifteen of 24 returned final candidates have masked returned-force RMS above the configured scale of 10. Several show OH HG1 collapsing toward the same THR's CG2, including a 1–4 electrostatic exception with zero LJ epsilon. Some other endpoints meet the global force scale but retain real loop defects.|Mixed failure: deeply colliding initial proposals, numerical polar-H pathologies, and residual bad geometry at the aggregate force criterion. Do not treat every reported endpoint as a meaningful local minimum or infer that more anchor motion is the missing ingredient.|
|2BE1|All 25 last full-FF states have very large mobile residual forces; all 20 reported final candidates are above the configured force scale. SER B161:HG is the dominant mobile-force site in every cell and approaches HOH B12:O to 0.124 Å or less (0.002877 Å in `w1-e3-k10000`). HG has charge +0.4275e and zero LJ epsilon; the water oxygen has −0.834e and no exception for this pair.|Unphysical polar-H collapse in the selected unconstrained full-FF setup, not evidence that the internal loops are inherently unreconstructable. A normal minimizer API return did not establish convergence.|
|1KTZ|No newly exposed severe final contacts in the topology-only control; returned-force diagnostic is below 10 in every cell, as in 6Z38. The existing 18 mask+motion passes remain unchanged by this particular control.|Useful zero-motion positive control, not a universal policy calibration or chemical-completeness/MD pass.|

The force diagnostic zeros massless-particle components and reports their RMS
over the full component axis, separately from active-particle and generated-heavy
RMS. Comparing it with the configured force scale is not an exposed OpenMM
termination code or proof of convergence; a future numerical guard must match the
engine's actual normalization. The enormous 2BE1 gradients are not marginal
rounding/tolerance excursions.

**Concrete fixed-context contradiction in 2BE1.** SER B161:OG and HOH B12:O
are 1.658382 Å apart, and both are outside the declared mobile windows. Direct
Gemmi readback of the retained deposited mmCIF confirms those exact coordinate
vectors in model 1, occupancy 1.00, no alternate location; neither the historical
31-atom preparation nor this grid created the contact. This establishes the
contact's provenance, not a new chemical interpretation of it. In the declared
SER + TIP3P-water model there is no O–O covalent bond and that contact cannot be
removed while both observed heavy atoms remain fixed. Any source-context change
needs its own explicit preparation/chemistry policy; do not delete or move the
water silently.

**Confirmed diagnostic exclusion bug.**
`diagnose.steric._build_neighbor_graph` currently augments its covalent graph with
every sufficiently short distance (1.90 Å for ordinary heavy-atom pairs), even
between canonical protein and water. It therefore invents the SER–water O–O edge,
then `_expand_exclusions` hides that pair and related 1–3/1–4 contacts. This is not
a valid H-bond-envelope exemption: 1.658 Å is below that envelope, and the actual
FF graph/annotations do not contain the bond. The existing 0.9 Å error overlap
would flag 1.382 Å overlap if the fake exclusion were absent.

A private **same-threshold, actual-FF-topology control** exposed at least one
hidden severe fixed contact in all 25 2BE1 endpoints and additional hidden severe
generated contacts in six states. It exposed no hidden severe final contacts in
the other three cases. In raw witness proposals it exposes 59 previously hidden
contacts in 1NON and three in 6Z38. This control retains actual 1–2/1–3/1–4,
same-residue, water-inclusion and heavy-H-bond rules; it does not tighten the
0.6/0.9 Å thresholds or silently replace the historical reports. Production clash
logic has **not** been changed by this investigation.

**Proposal/context mismatch, independently bound.** The retained seed-101
MODELLER template files for 1NON and 2BE1 contain no HETATM/water records and
match their historical declared source digests. They were not conditioned on the
explicit crystallographic waters later attached for common refinement. The very
large initial generated–water overlaps are therefore not evidence that the
generator already solved the full-context problem. This is a property of these
retained proposals/routes, not a claim about all MODELLER capabilities.

**D-Cα localization.** All six recorded D-Cα failures first occur in the full-FF
phase: LYS B78 in 1NON and VAL B261 in five 2BE1 cells. The initial and preceding
recorded phase states are L. All six centres are present in the active Cα wall;
the wall is finite, not an omitted force or a hard feasibility constraint. The
H collapse and inversions coexist in that phase, but endpoint inspection alone
does not establish that one caused the other. The post-phase rejection correctly
kept all six states rejected; no reflection or repair has been performed.

**What this changes about the earlier interpretation:** the original 18/100
count is retained, not inflated by dropping failures. But 2BE1 and many 1NON
cells are confounded by numerical/model failures, while 6Z38 is a terminal
counterexample and 1KTZ already passes at zero motion. These data cannot select
a physically justified universal anchor threshold. In particular, merely raising
the iteration budget or widening envelopes is not a supported remedy for an
unphysical near-singular H attraction.

Proposed next work, **not executed here**: first correct topology-based clash
exclusions and reject/report incompatible immutable source contacts; then add
explicit returned-state numerical/H-geometry guards and test a separately
declared H-safe refiner setup (e.g. FF-equilibrium constraints on covalent X–H
bonds, not fixed noncovalent H-bond distances or frozen raw sampled heavy-bond
lengths). Retest a small predeclared saved-proposal control
set before another 100-cell grid or a new cohort. Proposal regeneration with the
actual water context and THR sidechain stereochemistry protection are separate
changes to evaluate, not automatic post-hoc repairs of the frozen study.

Failure evidence:
- [All-100 force/contact scan summary](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/summary.json).
- [Exact per-cell forces, attractive pairs and hidden contacts](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/endpoint-scan.json).
- [2BE1 full phase history](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/2BE1-w1-e3-k10000-history.json).
- [1NON fixed-anchor history](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/1NON-fixed-history.json).
- [6Z38 terminal/stereochemistry history](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/6Z38-fixed-history.json).
- [Fixed-contact diagnostic mechanism and preparation provenance](../../tmp/anchor-failure-analysis-2026-11-17/source-contact-provenance.json).
- [Deposited mmCIF contact proof](../../tmp/anchor-failure-analysis-2026-11-17/deposited-contact-proof.json).
- [Historical protein-only proposal-template proof](../../tmp/anchor-failure-analysis-2026-11-17/generator-context-provenance.json).
- [Exact saved-state THR A1 Cβ sign-change proof](../../tmp/anchor-failure-analysis-2026-11-17/thr-beta-proof.json).

Retained local evidence (inside the repository working directory; transfer these
artifacts explicitly for another host):

- [Finite input/cell inventory](../../tmp/anchor-calibration-2026-11-14/preflight-attempt02/inventory.json).
- [Guarded declaration and exact implementation/input hashes](../../tmp/anchor-calibration-2026-11-14/preflight-attempt02/plan.json).
- [Worker outcome summary](../../tmp/anchor-calibration-2026-11-14/preflight-attempt02/preflight-summary.json).
- [Independent readback and OpenMM XML resource hashes](../../tmp/anchor-calibration-2026-11-14/preflight-attempt02/independent-readback/receipt.json).
- [Explicit-terminal preflight outcomes](../../tmp/anchor-calibration-2026-11-14/preflight-attempt04/preflight-summary.json).
- [Four-case independent particle/base-force-field readback](../../tmp/anchor-calibration-2026-11-14/preflight-attempt04/independent-readback-v2/receipt.json).
- [Frozen 100-cell execution manifest](../../tmp/anchor-calibration-2026-11-14/grid-attempt01/plan.json).
- [Final no-minimization adapter verification](../../tmp/anchor-calibration-2026-11-14/adapter-final-verification/receipt.json).
- [Successful complete fast-suite report](../../tmp/anchor-calibration-2026-11-14/fast-tests-path-fixed.xml).
- [Original complete grid ledger (scope caveat above)](../../tmp/anchor-calibration-2026-11-14/grid-attempt01/summary.json).
- [Corrected independent 100-cell readback](../../tmp/anchor-calibration-2026-11-14/grid-readback-attempt01/receipt.json).
- [Rejected-state native fixation and D-Cα confirmation](../../tmp/anchor-calibration-2026-11-14/rejected-state-readback-attempt01/receipt.json).
- [Paired trade-off map, 100 CSV rows](../../tmp/anchor-calibration-2026-11-14/tradeoffs-attempt01/tradeoffs.csv).
- [Bound map summary and limitations](../../tmp/anchor-calibration-2026-11-14/tradeoffs-attempt01/summary.json).

## 4. Continuation contract

1. Inspect repository changes and retained evidence. Preserve uncommitted work,
   original source bytes, frozen requests and rejected attempts. Start from the
   retained list, not a new collection campaign or a compulsory review backlog.
2. For each next task, identify the experimental accession, exact known construct
   sequence, missing residues, complete relevant context and computationally
   movable atoms. Resolve task-specific ambiguity without rewriting history.
3. Bind any new execution to an explicit input/backend/seed/budget declaration
   under the user's actual execution authorization. Do not infer blanket approval
   from an old autonomous-run note or from this documentation request.
4. Preserve native/raw outputs before source-frame materialization. Independently
   audit actual observation fixation before restoring source records; copying
   observations back does not prove that a backend kept them fixed.
5. Change only the declared missing-coordinate atoms. Any permitted observed-anchor
   motion requires a separate explicit scope and movement audit. Keep source links,
   waters and other relevant context; do not silently crop a frozen request.
6. Apply consistent full diagnostics to the original source and completed result:
   inventory, junctions, stereochemistry including sidechain centres, global
   L-chirality, declared links, Ramachandran/chi12 and steric findings. Separate
   unchanged source findings, newly introduced defects and later preparation needs.
   A narrow-validator pass is not a clean full-diagnostic result.
7. Retain native/raw, materialized and refined stages and all failures. Use
   `DatRecord` for exact completed-atom ownership. Do not tune historical thresholds
   to accept a result or infer accuracy from a geometry pass.

Continue on Mac; transfer to Linux/CUDA is not the current objective. Preserve
declared host/resource limits, one-ML-process concurrency, native-MPS labeling and
no silent CPU fallback. Use `--no-solvent` for development `minimize`/`zbs` runs.
Do not start new downloads, remote sequence uploads, resource-limit calibration,
checkpoint changes or broad job matrices under this handoff alone.

Public CLI/API/Homology/GUI capability changes, production-default changes and
broader scientific qualification remain separate work. The original handoff
recorded 2,848 fast tests passed, three optional skips and 33 slow tests deselected;
retain that historical result separately from the current continuation verification
linked above.

## 5. Historical evidence and implementation references

- [Real-gap history and work packages](mac-reconstruction-real-gaps.md): retained
  inventories, experiments and older proposals; current directions above take
  precedence over contradictory historical next steps.
- [Original multi-backend portfolio design](mac-reconstruction-portfolios-handoff.md):
  comparison methodology, candidate lineage and earlier bounded studies.
- [Private RFD3 source-link custody](rfd3-source-link-custody.md): native graph
  transport, link provenance and derivative checks. A custody pass does not erase
  the original native graph-transport failure.
- [Dated execution log](../research/mac-real-gap-continuation-2026-11-12.md):
  chronology, negative outcomes and readback references.
- [Structure-preparation context](../agent/contexts/structure-preparation.md):
  code ownership, contracts and entry points.

Historical approval, refusal, launch and pending-decision statements describe
their own stages. Keep them readable; do not treat each as a new current command.
