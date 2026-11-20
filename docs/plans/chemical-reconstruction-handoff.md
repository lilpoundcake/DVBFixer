# Chemical reconstruction — active development handoff

Date: 2026-11-20. Status: WP0–WP2E, bounded WP2F complete-selected-H, WP2G
joint-sp3-H, WP2H planar-N/thiol and WP2I declared PRO/HIP controls checked;
WP2J all20 source coverage and WP2K fragment-grid/steric ranking are case-tested,
not native-quality qualification. Final fast suite passed (3754), split slow boundaries
passed. WP2F's large combined slow-inclusive run remains incomplete after outer
timeouts. General
reconstruction and scientific acceptance remain
incomplete. This document is the **starting point for
future reconstruction work**, not a claim that a universal tool has shipped.

**Current continuation:** bounded WP2L heavy→H composition, strict research
readback and independent combined auditing are implemented. Seventeen of eighteen
selected parents pass; PRO19's incoming-peptide-plane refusal is retained. This is
not native-H/general-complex qualification. See the final section and validation.

## 1. Persistent branch and authorization

- **Continue on `feature/chemical-reconstruction`.** It was created from
  `a7d41f9` (`research/openmm-upstream-report`), preserving the existing prototype,
  research scripts/tests, negative results and handoffs without cherry-picking.
- User explicitly requested one new long-lived branch, reuse of existing work
  and this handoff. **Do not create another branch/worktree for each experiment,
  milestone, handoff or session.** Use commits and versioned prospective protocols
  on this branch. Create/switch to a different development branch only if the user
  explicitly requests it. Inspect unrelated local edits before switching; never
  reset/stash/discard them automatically.
- Following the user's "start work" and subsequent "push and continue" instructions,
  WP0 `374158d` was pushed and WP1 continues on this same branch. No new scientific
  research cohort or production integration follows from this bounded implementation.
- No main merge, release, default-backend switch, dependency removal, installed
  OpenMM patch or automatic upstream publication follows from this handoff.

Original handoff setup scope (six files): this handoff and its validation JSON; discovery
updates in `AGENTS.md`, the historical PDBFixer handoff and the existing agent
task record; one narrow `.gitignore` entry for the retained QA artifacts. User's
direct handoff request authorizes these bounded changes.
WP0 documentation plan: shared-contract decision, this handoff update, task-map
discovery and a separate WP0 validation JSON, plus a narrow test-artifact ignore
entry. Pending links: none. The shared DTO/admission boundary now exists privately;
the coordinate engine and later architecture below remain proposed.

## 2. Actual objective and non-negotiable distinctions

The motivation is recurrent incorrect amino-acid stereochemistry during structure
preparation, plus hydrogen placement with invalid local geometry. The deliverable
is an own, maintainable, extensible reconstruction system — **not merely a
one/two-atom patch or another tool that rejects every nonprotein input**.

1. Reconstruct arbitrary admitted missing subsets, including **whole sidechains**;
   later cover missing backbone/true termini explicitly. Whole sequence gaps and
   loop modeling remain distinct owners, not implicit sidechain completion.
2. Make correct stereochemistry and H placement first-class requirements. Correct
   heavy-atom construction does not excuse invalid H or downstream damage.
3. Work against the **complete retained complex**: neighbouring chains, ligands,
   glycans, ions, solvent and covalent links cannot be removed just to pass.
4. Architecture must support nonprotein chemistry, not encode all chemistry as
   canonical protein residue names. Protein strategy can be the first adapter,
   but is not the architectural boundary or the final acceptance scope.
5. "Any chemistry" means extensible graph/authority/strategy support, NOT guessing
   an unknown molecule from incomplete coordinates. Missing or ambiguous chemical
   authority must be reported explicitly.
6. Distinguish **retain/read**, **assess missing atoms**, **reconstruct**,
   **assign protonation/place H**, **parameterize** and **minimize**. Passing one
   does not imply the others; an unassessed component is not "complete".

Primary acceptance is correct specified stereo, nondegenerate bond/angle geometry,
identity/observation preservation, clash handling and chemistry/defect coverage.
Held-out conformational accuracy and cost remain important secondary measures;
they must not replace the stereo objective or be used to select/tune an algorithm
after reference coordinates have been inspected. Our own code gives control, not
automatic reliability or lower maintenance cost.

## 3. Reuse map: code exists, broad qualification does not

| Existing owner | Reuse | Current limit / caution |
|---|---|---|
| [`domain/structure_identity.py`](../../src/dvbfixer/domain/structure_identity.py) | `ComponentInstanceRef`, `ExactAtomRef`; chain case, model, insertion code, altloc and occurrence | Do not collapse to chain/resnum or truncate identifiers |
| [`protein_completion/contracts.py`](../../src/dvbfixer/protein_completion/contracts.py) | Immutable source-bound observations, exact missing mask, additions/attempts/outcomes | Residue/name-oriented bounded protein contract, not a generic chemistry graph |
| [`protein_completion/geometry.py`](../../src/dvbfixer/protein_completion/geometry.py) | `l_cb`, internal-coordinate `extend`, two-anchor cone solver | Illustrative fixed geometry, NOT a calibrated universal geometry library |
| [`protein_completion/service.py`](../../src/dvbfixer/protein_completion/service.py), [`validation.py`](../../src/dvbfixer/protein_completion/validation.py) | Generated-only proposals; independently measured geometry, stereo and ownership rejection | ALA/SER/free CYS/VAL/THR scope with explicit missing-backbone, occupancy, metal and link refusals; not all20 amino acids |
| `protein_completion/{source,links,cohort,branched_source,anchored_source,mixed_cohort}.py` | Exact source projection, declared links, fixed-partner solving and bounded joint alternatives | Small supported cohorts; neither whole-complex chemical inventory nor a general packer |
| `protein_completion/{evidence,mixed_evidence}.py` | Retained failures, semantic replay, identity/digest checks and manifest-last publication | Preserve old kinds/revisions; do not widen existing readers silently |
| [`domain/nonprotein_reconstruction.py`](../../src/dvbfixer/domain/nonprotein_reconstruction.py) | Chemical graph/atom/bond/state/authority contracts, external links, explicit outcomes | Current graph data include ideal positions; not yet a complete shared constraint engine |
| [`nonprotein_reconstruction.py`](../../src/dvbfixer/nonprotein_reconstruction.py), [`ccd.py`](../../src/dvbfixer/ccd.py), [`pdb_component_reconstruction.py`](../../src/dvbfixer/pdb_component_reconstruction.py) | Pinned-local/user authority, Class A isolated-component reconstruction and exact source-preserving materialization | Existing opt-in capability; linked B/C chemistry, flexible closure and metal/redox solving are not covered; not automatic prepare fallback |
| [`ffutils/geometry.py`](../../src/dvbfixer/ffutils/geometry.py) | Existing stereo/H diagnostics, seeded rebuild wrapper and fail-closed errors | Emergency reflection/Cartesian H repair is not the new construction algorithm or proof of complete geometry |
| [`domain/force_field_naming.py`](../../src/dvbfixer/domain/force_field_naming.py), `force_field_naming.py`, `ffutils/variants.py` | Canonical variant/name and chemistry-state handling | Do not create a second AMBER/CHARMM alias implementation |
| `ffutils/{dat,ligand_valence}.py`, `domain/parameterization.py`, `structure_input.py` | Generated provenance, ligand chemistry knowledge, separate parameterization policy and CIF boundary | No hand-written `.dat` JSON, duplicated ligand rules or new ad-hoc CIF parser |

The source inventory and replacement obligations are already documented in
[PDBFixer usage/scope](../research/pdbfixer-usage-and-replacement-scope.md).
Its audited snapshot has11 constructors/nine seeded-rebuild sites, plus discovery,
mutation, state-carrier and readiness consumers. These are historical source
counts, NOT a fresh inventory made for this handoff. Refresh the existing scanner
when editing consumers; do not rediscover the entire research history first.

## 4. Preserved scientific findings and remaining gates

- Protein prototypes can preserve observations and construct useful constrained
  candidates. Their narrow geometry successes do not prove broad accuracy or
  chemistry coverage. On new inputs, library priors/means/ref2015 did not provide
  a reliable general ranking improvement; do not adopt them as established winners.
- Existing coarse clash gates are not a calibrated universal steric policy:
  [retained native 1–4 controls](../research/protein-completion-unfiltered.md)
  exposed false rejections. Reuse measurement/identity machinery, but version and
  independently qualify topology-aware contact rules; do not silently loosen old
  cutoffs or assume every old accepted/rejected geometry label is chemical truth.
- Full legacy comparisons and isolated downstream experiments show that structure
  defects can arise during H addition/relaxation, not only PDBFixer heavy repair.
  Removing PDBFixer alone is NOT a complete solution.
- [Compatible-noop post-H control](../research/post-hydrogen-local-noop.md): both
  saved states passed offline proposal audits, but **0/2 projected paths passed
  native all-node geometry**; HA286 became incompatible again. Transient negative
  Cα nodes occurred; final States had zero detected D-Cα. Normal C-H length and
  zero final D-Cα were insufficient. Preserve every refusal, bad node and original
  numerical envelope; no post-hoc reflection or threshold widening.
- [Synthetic OpenMM reproducer](../research/openmm-cpu-sign-reproducer.md): first
  oversized series negative; separately frozen representable-force trigger showed
  a reporter/objective gradient inconsistency with11 passing controls. This is
  not direct internal-fallback/binary attestation or proven cause of HA286.
- [Upstream readiness](../research/openmm-upstream-report-readiness.json): compact
  native entry subsequently verified once; issue remains UNFILED, requiring the
  user's own opening/review under OpenMM policy. This side task must not become
  the whole reconstruction project. No automatic maintainer replies or dependency
  patch. Existing scientific counts and later packaging counts stay separate.
- [Historical handoff](pdbfixer-reconstruction-handoff.md) is the chronology and
  evidence index, not the active instruction to create another research branch.
- Installed PyRosetta access is **research-only**. Original source/version/license
  permissions and explicit production approval remain mandatory before production
  integration/distribution of that route. Do not bundle tables/raw empirical
  arrays. Missing original tables do NOT block already authorized local research
  or table-independent contract/geometry development.

## 5. Target architecture — proposed, not implemented

```text
source + exact identity + explicit chemical authority
                  |
        graph/state/ownership request
                  |
      capability admission and explicit routing
                  |
   +--------------+----------------+----------------+
   | protein      | nonprotein     | H geometry     | further strategies
   | construction | graph/closure  | placement      | nucleotide/glycan/metal
   +--------------+----------------+----------------+
                  |
       shared geometry/constraint primitives
                  |
    independent chemistry/stereo/geometry audit
                  |
       source-preserving result and provenance
                  |
   separate parameterization / optional refinement
                  |
       stage/serialization/publication audits
```

### Shared contract requirements

- Complete intended graph and external links, including bond order/aromaticity,
  formal charge, explicit H/state and stereo constraints with authority provenance.
  Chemistry must not be guessed from a residue name, proximity or incomplete atoms.
- Exact source binding, explicit coordinate units, atom identity and separate sets
  for observed, requested-generated, temporarily materialized and explicitly mobile
  atoms. No raw topology indices across rebuilds. No silent promotion of temporary
  FF atoms to output or `.dat`.
- Constraints can express lengths/angles/planarity, specified tetrahedral centres,
  torsional freedom, ring/fragment closure and environment exclusions. L is a
  canonical-protein requirement, NOT a universal stereochemical label; ligand,
  sugar, nucleotide and intentionally non-L chemistry follow their authority.
- Fixed observations remain fixed. Inconsistent observations produce a retained
  invalid-input/ambiguity/refusal result, not an observed-atom reflection or a
  concealed restore after an optimizer moved them. Existing repair-enabled
  workflows require separately explicit movement masks and provenance.
- Distinguish exact source-record preservation from explicitly bounded native
  precision transport. Bind before/raw/returned coordinates and movement masks;
  allowed device rounding is not an intentional geometry edit, and restoring
  moved observations afterward cannot manufacture an observation-preserving run.
  Declare numerical/serialization bounds before execution, never after outcomes.
- Explicit operation-specific capabilities and statuses: unsupported, ambiguous,
  invalid source, no feasible geometry, exhausted search/numerics, unchanged and
  successful reconstruction cannot be merged into one fallback/"clean" state.
- Admission, generation, hard validation, ranking and refinement are separate.
  Reuse full-complex context; candidate scores cannot override stereo/identity
  failure. Preserve all candidates/refusals and bounded-search metadata.
- Version new graph contracts/adapter policies. Keep old `CompletionRequest`,
  Class A contracts, evidence kinds and replay semantics intact behind adapters.

WP0 uses the private `chemical_reconstruction` namespace, as specified in the
[shared-contract decision](chemical-reconstruction-contract-v1.md). It reuses
domain/nonprotein graph concepts rather than inventing a second ligand authority
or renaming the old residue contract "universal". This is a planning boundary,
NOT an implemented generic coordinate engine or publicly registered capability.

### Hydrogen requirements are a parallel workstream, not a final afterthought

Resolve protonation/tautomer/microspecies through existing chemical authorities
and PROPKA/Reduce/variant owners; geometry placement must not silently change the
chemical state or pH-adjust authoritative ligand SMILES. Build H using specified
sp2/sp3/amide geometry and correct atom counts/stereo. Search OH/methyl and other
admitted orientations in geometry-preserving degrees of freedom against the
retained complex. Observed H ownership and explicit re-protonation/removal are
distinct from "new H". Jointly check heavy/H tetrahedral geometry, not just Cα
heavy-volume sign or C-H distances. This is a requirement, not a qualified solver.

Cartesian minimization with finite penalty walls is not a mathematical guarantee
of stereo preservation. Future refinement must either preserve the admitted
geometric constraints by its parameterization or reject invalid returned outputs.
Do not silently invent new native forces, freeze production to Reference, remove
constraints, or assert the OpenMM issue is fixed. Validate after each public
materialization/refinement boundary and again after PDB rounding/publication.

## 6. Capability expansion and integration gates

| Area | Required development / acceptance before claiming support |
|---|---|
| Canonical proteins | All20 residues, complete and partial sidechains, multi-χ, THR/ILE centres, PRO/aromatic rings, disulfides and external anchors; no arbitrary one/two-atom limit |
| Backbone, termini and caps | Explicit peptide closure and true/temporary terminus ownership; OXT/ACE/NME; no fictitious termini at gaps; whole-gap modeling remains separately routed |
| H | Exact state/H count/naming, HA and sidechain stereo, sp2/sp3/amide geometry, permitted orientation search, source-H preservation and downstream validation |
| DNA/RNA | Authoritative nucleotide graph, sugar stereo/pucker, base planarity, phosphate and polymer links; preserving a chain is not rebuilding it |
| Isolated ligands | Reuse Class A authority/identity/publisher; add qualified flexible/closure strategies where required; no silent SMILES matching or state guesses |
| Glycans, PTMs, covalent ligands | Graph spans residue boundaries; anomers, linked anchors and ring closure; explicit Class B qualification, not a Class A fallback |
| Metals/cofactors | Declared coordination/redox authority and suitable strategy or explicit not-assessed/refusal; no generic tetrahedral protein rules |

During transition preserve valid complexes and current public behavior through
explicit adapters/routing. Unsupported *reconstruction* is not permission to drop
a component or claim it was completed. Predeclared fallback to an existing engine
must distinguish lack of capability from scientifically rejected geometry; never
catch arbitrary geometry errors and quietly retry different physics.

Do not remove PDBFixer until discovery, mutation/normalization, state carriers,
every materializing/read-only consumer, package/readiness and broader chemistry
obligations are replaced or explicitly routed and tested. Do not equate selecting
a preparation backend with a PDBFixer-free end-to-end pipeline.

## 7. Work packages and concrete next action

**WP0 — first coding slice; implemented and checked on this SAME branch.**

1. Read this handoff, existing agent research task, protein contracts and Class A
   domain/service. Rebind the existing usage inventory if touching call sites.
2. Write a short versioned shared-contract decision: which graph/authority concepts
   adapt existing domain types, which geometry constraints are new, their units,
   ownership, capabilities, statuses and source/evidence semantics. Proposed types
   must not be advertised as shipped. Do not change old schemas in place.
3. Implement the private pure contract/admission slice and explicit adapters from
   at least one supported protein request AND existing isolated-component graph.
   No coordinate generation, optimizer, network/CCD fetch or public registration
   merely to inspect/admit a request. Both heavy and H operation requirements must
   fit the contract from the start.
4. Tests: exact mixed-complex context and links survive; case/icode/occurrence and
   units stay distinct; graph/state ambiguity and conflicting fixed stereo refuse;
   old protein/Class A evidence readers behave unchanged; no residue-name-only
   graph fabrication; unsupported class cannot become a false success. Pure tests
   use synthetic authorities, not licensed PyRosetta data.
5. Commit the coherent contract/adapter/tests slice here and update this handoff.
   Do not start another native sampling cohort or branch to finish WP0.

**WP1 — shared geometry primitives / first cross-chemistry slice; bounded private
implementation checked.** Adapt tested
internal-coordinate/anchor primitives and a Class A example to the agreed graph
contract; independent validators must measure rather than call the generator.
Close source and rounded geometry before selecting conformers. Existing ideal
lengths/angles are provisional until a documented geometry authority is chosen.

**WP2 — expand protein reconstruction and H in parallel.** Whole/partial sidechains,
multi-χ and rings, then explicit backbone/terminal cases; authoritative geometry
and stereo. Geometry-preserving H construction and orientation solving are an
equal exit condition, not a promise after heavy integration. Keep empirical prior
and physical/context scoring separate from hard validity.
WP2A starts with pinned all-atom CCD authority and ideal roundtrip controls below;
it does not close polymer residue/state, distorted-source, ranking or general H gates.

**WP3 — cross-chemistry and whole-complex qualification.** Freeze new evidence
cohorts/capability strata before execution: protein, isolated nonprotein and linked
chemistry; exact source/candidate/movement accounting. Add nucleotide/glycan/metal
strategies under their own declared gates; do not substitute passage through the
router for successful reconstruction. Retain unsupported/refused cases.

**WP4 — opt-in production adapters, then dependency retirement.** Only after
stereo/H/coverage evidence and explicit approval, integrate consumers one at a
time, preserving variant/link/dat semantics and fallback policy. No default or
dependency change just because WP0/WP1 tests pass.

Maintain objective stop/go criteria before experiments: worthwhile improvement
in stereo/geometry reliability, observed-atom preservation, coverage or cost at
acceptable conformational quality. Keep negative results; do not tune on known
native RMSD, call synthetic coverage general accuracy, or bury side effects.

## 8. Checks, environment and resume discipline

Existing interpreter used in the retained research is the `dvbfixer` scientific
environment. Do not repair/install/swap it merely for this handoff; platform and
resource readiness must be verified for the specific next native protocol. No
permission-denied external directory discovery may be bypassed. Private raw
coordinates/RNG/Systems and licensed tables stay under ignored project paths.

Useful existing focused checks (run applicable tests, not every historical study):

```bash
PYTHONPATH=src:. python -m pytest -q tests/test_protein_completion.py \
  tests/test_protein_completion_source.py tests/test_protein_completion_links.py \
  tests/test_nonprotein_reconstruction.py tests/test_ccd_component_reconstruction.py \
  tests/test_pdbfixer_usage_inventory.py
python scripts/check_agent_docs.py
git diff --check
```

Heavy/H source changes additionally need relevant stereo/hydrogen and seeded-
rebuild regressions before the fast suite. Every dev `minimize`/`zbs` invocation
must include `--no-solvent`; no native runtime settings may be changed solely to
green a test. Last inherited full fast suite:3305 passed,3 optional Torch/Biotite
skips,33 slow deselected — historical unchanged-source result, not new generic-
engine qualification or a fresh suite run in this documentation task.

Continue the WP1/WP2 geometry and H workstream on this same branch. Do not restart completed PDBFixer/Rosetta/OpenMM investigations,
silently retry frozen scientific cohorts, or create per-phase branches. Report
actual test outcomes and update per-work-package status here; historical handoffs
remain evidence, while this document and the user's later instructions own the
forward development direction.

### Documentation completion checklist

- [x] Persistent branch created from the complete current work.
- [x] Scope, reuse map, limitations, architecture and first coding slice written.
- [x] AGENTS/task/historical-handoff discovery updated.
- [x] Local links, source preservation and focused checks verified.
- [x] Validation/quality record created; no pending links.

[Handoff validation/quality record](chemical-reconstruction-handoff-validation.json):
182 focused regressions passed with no skips; agent maps and whitespace checks
pass. The first invocation had170 passes/12 setup errors from a missing basetemp
parent; its JUnit is retained. Only the launch directory was corrected for the
second invocation; no source/scientific threshold/environment change. This is
documentation/regression verification, not new chemistry qualification. No new
native scientific cohort or production implementation was run for this handoff.

## 9. WP0 implementation continuation

The private [contract/admission decision](chemical-reconstruction-contract-v1.md)
now maps to executable DTOs, a pure admission function, independent fixed-constraint
measurement and explicit legacy protein/in-memory Class A adapters. No old APIs,
evidence readers, scientific solver, source parser or production consumers changed.

- All supplied observations and link endpoints remain bound to the request digest.
  Case, insertion code, occurrence, coordinate units and source SHA are distinct.
  Protein links without authoritative order remain unresolved instead of guessed.
- Admission distinguishes invalid masks/fixed constraints, ambiguous state/link
  chemistry, unassessed target chemistry and unsupported declared capability.
  Missing/mobile constraints defer explicitly; all readiness/generation flags stay
  false, including on contract admission. No full-source inventory claim.
- H inventory is explicit; new/temporary H cannot be admitted from an unspecified
  H inventory. Heavy, H and repack ownership remain separate. The synthetic
  12-missing-atom control tests contract breadth, not whole-sidechain construction.
- The first Class A compatibility adapter is intentionally the connected CHNO
  in-memory subset, not all legacy Class A chemistry. Non-CHNO organic graphs
  refuse at this adapter only. Generic DTOs/capability descriptors are not limited
  to CHNO, protein or one/two missing atoms. Extract a shared Class A classification
  policy before expanding this adapter; do not duplicate the old solver's rules.
- Pure-boundary tests guard imports and execution against scientific toolkit
  loading, data-file reads, subprocesses and network. Standard Python source
  discovery is permitted during imports only, not during admission/adaptation.

[WP0 validation record](chemical-reconstruction-wp0-validation.json): 61 new
synthetic controls; 335 final focused tests passed, including unchanged protein,
Class A, CCD, usage-inventory and evidence readers. Full fast run: 3366 passed,
three existing optional Torch/Biotite skips, 33 slow deselected. The final small
follow-up preserving Class A policy values explicitly in `Origin.settings` was
checked by the 335-test focused run and static checks, not another full-suite run.
Ruff passes for all source and the new tests; strict mypy passes for all five
new package modules. The broader existing mypy command reports 41 errors in 14
files, reproduced with byte-identical diagnostics on an unmodified `21a7162`
source archive. Those unrelated failures were neither fixed nor suppressed.

The first full invocation retained 3353 passes/six failures/three optional skips:
the agent omitted the selected environment's bin directory from PATH, so six
unchanged subprocess tests could not start `dvbfixer`. The corrected command only
added that existing bin directory to its child PATH; no package installation,
environment patch, threshold change or test skip. Both original and corrected
reports remain under the narrowly ignored `tmp/chemical-reconstruction-wp0-v1/`.
No new native reconstruction research cohort, production consumer or licensed
data was added. The historical handoff validation JSON remains unchanged.

WP1 starting decision (now implemented below): choose a documented geometry authority and adapt the existing protein
internal-coordinate and Class A construction primitives through the private graph
contract with independent raw/rounded audits. The current Cartesian measurement
module is not a constructor, full chemical validator, CIP implementation or
calibrated steric policy. Do not expose a generic CLI/GUI/default or reinterpret
these synthetic controls as native stereochemical/H qualification.

## 10. WP1 geometry continuation

The [geometry decision](chemical-reconstruction-geometry-v1.md) uses the selected
state's explicitly supplied ideal graph/provenance/units as the authority for
this bounded slice. This is data-driven geometry transport, not a calibrated
all20 geometry library. Tests supply labelled synthetic authorities, never held-out
native positions or licensed rotamer data. Real-chemistry data/qualification gates
remain open.

- `construction.propose` reuses `protein_completion.geometry.extend` and
  `two_angle_tip` through explicitly ordered plans, and the existing Class A
  `_fit_coordinates` rigid-fit owner through a separate lazy adapter. There is no
  coordinate optimizer, ranker, emergency reflection or engine fallback.
- Arbitrary admitted sequential missing subsets fit the plan: the synthetic
  12-atom chain now exercises construction, not merely WP0 mask admission. A
  separate explicit-H smoke uses the same authority; it is not general H chemical
  count/protonation/geometry qualification.
- The Class A route retains its connected CHNO, isolated, in-memory envelope,
  requires the original explicit pH metadata without using it to choose a state,
  and refuses any unowned missing graph atom. Temporary/mobile atoms and generated
  external-link endpoints require future explicit strategies and currently refuse.
- `candidate_audit.audit_candidate` measures rather than reconstructs: request/
  plan bindings, exact generated inventory, graph bond/angle authority, every
  declared constraint and severe generated/context separation. The whole shared
  request remains intact around the target-only rigid fit. Unknown context remains
  listed unassessed; no VDW/CIP/global-L/valence/MD capability is inferred.
- Both raw and generated-only PDB-precision-rounded proposals must pass. Observed
  coordinates, original units/order and context links are never rounded or moved.
  Wrong or partial proposals and both failed audits remain in the private result.
  No PDB bytes, workspace artifact, `.dat` or public bundle is written.

[WP1 validation record](chemical-reconstruction-wp1-validation.json): 34 new
synthetic geometry controls; 465 focused regressions passed with no skips, including
WP0, legacy protein/Class A/evidence and existing chirality/seeded-rebuild/H guards.
Strict mypy passes for all nine private package modules; Ruff passes for all source
and new tests. Full fast: 3400 passed, three existing optional Torch/Biotite skips,
33 slow deselected against the same frozen source snapshot. Agent maps, local
links, source/evidence hashes and whitespace checks pass; all four planned WP1
documentation/map files exist, with no pending links. No production source or
existing scientific threshold changed.
WP0's baseline-only 41 legacy mypy errors remain retained and out of scope.

Next development after this slice: independently sourced/documented real-chemistry
geometry authorities and prospective multi-residue/cross-chemistry controls; expand
whole/partial sidechains, closure and geometry-preserving H together under WP2.
Do not call a synthetic chain a realistic sidechain benchmark or enable production
because its hard geometry passes. Prior native post-H/ranking failures and all
license/source/evidence gates remain unchanged.

## 11. WP2A pinned CCD authority and chain compiler

The [protocol](chemical-reconstruction-ccd-v1.md) was committed as `67c99e8`
**before acquiring/inspecting the pinned control CIFs or executing controls**.
An earlier webfetch endpoint probe returned only an unsupported-content-type error,
not coordinates. It is immutable: the private runner
checks its SHA and the post-acquisition source-pin manifest SHA. Do not edit that
protocol to change outcomes; declare a separate version before future experiments.

- Seven unmodified wwPDB CCD source CIFs and per-file pins are tracked below
  `tests/fixtures/chemical_reconstruction/ccd/`. The existing CCD downloader/parser/
  local pin checker are reused, not a second CIF reader or runtime network fallback.
  Fixture provenance cites official CCD and CC0 usage-policy pages. PyRosetta/
  Dunbrack permissions and public/production gates remain unchanged.
- `ccd_to_all_atom_graph` is a data-conversion helper alongside the unchanged
  heavy-only converter, not a public reconstruction workflow. All source H,
  leaving atoms, explicit bonds/charges and ideal positions are retained. Free
  amino acids keep OXT/HXT/H2; **these are not peptide-residue state authorities**.
- The private `ccd_authority.load_geometry_authority` requires an actual SHA pin
  even at runtime, then sanitizes a temporary explicit-H RDKit graph. Charge changes
  or implicit-H additions during RDKit cleanup, unsupported radicals, unrepresented
  bond stereo and CCD-versus-3D CIP contradictions refuse. There is no state/pH
  selection, AddHs, experimental-coordinate import or silent source repair.
- CIP agreement is established before deriving ordered stereo constraints. Full-H
  controls also receive four-ligand and sp2/aromatic plane checks. The source graph
  remains unchanged; RDKit version `2025.09.6` is retained in evidence.
- `plan_compiler.compile_chain_plan` chooses ordered chain references from graph
  adjacency, source ideal nondegeneracy and the observed mask, not from deleted
  experimental coordinates or conformer scores. It does not change frames to repair
  invalid observations and does not implement a new cone/planar/closure solver.

All **28 predeclared entries** executed: seven components × heavy/H × two proper
rigid transforms. **24 passed both raw and rounded declared-geometry audits**;
four expected compiler refusals remain: EOH heavy-O deletion has fewer than three
observed anchors, and ACT's terminal O has no supported chain frame. All 14 H
entries (10 chiral, four achiral) passed independent raw/rounded CIP readback. Every observation inventory/
coordinate tuple remained unchanged. Request/plan JSON hashes were independently
recomputed; raw proposals and failures remain in the project, not system `/tmp`.

These are **same-CCD-authority ideal roundtrips**, not experimental holdouts,
independent conformational predictions or all20/protein-polymer/H qualification.
PHE's ring passes only for the declared ideal frames; no distorted ring closure
claim follows. The H arm uses the CCD heavy observations directly, not heavy-arm
generated atoms silently promoted to observations. No PDB bytes, `.dat` or public
bundle is emitted. A caller-supplied ideal graph remains a geometry authority,
not a calibrated context-dependent protein geometry library.

After the first formal run, nonnumerical hardening added strict invalid-pin refusal
and refusal of RDKit source-charge cleanup. A labelled deterministic regression
replay produced **all 29 JSON files byte-identical** to the retained first run;
it is not additional independent evidence or another sampling cohort. Original
reports are never overwritten. Tests use temporary mutated copies for negative
source/charge/stereo controls, not edits to tracked inputs or wider thresholds.

Next step: explicit polymer/variant/link/H state authority and prospective observed-
complex controls, then planned planar/closure and conformer/orientation expansion.
Do not simply drop every CCD leaving atom to make a peptide, or apply free-amino-
acid H counts to an internal residue. Reuse existing variant/state/identity owners;
complete heavy and H geometry together without weakening source preservation.

[WP2A validation](chemical-reconstruction-ccd-validation.json): 29 new regression
tests; 495 focused passes and 3429 full-fast passes, three existing optional Torch/
Biotite skips and 33 slow deselected. Strict mypy passes for all 11 private modules;
Ruff passes for all source and the new script/tests. The first full run retained
3427 passes/one fixture-inventory failure/three skips: the inventory recognized only
PDB/FASTA. It now explicitly includes CIF plus the CCD source-manifest allowlist,
retaining exact set equality, duplicate detection and all 63 checksums. Original
upstream CIF padding is preserved through a narrow `.gitattributes` rule, not by
changing source bytes or pins. Existing scientific tests/thresholds, production
paths, environment and prior failure/license gates are unchanged. WP0's inherited
41 legacy mypy errors remain documented, not suppressed or relabelled as passes.

## 12. WP2B conditional peptide chemistry and observed-source probes

[Protocol](chemical-reconstruction-peptide-v1.md) was frozen as `3b19d8d` before
projection or probe execution. It is not edited after outcomes. Existing source
CCDs and the pinned Amber19SB RTP parser/inventory are reused; partial charges
are never converted into formal charges. Canonical naming maps come from the
existing domain owner, not another variant/alias table.

`peptide_state.assemble_witness` applies exact declared C/N substitutions and atom
removals, checks leaving-branch ownership, consumed ports and connectedness, and
proves the joined valence/CIP with explicit chemical extra bonds. Cross-unit
geometry/ideal bond lengths are not invented. Source frames remain independent.
`derive_local_peptide` uses two virtual neutral ALA partners; they are **chemical
witnesses, not actual observed neighbours**, and never become request observations,
generated atoms or output. Parent graphs, removals, source/state inventory and the
complete proof are retained. The four accepted projected inventories have a peptide
N–H count, but its geometry is still unqualified: the H operation and target measured
H refuse explicitly. There is no pH-based state choice or silent LYS→LYN switch.

`observed_peptide` reuses the existing whole-PDB atom/link parsers and peptide
continuity owner. It retains every source chain, heterogen and water, original
source SHA, exact identity and the explicit withhold mask. Only retained observation
DTOs—not withheld experimental coordinates—reach the constructor. Target CB uses
the existing two-angle cone against CA/N/C via a cyclically reordered, CIP-proven
stereo frame; the remaining dependency plan reuses the chain compiler. Carbonyl
phase cannot silently determine CB handedness. Old chain-only defaults and the
frozen CCD roundtrip profile remain unchanged. No new planar/closure solver exists.

### Actual outcomes, including incorrect pre-run expectations

All **12 declared probes** completed. Chemistry: ALA, THR, **LYS** and PHE passed
the conditional joined proof and inventory; ILE and LYN refused.

- The protocol's LYS/LYS refusal prediction was wrong. The pinned CCD LYS is
  already **+1 charged**, not neutral, and its retained inventory matches RTP LYS.
  No charge/H correction or state substitution was made. This is an unexpected
  valid proof, not justification to rewrite the frozen hypothesis.
- ILE unexpectedly failed the fixed canonical inverse-name policy: the comparison
  lacks CCD `CD1/HD11/HD12/HD13` and has RTP `CD/HD1/HD2/HD3`. This is a naming/
  inventory boundary refusal, not evidence of impossible ILE chemistry. No
  opportunistic alias fallback was added. A future qualified naming policy must
  preserve this result and declare its source/name mapping explicitly.
- LYS/LYN refused unsupported explicit variant projection as declared. The +1 CCD
  is not silently deprotonated to satisfy a neutral template.

Observed heavy probes produced **zero coordinate proposals**: 1UBQ A/PHE45 retained
all context and refused partial occupancy (660 source atoms, 653 retained after
seven explicit withholds); 1FR2's selected A/THR11 refused dangling/empty CONECT;
1A3N's selected A/ALA5 refused ambiguous endpoints in the existing declared-link
parser (eight source altloc atoms, four LINK/SSBOND records). These are input/
admission refusals, not failed or passing native geometry trials. All three H arms
gave the intended amide-H-geometry refusal. No context was cropped, heme stripped,
source repaired or tolerance widened. The source files remain byte-identical.

Original evidence is under `tmp/chemical-reconstruction-peptide-v1/probes-01/`.
A separately labelled diagnostic replay adds parent/inventory and admission
payloads; its summary, dispositions and authority digests are identical, and it
is not independent evidence. Synthetic transport tests—not observed structures—
exercise the new CB cone, three-atom THR completion, context retention and source
guards. Those must not be reported as successful observed-complex reconstruction.

### Typing and remaining work

The [WP2B validation record](chemical-reconstruction-peptide-validation.json)
binds the frozen protocol, final source hashes and retained reports: **513 focused
passes; 3448 fast passes, three optional Torch/Biotite skips, 33 slow deselected**.
All 19 new regressions pass. Ruff passes; no observed coordinates were generated.

The private package's strict scoped command is now
`mypy --strict --follow-imports=silent src/dvbfixer/chemical_reconstruction`:
13 modules pass without new ignores/config suppression. The unsilenced strict
command reports **71 errors in four unchanged legacy files** reached through the
reused PDB/continuity adapters. The exact 71 error lines reproduce on the unmodified
`21a7162` source archive with the old sequence-evidence entry point; logs remain
retained. This is separate from WP0's 41-error broader non-strict baseline command,
not a claim that repository-wide typing passes.

Next: qualified explicit residue-name/H alias transport, independently specified
source admission/preparation (preserving originals/context), and amide-H geometry
with external anchors. Do not rerun/repair these frozen observed sources until a
separate prospective policy is declared. Actual neighbour chemistry, full polymer
state/CIP, H geometry, distorted closure and real-complex reconstruction remain
unqualified; production/default/dependency changes remain disabled.

## 13. WP2C explicit inventory naming transport

[Protocol](chemical-reconstruction-naming-v1.md) was frozen as `340d49b` before
the new six-control execution. It is a separate prospective naming policy, not
a repair or relabelling of WP2B's ILE refusal. The old protocol, runner, tests,
CCDs, RTP and canonical naming owner are byte-identical to `7802240`.

The private `inventory_naming.transport_inventory` reuses the domain
`atom_name_policy` with an explicit `AMBER19SB` / `STANDARD` selection. No alias
table was added or edited. It records the raw inventory, complete simultaneous
mapping/rule IDs, normalized inventory, exact source pin and policy selection.
It refuses already-aliased input, implicit-identity collisions, malformed or
unknown local tokens, mutable external rows and unsupported dialect/profile.
The transport DTO verifies that its stored inventories/map agree with the owner
policy, including the source-pin binding. Naming does not choose a chemical state.

External row order and repetitions survive; `-C`/`+N` signs are preserved. Only
external backbone C/N names are qualified here: neighbour residue chemistry is
unknown, so the target residue's sidechain naming rules must NEVER be applied to
an external sidechain atom. This is an explicit refusal, not a reason to discard
the reference. No source PDB, observation identity or generated coordinate is
renamed by this inventory adapter.

### Actual outcomes and limits

All six controls executed. ALA, THR, **ILE**, LYS and PHE passed the transported
exact atom/bond inventory plus existing conditional peptide witness proof. ILE
uses the existing canonical owner for `CD→CD1`, `HD1/2/3→HD11/12/13` and
`HG11→HG13`; no graph search or post-result fallback was used. Renamed atom
counts were 0/0/5/4/1 respectively. LYS/LYN still refused
`explicit-variant-projection-unqualified`; no +1→neutral charge/H edit occurred.
This is not a new formal-charge, bond-order or stereospecific H-coordinate
authority. Equivalent-H naming is deterministic convention, not an independent
pro-R/pro-S geometry qualification.

Original evidence is in `tmp/chemical-reconstruction-naming-v1/controls-01/`.
A labelled regression replay after name-shape hardening preserved all seven JSON
files byte-for-byte and is not independent evidence. No observed WP2B input was
repaired or rerun as a scientific trial. No coordinate constructor ran.

The 36 new regression cases and focused integrations pass (**549 passed**).
The final fast suite passed **3484 tests**, with three optional Torch/Biotite
skips and 33 slow tests deselected; no failures or errors.
Ruff passes. Scoped strict mypy passes 14 private modules; the unsilenced command
still reports the same 71 legacy error lines as the unchanged `21a7162` archive.
No error was hidden with a type ignore or config change.
The [WP2C validation record](chemical-reconstruction-naming-validation.json)
binds the pinned inputs, exact outcomes, final code and retained reports.

Next independent work items remain prospectively specified source admission/
preparation and external-anchor amide-H geometry. Naming transport is now a
reusable boundary, but actual neighbour chemistry, full polymer state/CIP,
observed reconstruction, H geometry and general complex qualification remain
open. Public/default/dependency and data-license gates are unchanged.

## 14. WP2D explicit geometry metrics and external-anchor amide H

[Protocol](chemical-reconstruction-amide-v1.md) was frozen as `8fa8cd9` before
new construction/control execution. Planning had inspected the three pinned
ff19SB parameter rows, explicitly disclosed in that protocol. No new source was
downloaded and no PyRosetta data accessed.

`metric_geometry.MetricGeometryPlan` is a separate private plan type, with explicit
provenance-bound bond/angle metrics involving generated atoms only. It does not
replace the original template geometry or modify CCD chemistry/coordinates.
Duplicate/reversed metrics, observed-only overrides, missing connectivity and
incomplete planar authority refuse. Its new digest binds metric values and
provenance; old `GeometryPlan` fields/digests stay unchanged. The first metric
placement is a planar exterior equal-stiffness angular compromise from three
observed anchors. One anchor may be across an explicit covalent connection;
unknown/absent linkage authority is not guessed from proximity. Other metric
placement families are unqualified, not silently routed to template defaults.

`amide_authority` selects three unique exact rows from the pinned existing
`ffbonded.itp`: H-N length 0.101 nm, C-N-H 120 degrees and XC-N-H 118.04 degrees.
The two angular stiffnesses are equal (418.4). Full-file and exact selected-token
digests bind the source and units; wrong type/function, changed selection, absent/
duplicate rows, unequal stiffness or pin mismatch refuse. The RTP must explicitly
identify N/H/CA as N/H/XC. This is a narrow parameter adapter, not general wildcard
force-field typing, protonation or a calibrated experimental geometry library.

The constructor uses the actual observed neighbour plane and the declared
analytical angular compromise. The measurement-only auditor independently checks
every explicit metric, including cross-component angles, plus retained template
metrics, stereo/planarity constraints and context clashes. Explicit N-H/CA-N-H
metrics replace the incompatible free-CCD H geometry; no inter-template external
ideal frame or bond length is manufactured. Raw failures and generated-only
PDB-rounded audits survive. No observation is moved, optimizer called or second
engine selected. No public API/CLI/GUI/default backend changed.

### Actual outcomes and evidence boundaries

Six synthetic controls executed: ALA, THR, ILE, LYS and PHE passed both raw and
rounded geometry; ALA with incoming angle 100 degrees was rejected in both.
The positive anchor angle 121.96 degrees is intentionally the exact complement
of the two parameter equilibrium angles: these are **same-parameter synthetic
consistency controls**, not independent experimental validation. Positive raw
metric residuals are numerical zero; largest rounded length/angle deviations are
about 0.00061 A / 0.0168 degrees. The distorted case retains 10.98-degree raw
residuals on both angles against the fixed 5-degree tolerance. No source repair,
tolerance change, anchor movement or coordinate ranking was used.

Original evidence is in `tmp/chemical-reconstruction-amide-v1/controls-01/`.
Post-control hardening binds exact parameter-row tokens, rejects a plain plan
trying to use metric placement without metric authority, and requires all three
anchors to be actual observations rather than merely referenced link metadata. Its labelled regression
replay preserved all seven JSON records byte-for-byte. Additional old CCD,
peptide and naming replays preserved all 29/13/7 prior JSON files byte-for-byte,
including old defaults/digests, WP2B ILE refusal and all H refusals. Replays are
regression checks, not independent scientific evidence.

All **31 new regression cases** pass; final focused integrations passed **580 tests**.
The earlier 579-case focused and 3514-pass fast runs are retained but superseded
by the final observed-anchor gate and its extra test; they are not labelled as
final-source verification.
The final-source fast suite passed **3515 tests**, with three optional
Torch/Biotite skips and 33 slow tests deselected; no failures or errors.
They include proper rigid-frame equivariance, a deliberately poisoned free-CCD H
coordinate (no effect on output), audit corruption, external-link/metric/plane
gates and a retained opaque-context clash. One initial test fixture incorrectly
attempted a coordination record with covalent order and was rejected by the
existing contract before the intended assertion; the fixture now uses the wrong
covalent endpoint. The shared connection guard was not weakened. Ruff passes;
scoped strict mypy passes 17 private modules, while the unsilenced strict run
retains the exact same 71 legacy diagnostic lines as baseline `21a7162`.
The [WP2D validation record](chemical-reconstruction-amide-validation.json)
binds the source/protocol/code/report hashes, exact outcomes and qualification
limits; CLI reference and agent map checks pass with explicit source-path isolation.

No observed PDB source was repaired/rerun, and the existing observed adapter still
refuses H operations. This is a reusable metric boundary plus a first planar-H
primitive, **not** full H placement or a one-atom product definition. Measured H
transport, stereospecific methylene/methyl H geometry, arbitrary missing subsets,
actual neighbour chemistry, full-polymer state and real retained-complex H
qualification remain open. Source admission/preparation is the next separately
declared slice; preserve original full context and the frozen WP2B refusals.

## 15. WP2E source-bound observation preparation

[Protocol](chemical-reconstruction-source-v1.md) was frozen as `0d83810` before
the three new source controls. It explicitly declares metadata preparation and
conservative opaque fractional occupancy, rather than repairing/rerunning WP2B
under its old policy. Original fixtures, FF files and old protocols stay exact.

`source_preparation.SourcePreparation` retains original/prepared ASCII/LF text,
both hashes and exact one-based removed-record evidence. Its constructor recomputes
the projection: only whole dangling/empty CONECT records can disappear. Valid
records, ATOM/HETATM/ANISOU, coordinates, occupancy, names/order, headers and
linkage metadata otherwise remain byte-identical. Malformed/self CONECT, altloc
or multiple-model ambiguity refuse, rather than being sanitized opportunistically.
The library accepts bytes and never overwrites a source path.

`PreparedContextRequest` is a separate digest-tagged private contract; the original
request fields, digests and occupancy refusal remain unchanged. Its explicit policy
binds prepared-source hash, preparation digest, retained-observation digest and
the exact fractional set. Only unselected unresolved UNKNOWN components qualify;
fractional targets, resolved chemistry, owned/mobile atoms, constraint atoms and
connection endpoints refuse. Target occupancy is checked BEFORE withholding.
Fractional context coordinates/occupancies remain unchanged and fully participate
in conservative clash checks: they are not normalized, moved or probabilistically
discounted. No altloc exemption or observed-H capability was added.

`audit_source_transport` independently reparses the prepared source and checks the
exact ordered retained observation DTOs, counts, mask scope and origin/evidence
bindings before compilation/construction. A rehashed but source-inconsistent DTO
still fails this audit. The full original/prepared payloads and withheld reference
coordinates do not enter the coordinate constructor. Existing admission, graph
compiler, canonical naming, peptide proof, proposal and raw/rounded audits are reused.

### Actual real-source outcomes, not broad qualification

All **three declared controls** executed; preparation completed for two, requests
and constructors executed for two, and both raw/rounded geometry audits passed:

- **1UBQ A/PHE45:** seven sidechain heavy atoms proposed. Original and prepared
  source bytes are identical; all 660 source atoms remain except the seven explicit
  in-memory withholds (653 observations). All **57** fractional opaque-context
  observations retain their original occupancy and coordinates; no metadata changed.
- **1FR2 A/THR11:** three sidechain heavy atoms proposed, 1686 observations from
  1689 source atoms. Exactly **nine dangling CONECT records** were dropped as whole
  records; their original bytes and line numbers remain in evidence. Every structural
  record and every other metadata line is unchanged. No chemistry/bond was invented.
- **1A3N:** original full source retained, refused
  `explicit-altloc-selection-required`. No prepared file, request or constructor
  exists for that entry. No occupancy winner, alternate-atom deletion or heme strip.

These are explicit observed-sidechain masks and fixed predeclared geometry checks,
**not native-conformer accuracy**, template superiority, full-source chemical
assessment, parameterization or MD readiness. No held-out RMSD scoring/ranking,
geometry-bound adjustment, alternative engine or coordinate repair was performed.
No PDB candidate, `.dat` or public bundle was written; retained `prepared.pdb`
files are source metadata projections only. All source-coverage/MD/public flags
remain false, and unresolved retained context remains visibly unassessed.

Original evidence: `tmp/chemical-reconstruction-source-v1/controls-01/`.
Post-control hardening adds observation-digest binding and independent transport
audit. Its labelled `binding-replay-02/` has an identical summary, all five retained
source PDBs and all raw/rounded additions, but request/plan digests intentionally
change with the stronger binding; **the full evidence JSON is not byte-identical**.
Independent saved JSON/source readback re-audited both candidates in raw/rounded
form without calling the proposal constructor, and verified exact source changes,
observation retention and every request/plan/source digest. Replays/readback are
regressions, not additional independent scientific controls.

### Shared helper extraction and checks

The first adapter imported the existing `renumber.update_conect` owner. Static
mypy traversal then reached 795 legacy errors in 51 CLI-dependent files; all 795
error lines exactly reproduced on unchanged `21a7162` with the same legacy entry
points. This was not a private-package typing failure, but unnecessary coupling.
The existing helper was moved, **not reimplemented**, into `pdbutils.io`, with an
identical executable-body/docstring AST and a typed signature. `renumber.update_conect`
remains the same callable via import alias; `pdbutils` re-exports it too. A clean
subprocess verifies source preparation does not import the renumber command module.

Final scoped strict mypy passes **19 files** (18 private modules plus the extracted
pure I/O owner). Unsilenced checking is back to the same 71 diagnostic lines in
four untouched legacy files; both the initial 795-error and final 71-error logs
are retained, not suppressed or labelled passes. Ruff passes. All **23 new
regression cases** pass; final focused integrations passed **610 tests**, including
renumber/SEQRES controls. The earlier 609-case focused and 3537-pass fast run are
retained but superseded by extraction plus its import/alias regression.
The final-source fast suite passed **3538 tests**, with three optional Torch/Biotite
skips and 33 slow tests deselected; no failures or errors.
Old CCD/peptide/naming/amide replays preserved all 29/13/7/7 JSON files byte-for-byte,
including WP2B input and H refusals and original request/plan digests.
The [WP2E validation record](chemical-reconstruction-source-validation.json) binds
protocol/source/code/report hashes, observed outcomes, source readback and limits.

Next: prospective broader retained-complex and missing-subset evidence plus explicit
H inventory/geometry transport on actual sources. Full-polymer/actual-neighbour
chemistry, measured-H support, complete H families, distorted ring closure,
qualification and public/default/dependency/data-license gates remain open. Do not
promote two geometry successes to a public reconstruction or whole-complex claim.

## 16. WP2F complete selected-residue H with fixed heavy observations

[Protocol](chemical-reconstruction-hydrogens-v1.md) was frozen as `3b85d5c`
before the twelve entries executed. Its prose incorrectly says five 1UBQ targets;
the explicit table has **four** and twelve total entries. That typo was disclosed
before execution; the frozen file is unchanged, the table followed exactly and
no extra target substituted. Four targets in one source are not four independent folds.

`hydrogen_completion` composes existing constructors for a complete explicit H
inventory, never just the convenient atoms. Heavy observations, all context and
occupancies remain fixed. It reuses source preparation/audit, canonical naming,
conditional peptide chemistry and pinned backbone N-H metrics. Carbon CH/CH2
groups use ideal-name signed cones, terminal C-H3/+1 N-H3 and hydroxyl groups use
one deterministic heavy chain frame, and aromatic CH uses explicit CCD scalar
planarity metrics. All frames are observed-heavy and selected from topology/ideal
authority only; no observed-success frame retry, phase search, pH inference,
rotamer optimization, reflection or heavy-atom repair exists. Unlisted groups,
missing heavy atoms, measured target H/D and partial target occupancy refuse.

`GeneratedMetricGeometryPlan` independently binds stricter generated-atom bond/
angle bounds (0.02 A / 5 degrees) while retaining fixed-heavy bounds (0.10 A /
12 degrees). Both use the same 0.7 A clash limit. H planes use 0.005 A; heavy
planes keep 0.05 A. Old plan types/digests and default authority constraints are
unchanged. Both fixed-heavy stereo frames and full four-ligand frames survive.
RDKit also checks the complete joined witness with actual target heavy + generated
H coordinates. This is conditional target CIP with virtual ALA chemistry, NOT
actual-neighbour or full-polymer chemistry, H-position accuracy or pro-R/pro-S
qualification. The old observed adapter's H refusal remains unchanged; this is a
separate private research route, not CLI/Homology/GUI/public capability.

### Actual outcomes: success and negative evidence remain separate

Twelve declared entries completed. Five nominal synthetic complete-H views
(ALA/THR/ILE/LYS/PHE: 5/7/11/13/9 H) passed raw/rounded geometry and conditional
CIP. The synthetic ALA incoming-angle-100 case generated all five H, retained its
amide-angle geometry rejection, and still passed conditional CIP.

Six observed-source entries, selected before geometry inspection:

| Target | Outcome |
|---|---|
| 1UBQ A/ALA28 | All 5 H, raw/rounded geometry and conditional CIP passed |
| 1UBQ A/ILE3 | All 11 H retained; raw/rounded angle rejection, conditional CIP passed |
| 1UBQ A/LYS6 | All 13 H retained; raw/rounded angle rejection, conditional CIP passed |
| 1UBQ A/PHE45 | Fixed heavy-ring plane failed admission; no H constructor/CIP executed |
| 1FR2 A/THR11 | All 7 H, raw/rounded geometry and conditional CIP passed |
| 1A3N | Original full source retained; explicit altloc-selection refusal before request |

There are **two observed geometry/CIP successes**, two geometry rejections, one
fixed-input refusal and one preparation refusal, not six successful H trials.
Ten total candidates (six synthetic/four observed) were constructed; all ten
passed conditional CIP in raw and rounded forms, but only seven passed geometry.
Nothing was filtered, relabelled, cropped, repaired or loosened to increase that count.

Post-hoc diagnostic measurements identify ILE HA-CA-N and CG2-CB-HB residuals
of about -5.75/-6.72 degrees, and LYS HB2-CB-HB3 / HD2-CD-HD3 residuals of about
-9.73/-5.61 degrees, against the unchanged 5-degree bound. Independent exact-cone
H placement is not adequate for all observed tetrahedral distortions. PHE45's
already-observed ring deviates 0.06216 A from its declared plane versus 0.05 A;
this is not a failed H construction. The earlier WP2E PHE mask rebuilt the ring
instead of treating all its deposited coordinates as fixed, so that success is
not contradicted. No new repair policy is inferred from either result.

An initial unit hypothesis that moving HA alone would force RDKit CIP failure
was wrong: the test retained S with unchanged heavy orientation while geometry
correctly rejected the H. That case is now explicit negative geometry evidence;
a mirrored complete synthetic stereocentre is separately rejected by CIP.
Swapping identical methylene H can likewise retain CIP while failing ordered
name-frame geometry. Do not use CIP alone as H-quality or stereospecific-name proof.
The audit was not weakened or changed to produce an expected result.

Original evidence is under `tmp/chemical-reconstruction-hydrogens-v1/controls-01/`.
All 24 retained files replay byte-for-byte on final source. Saved source/request/
plan/candidate readback independently remeasured 20 raw/rounded geometry audits
and 20 conditional CIP audits **without proposal construction**, retaining every
negative outcome. No experimental H reference, native-H score, phase ranking or
candidate PDB/.dat/public bundle exists. Source/FF bytes and prior evidence remain
unchanged; old CCD/peptide/naming/amide/source replays preserve 29/13/7/7/9 files.

### Verification and unresolved integration behavior

All **27 new regression cases** pass. Final focused fast lane: **634 passed,
three existing slow tests deselected**. The 27 new cases plus both seeded rebuild
slow boundaries passed together (**29 passed**, 14.6 s); renumber's slow boundary
passed separately (one pass). The full fast suite passed **3565 tests**, three
optional Torch/Biotite skips and 33 slow deselected. Ruff passes; scoped strict
mypy passes 20 files (19 private modules plus I/O owner). Unsilenced mypy retains
the same 71 diagnostic lines in four unchanged legacy files, not suppressed.

The large 637-case slow-inclusive focused run was NOT completed after the unit
correction: outer limits of 120 and 600 seconds expired while the seeded prepare
test was displayed. Logs are retained. That legacy test separately passed in
11 seconds and in the 29-case boundary run; the cause of the combined-run stall
is unresolved. Do not relabel either timeout as a pass, attribute causality to
the new H code without evidence, add skips, weaken assertions or change the
production numerical platform. This remains a separate integration investigation.
The [WP2F validation record](chemical-reconstruction-hydrogens-validation.json)
binds protocol/source/report hashes, all positive/negative outcomes, readback and
the unresolved combined-run timeout evidence.

Next at WP2F (completed in WP2G below): declare a joint generated-only sp3-group fitting policy to distribute angular
strain instead of independently satisfying two cone angles, preserving current
negative results and the same observation ownership. Fixed-heavy ring defects
require a distinct source policy, not a hidden H workaround. Broader complex/
chemistry evidence, measured/partial H transport, all20 and stereospecific H names,
actual neighbour/polymer chemistry, hydrogen-bond quality and production/default/
dependency/data-license gates remain open.

## WP2G — joint generated-only sp3 H orientation (2026-11-20)

The [frozen joint-H protocol](chemical-reconstruction-joint-h-v1.md) was committed
as `7d8267b` on the persistent branch before numerical execution, with baseline
`8daa3b6`. No branch/worktree, public capability, default replacement, chemistry
change, mask/phase search or relaxed threshold was introduced.

`joint_hydrogen_geometry.py` adds explicit `DirectionGroup` and digest-tagged
`JointMetricGeometryPlan` records. `compile_joint_hydrogens` upgrades only neutral
carbon 3-heavy/1-H and 2-heavy/2-H cone groups. All heavy-neighbour directions
participate equally in one proper rotation, reusing the unchanged
`model.diffusion.geometry.weighted_kabsch` owner on normalized antipodal clouds.
Translation is checked but never applied; observations stay fixed. Each group's
generated H shares that rotation and keeps source lengths/internal H-H angles.
Backbone amide, aromatic CH, hydroxyl and terminal-three-H strategies are unchanged.
The constructor caches one fit per group; the independent candidate auditor checks
authority, complete membership, fixed-heavy ownership and geometry without fitting.

All twelve declared comparison entries executed with the EXACT WP2F requests,
constraints, states, masks, names and fixed/generated audit policies:

- Five nominal synthetic controls passed raw/rounded geometry and conditional CIP.
- Observed 1UBQ ALA28 (5 H), ILE3 (11 H), LYS6 (13 H) and 1FR2 THR11 (7 H) passed
  raw/rounded geometry plus conditional CIP. ILE/LYS improve over the unchanged,
  retained WP2F angular failures; the old algorithm still rejects those inputs.
- Synthetic ALA incoming angle 100 remains geometry-rejected despite CIP pass.
  PHE45 still refuses its fixed-heavy ring plane before construction; 1A3N still
  refuses unresolved altloc. No fixed-heavy source workaround was added.

These are known retrospective cases with a prospectively frozen algorithm, not
independent folds or native-H accuracy evidence. Of six observed declarations,
five produce prepared requests, four run coordinate construction and four pass;
the other two are fixed-input/preparation refusals. Conditional CIP remains
necessary but does not validate H coordinates or stereospecific H naming.

Saved old/new views were decoded independently: request/constraints/policies and
plan bindings match; saved observed originals/prepared bytes were rechecked.
Twenty raw/rounded geometry and twenty conditional-CIP audits match the retained
reports with constructor/fitter calls forbidden. All 113 files across joint,
complete-H, source, amide, naming, peptide and CCD replays are byte-identical.
Readback and replays are regressions, not additional scientific trials.

Verification: 34 new cases; focused fast **668 passed**, three slow deselected;
separate joint/H/seeded/renumber boundary process **67 passed**; full fast
**3599 passed**, three existing optional Torch/Biotite skips, 33 slow deselected.
Ruff passes and scoped strict mypy passes 21 files. Unsilenced mypy exposes 793
legacy diagnostics in 51 files through the reused helper's dependency graph:
the exact diagnostic multiset matches an unchanged `8daa3b6` source archive,
with zero new errors. The helper was not copied/refactored to hide those errors.
CLI/GUI generated checks pass when bound to this checkout (`PYTHONPATH=src:.`).
The old WP2F large combined slow-inclusive timeout is still unresolved, not
closed by this small boundary run. See the
[WP2G validation record](chemical-reconstruction-joint-h-validation.json) for
source/protocol/report hashes and retained evidence paths.

Next at WP2G (bounded extension completed in WP2H below): prospectively declare broader explicit-chemistry/H coverage and independent
quality evidence before extending this algorithm. Keep heavy reconstruction and
H geometry as equal tracks within the shared retained-complex architecture.
Measured/partial H, all20/linked chemistry, actual-neighbour/polymer witnesses,
stereospecific names, native-H/hydrogen-bond quality and broad complex qualification
remain open. Fixed-heavy ring repair requires a separate explicit source policy.
Public/default/dependency replacement and production/data-license approval remain
disabled; local PyRosetta research permission is not production permission.

## WP2H — explicit planar-N and thiol H groups (2026-11-20)

[Protocol](chemical-reconstruction-extended-h-v1.md) frozen as `c6c13b0` before
acquiring eight unmodified RCSB CCD sources (ASN, GLN, ARG, CYS, TRP, GLY, PRO,
HIS). The separate CC0 fixture manifest was pinned in `c3f8c0c` before numerical
controls. The original CCD manifest, force fields, naming owner and all previous
protocols/pins remain unchanged. This is still the same persistent branch.

The separate `compile_extended_hydrogens` entry point admits N-H/two-heavy and
N-H2/one-heavy groups only when their exact parent/neighbour set has the existing
CCD authority's sanitized SP2 planarity proof and supported charge/H bonds. It
also admits neutral single-bond thiol S-H, not thiolate or disulfide chemistry.
Planar N-H reuses metric planar-exterior placement; named N-H2 and S-H use pinned
CCD chain torsions in deterministic observed-heavy frames. N-H2 has an additional
0.005-A complete-group/heavy-anchor plane. The unchanged joint-carbon compiler is
then applied. No numerical engine, phase search, protonation inference, measured-H
replacement, source repair or observed movement was introduced. The explicit origin
policy is `complete-H-groups-v2-explicit-sp2-N-thiol`; old defaults still refuse
the new group chemistries.

All sixteen declared entries executed:

- Six synthetic ASN/GLN/ARG/CYS/TRP/GLY cases passed raw/rounded geometry and
  conditional CIP (6/8/13/5/10/3 H). These are ideal-authority consistency controls,
  not native-H holdout evidence.
- Observed 1UBQ A/ASN25 (6 H), A/GLN2 (8 H), A/GLY10 (3 H), and 1FR2 A/CYS23
  (5 H), A/TRP74 (10 H) passed raw/rounded geometry plus conditional CIP with every
  observation/context record retained under the existing source policy.
- Observed 1UBQ A/ARG42 was refused before coordinate construction: its fixed-heavy
  guanidinium plane deviation is 0.099886509 A against the unchanged 0.05-A bound.
  No target substitution or implicit heavy-geometry correction was attempted.
- Synthetic and observed PRO refuse the existing H2-based peptide substitution
  (`removal outside named source components`). Both HIS entries refuse because
  the exact same-named RTP state is unavailable; no HID/HIE/HIP substitution.

Thus eight observed declarations yield six prepared requests, five constructors,
five raw/rounded geometry/CIP passes, one fixed-input and two chemical-authority
refusals. No all20 or missing-heavy-sidechain support follows from this H slice.
Chemistry still uses virtual neutral-ALA witnesses, not actual-neighbour polymer
chemistry. Conditional CIP still does not establish H-coordinate/pro-R/pro-S quality.

Independent saved readback decoded twelve views, rechecked all eight original/
prepared observed sources, repeated 22 geometry and 22 conditional-CIP audits with
constructor/fitter calls forbidden, and remeasured the ARG fixed-plane refusal.
All 146 files across the new controls and seven old arms replay byte-identically;
these are regression replays, not independent trials.

Verification: 39 new cases; final focused fast **708 passed**, three slow deselected;
split joint/H/seeded/renumber/fixture boundaries **107 passed**; full fast
**3638 passed**, three existing optional Torch/Biotite skips and 33 slow deselected.
Ruff and scoped strict typing pass. Source-bound, nonincremental strict checking
also passes 22 explicitly selected files including the shared geometry/sequence
dependencies. The ambient-import invocation still reproduces the exact 793 legacy
diagnostics: retain the import-context distinction, not a global typing-pass claim.
An initial test expected a nonexistent specific planarity finding code; it now
checks the real constraint-conflict code AND directly measures failure of that
same plane. The fixture test initially lacked the new metadata allowlist entry;
the explicit source manifest was added, without changing checksums/assertions.
Initial failures remain in the [validation record](chemical-reconstruction-extended-h-validation.json).
WP2F's large combined slow timeout remains unresolved.

Next at WP2H (bounded PRO/HIP declarations completed in WP2I below): declare explicit leaving-atom/variant state proofs (not a residue-name guess)
before PRO or histidine can progress; continue toward all20 heavy/H coverage and
full retained-complex sidechain reconstruction. Native-H/name/orientation evidence,
measured/partial H, actual-neighbour/polymer chemistry, linked/ring heavy repair,
broader complex qualification and production/public/default/dependency/license
gates remain separate and open. Do not repair ARG's fixed source or infer a state
merely to increase this protocol's pass count.

## WP2I — declared PRO/HIP states and source-only metrics (2026-11-20)

The [state declaration protocol](chemical-reconstruction-declared-state-v1.md)
was frozen as `19dc275` from `1d52019`. The new private
`PeptideStateDeclaration` binds exact target/state, parent/ALA/inventory digests and
both directed substitutions. `DeclaredPeptideAuthority` subclasses the existing
authority but re-proves its whole projection on construction and uses a separate
tagged digest/provenance. The shared witness assembler, inventory transport,
source adapter, variants owner and independent validators remain canonical.

PRO explicitly removes its CCD-marked leaving N-H, not absent H2. HIP is an
explicit caller-requested state whose name maps to HIS through the existing
`VARIANT_TO_PARENT`; the unchanged +1 CCD HIS graph must exactly match the HIP RTP
inventory after peptide leaving-atom removal. No sidechain atom, formal charge or
bond order is edited. This does NOT establish the native protonation of 1UBQ HIS68
or select a state from pH. Unchanged-graph HID/HIE declarations refuse their exact
inventory mismatch; neutral transformations remain unimplemented.

The H compiler requires a real fixed/unit-occupancy incoming C with an explicit
order-1 peptide bond. PRO uses no fictitious backbone H/type and adds the fixed
N/Cpre/CA/CD plane at the existing 0.05-A bound. HIP reuses exact ff19SB amide
typing. Virtual ALA chemistry never substitutes for real observed neighbours or
enters output; actual-neighbour/polymer heavy-angle qualification is still open.

Two initial PRO integration refusals are retained, not retrospectively passed:

1. `controls-01`: chemistry was proved, but no backbone H meant no amide override
   metrics and the nonempty metric-plan contract refused before construction.
2. Before repeating, [the integration addendum](chemical-reconstruction-declared-state-metrics-v1.md)
   was frozen as `50fb041`, materializing complete source H metrics rather than
   weakening that contract. `controls-02` then hit the old nonplanar-override guard,
   again without constructing any PRO coordinates.
3. Before another run, [the source-only ledger contract](chemical-reconstruction-source-metrics-v1.md)
   was frozen as `96da6bb`. Separate `SourceMetricGeometryPlan` and
   `JointSourceMetricGeometryPlan` now bind complete H bond/angle coverage to the
   exact canonical graph, provenance, protium mask and unchanged numerical source
   values. The auditor checks that equality without generation. Old planar
   override and nonempty guards remain unchanged; arbitrary cone/joint overrides
   remain refused. PRO carries 25 genuine source metrics, not dummy N-H rows.

Final `controls-03` executes the SAME eight declarations with unchanged states,
observations, masks, synthetic anchors and bounds: synthetic PRO (7 H) and HIP
(8 H), observed 1UBQ A/PRO19 (7 H) and A/HIS68 explicitly requested as HIP (8 H)
all pass raw/rounded geometry and conditional CIP. HID/HIE unchanged-graph proofs,
wrong-PRO-H2 and wrong-source-pin remain refused. No state/target substitution,
phase search, heavy repair, optimizer, reflection or tolerance change was used.
These are bounded controls on known cases, not independent folds/native-H quality.

Saved readback decodes all four final views/declarations/source proofs and repeats
eight geometry and eight conditional-CIP audits with coordinate construction and
fitting forbidden. The two observed original/prepared sources match; HIP and all
four chemistry-negative evidence files match BOTH earlier runs byte-for-byte.
PRO input declarations, chemical graphs, observations, masks and typing match both
earlier runs. Final replays across nine arms reproduce all 159 files byte-for-byte,
including the original WP2H PRO/HIS refusals; no repeat is counted as new science.

Verification: 36 new cases; final focused fast **744 passed**, three slow deselected;
split declared/extended/joint/H/seeded/renumber/fixture boundaries **143 passed**.
Full fast **3674 passed**, three existing optional Torch/Biotite skips, 33 slow
deselected. Ruff passes;
scoped strict typing passes 22 selected files and source-bound nonincremental strict
checking passes 23. Ambient-import typing still exposes the exact 793 prior legacy
diagnostics, not suppressed or called a global pass. An initial wrong-kind test
fixture retained covalent order on an UNKNOWN link and correctly failed DTO
construction; its valid unknown-link fixture now tests the compiler's refusal.
See the [validation record](chemical-reconstruction-declared-state-validation.json)
for all source/report/protocol hashes and initial/final outcomes.

Next at WP2I (bounded coverage completed in WP2J below): retain the separation between explicit state selection and graph-changing
microstate transformations. Continue all20 heavy/H and whole-sidechain/closure
work with declared chemistry, not inferred neutral histidines or a PRO-only
product. Measured/partial H, actual-neighbour/polymer evidence, external heavy
angles, native-H/name/orientation quality, broader linked/complex qualification,
the unresolved WP2F combined slow timeout and production/public/default/dependency/
data-license gates remain open.

## WP2J — all20 source coverage, separate heavy/H arms (2026-11-20)

The [all20 protocol](chemical-reconstruction-all20-v1.md) was frozen as `cbb28b4`
from `91ef103`, before acquiring the seven remaining unmodified CC0 CCD sources
(ASP/GLU/LEU/MET/SER/TYR/VAL). Their separate manifest and global fixture inventory/
metadata allowlist were committed as `6c58855` before numerical controls. The
original sources/manifests, FF, naming and variant owners remain unchanged.

The exact declared state table covers twenty canonical SOURCE residue types;
ASP requests ASH, GLU requests GLH, HIS requests HIP. These are explicit experimental
states proved against unchanged CCD graphs, NOT native-state/pH inference or all
default-state coverage. Every row uses the existing declared-state witness,
inventory transport and validators. No new coordinate strategy, ring solver,
rotamer search, phase ranking, source repair or public capability was introduced.

`compile_peptide_request_plan` extracts the old probe helper's existing CB-cone /
chain behavior, avoiding fake source-probe metadata for synthetic requests.
`heldout_scoring.py::score_heavy_holdout` adds exact-identity, complete-mask,
unit-occupancy heavy scoring AFTER proposals: unchanged-frame RMS/max displacement,
no alignment, symmetry relabelling, acceptance gate or constructor input. Reference
coordinates are retained separately; neither request nor compiler accepts them.
Complete geometry-rejected candidates can be scored; incomplete ones cannot.

All **80 declared entries** executed (20 sources x synthetic/observed x heavy/H):

| Arm | Declared | Raw/rounded geometry passes | Other outcomes |
|---|---:|---:|---|
| Synthetic whole-sidechain heavy | 20 | 19 | GLY has no sidechain-heavy atoms: not applicable |
| Synthetic complete H | 20 | 20, also conditional CIP | Same-authority consistency, not native quality |
| Observed whole-sidechain heavy | 20 | 18 | GLY not applicable; no eligible interior MET |
| Observed complete H | 20 | 16, also conditional CIP | Three fixed-plane refusals; no eligible interior MET |

Heavy withholds every canonical sidechain heavy atom, not just one/two atoms;
H keeps all original selected heavy atoms fixed. **These are separate arms, not a
heavy-then-H workflow.** Heavy uses the unchanged WP2E 0.10-A / 12-degree policy;
H uses fixed-heavy 0.10-A / 12-degree and generated 0.02-A / 5-degree bounds.
Do not pool counts as one identically gated scientific success statistic.

Observed H refuses 1UBQ ARG42, explicitly ASH/ASP21 and TYR59 before generation:
fixed-plane deviations 0.099886509 / 0.141325012 / 0.064521902 A exceed the unchanged
0.05-A bound. Heavy reconstruction deliberately removes those sidechains; its
passes do not repair or overturn the fixed-input H refusals. MET1 is present only
at the source terminus, not the required interior target; no later/different file
was substituted. Both synthetic MET arms pass, so this is a selection limitation,
not failed MET chemistry. GLY empty masks are never counted as reconstruction
passes. PHE here is the declared first interior PHE4, NOT historical PHE45; THR is
1UBQ THR7, NOT historical 1FR2 THR11. Previous refusals still replay unchanged.

**Geometry acceptance is not native conformation recovery.** All eighteen complete
observed heavy candidates have raw/rounded exact-name, unsymmetrized holdout
scores. Raw RMSD median is **1.4313 A**, range **0.0867–5.9047 A**; worst is the
ten-heavy TRP74 sidechain in 1FR2, despite its declared-geometry pass. The engine
uses one pinned ideal phase, not complex-aware native conformer selection. Scores
never influenced generation/ranking/acceptance and do not establish superiority
over a baseline. Native-H accuracy remains unmeasured. These source families and
synthetic consistency repetitions are not independent folds.

Independent saved readback verifies 40 original/prepared source pairs, 78 requests,
76 plans, 146 geometry audits, 72 full-H conditional-CIP audits and 74 heavy scores,
with proposal/fitter calls forbidden. It separately remeasures the three fixed
planes and checks native heldout identities/coordinates against saved originals.
All 320 files across ten final replay arms match byte-for-byte, including old WP2H
refusals and WP2I's final control outputs. Earlier WP2I PRO interface failures remain
retained separately, not reclassified. Readbacks and replays are regression checks,
not new scientific trials.

Verification: **38 new cases** pass; focused fast **782 passed**, three slow
deselected; split H/heavy/seeded/renumber/fixture boundaries **181 passed**. Full
fast **3712 passed**, three existing optional Torch/Biotite skips, 33 slow
deselected. Ruff and scoped
strict typing pass 23 selected files; source-bound nonincremental strict checking
passes 24. The ambient invocation retains the exact 793 legacy diagnostics, not
hidden or counted as a global pass. Source fixture integrity is checked with exact
set/checksum equality. See the
[validation record](chemical-reconstruction-all20-validation.json) for hashes,
denominators, native-score evidence and retained refusal details.

Next at WP2J (bounded nonempirical baseline completed in WP2K below): freeze a genuinely complex-aware conformer proposal/ranking policy and its
heldout comparison before changing geometry or claiming native recovery. Keep
full-sidechain and H correctness as equal tracks; their combined lineage/ownership
and independent final audit remain a separate integration gate. No reduced masks,
relabelled native references, permissive clash/plane thresholds or chosen-success
targets to improve the observed table. Broader chemistry, actual-neighbour/polymer
evidence, measured/partial H, native-H/name/orientation quality, the old unresolved
WP2F slow timeout and public/default/dependency/production/license gates stay open.

## WP2K — graph-fragment grid with retained-context ranking (2026-11-20)

[Protocol](chemical-reconstruction-conformer-grid-v1.md) frozen as `802c677` from
`19025c8`; the path-independent
[baseline signature ledger](chemical-reconstruction-conformer-baselines-v1.json)
was committed as `cc43cca` BEFORE new search. It records only WP2J source/state/
mask/status and rounded GENERATED-coordinate hashes, not native coordinates.
The same twenty observed declarations are used; GLY remains N/A and the eligible
interior MET is still absent. All baseline signatures match before any rotation.
No new source, PyRosetta data, empirical prior, native engine or default change.

`torsion_geometry.py` compiles graph-cut axes from the existing sanitized chemical
witness, not residue-specific chi tables. Eligible single nonaromatic/nonconjugated
SP2/SP3 bridges point away from exact observed roots; all moving downstream heavy
atoms must be generated. Cycles, conjugated axes, observed distal anchors, empty
terminal motions and external generated links have explicit exclusion reasons.
PRO ring puckering is not sampled. Whole branches and aromatic/fused rings remain
rigid fragments. The rotator reuses RDKit Get/SetDihedralDeg; source observations
must remain bit-identical after EVERY nonzero operation, never restored afterward.
Scratch witness/H coordinates never enter output or become qualified H geometry.

`conformer_search.py` has a separate digest-bound search plan/strategy. The ordinary
single-proposal constructor explicitly refuses it. Grid offsets are `(0,-120,+120)`
relative to baseline torsions, complete Cartesian products in root/identity order,
at most 256 combinations (no truncation). Index 0 preserves baseline additions
exactly. All attempts retain raw/rounded geometry, all-path torsion trace and steric
score evidence. The measurement-only trace checks every generated-involving local
heavy four-path, including UNSAMPLED central bonds, at raw 1e-6/rounded 0.25 degrees.
These are provenance tolerances, not changed geometry acceptance bounds.

`steric_ranking.py` reuses the existing diagnostic VDW radii. It counts every unique
eligible pair involving a generated atom against the FULL retained context and
other generated atoms, excludes authoritative 1-2/1-3 paths, never downweights opaque
fractional atoms, and refuses missing radii/over-budget contexts rather than cropping.
After unchanged raw+rounded heavy geometry AND trace gating, selection minimizes
rounded sum of squared overlaps, then worst overlap, then index. This is NOT a force
field, calibrated energy/posterior, H-bond/solvation model or new VDW pass threshold.
No reference/RMSD enters construction or ranking. The driver only extracts native
holdouts after search/selection has returned; tests enforce that boundary.

All twenty entries executed: eighteen selected searches, one GLY N/A, one missing
eligible interior MET. Across the eighteen searches all **302** conformations were
constructed, **265** passed the combined raw/rounded geometry+trace gate, **37**
failed geometry, and **zero** failed the all-path trace. No new geometry, clash,
plane or source-preparation threshold was tuned. This is not an all20 native pass rate.

On the SAME known cases, rounded exact-name heavy RMSD median changes from
**1.431404 to 0.530263 A**: nine improved, eight unchanged and **one worsened**.
1FR2 TRP74 improves 5.904726 -> 0.263223 A. The important retained negative is
explicit ASH/1UBQ ASP21: steric-overlap sum improves 1.000021 -> 0.690350 A^2 while
native exact-name RMSD worsens **0.530115 -> 1.535525 A**. It is NOT rescued with
the reference, a new objective term, atom-name swap or another grid. ASH is a
declared experimental state, not inferred native protonation; unsymmetrized RMSD
and steric overlap are separate measures. WP2J already recorded this source's
fixed-plane defect, so RMSD alone is not chemical correctness; the case stays in
the denominator. Old WP2J bad native scores stay intact.
These retrospective source families are not independent folds or broad/native-H
qualification; a nonempirical relative grid is not a Dunbrack library.

Constructor-free saved readback (including SetDihedral disabled) verifies twenty
original/prepared source pairs, all topologies/grids/rounding/baselines/selection,
36 baseline geometry audits, **604** attempt geometry audits, **604** all-path trace
audits, **604** steric-score recomputations and **640** baseline/attempt native-score
measurements. Baseline native scores equal the saved WP2J scores. Every raw/rounded
candidate, including rejected geometry, retains its postselection reference score.
All **381 files** across eleven final replay arms match byte-for-byte. Replays,
readbacks, mathematical controls and repeated native requests are not new trials.

Verification: initial four source-math controls passed before native search; 41
extended tests passed, then an additional strict-grid metadata case brought the
final new set to **42**. Focused fast **824 passed**, three slow deselected; split
H/heavy/grid/seeded/renumber boundaries **223 passed**; full fast **3754 passed**,
three existing optional Torch/Biotite skips, 33 slow deselected. Ruff passes;
scoped strict typing passes 26 selected files and source-bound nonincremental
strict typing passes 27. The ambient-import run retains the same 793 legacy errors.
WP2F's large combined slow timeout remains unresolved, not causally closed here.
See [validation](chemical-reconstruction-conformer-validation.json) for pinned
source/report hashes, per-case comparisons and the negative ASH outcome.

Next: independently qualify additional structures/folds and ranking physics before
inferring broad native recovery. Empirical backbone-conditioned priors, continuous
orientation/refinement and H-bond/electrostatic/solvation terms need separately
declared evidence, not retrospective tuning to remove ASP21. Heavy/H ownership and
lineage integration, measured/partial H, actual-neighbour/polymer/link/metal/ring-
puckering chemistry, native-H/name quality and public/default/dependency/production/
data-license gates remain open. Production PyRosetta/data use still requires the
user's documented original-source/license permissions and explicit approval.

## WP2L first slice — rounded-heavy conditioning contract (2026-11-20)

[Protocol](chemical-reconstruction-heavy-h-v1.md) frozen as `bb292e8` from
`8d01121`; [parent ledger](chemical-reconstruction-heavy-h-parents-v1.json)
committed as `0475f76` BEFORE new staged-H scientific controls. The ledger pins
all twenty saved WP2K artifact hashes and eighteen selected rounded-heavy addition
signatures. GLY N/A and MET's absent eligible interior target remain unchanged.
No parent search was rerun for NEW scientific evidence. Old searches were replayed
only to check regression identity.

`heavy_conditioning.py` adds a distinct `HeavyConditionedRequest`, not a prepared
source subclass. It retains the original request and complete rounded candidate,
with a new tagged snapshot digest/origin. Numeric H inputs contain original
observations plus committed generated-heavy points; the original generated-heavy
mask remains explicitly available in outer provenance, never source-owned.
The DTO rejects nested/foreign source requests, source observation changes (even
numerically equal but digest-distinct scalar transport), wrong/missing/extra heavy
atoms, wrong elements/units/occupancy, non-PDB-rounded coordinates, changed original
chemistry/links/pH/constraint prefix, temporary/mobile atoms, and partial/mixed H
generation masks. Planning mode has EMPTY inner generation and cannot pass ordinary
completion admission; compiled H mode owns exactly the authoritative H mask.

The existing declared-state/extended/joint H compilers preserve this DTO through
`dataclasses.replace`, including PRO's intermediate plane constraint and existing
versioned compiler-origin settings. Those additions are transportable, not independently
qualified merely by the DTO: final staged auditing MUST recompile the declarative
plan and bind all authority/parameter settings before accepting any envelope.
This contract alone proves neither parent selection nor geometry.

The old fractional-scope checks were extracted verbatim into
`contracts.validate_fractional_context`; PreparedContextRequest fields/digests and
behavior stay unchanged. Admission recognizes fractional atoms for conditioning
ONLY through its exact original PreparedContextRequest proof. The shared scope
check still forbids selecting/owning/constraining/linking those opaque atoms.
`audit_source_transport` remains unchanged and rejects conditioning with
`prepared-context-request-required`, rather than certifying generated conditioning
as source observations. Generic HEAVY/HYDROGENS separation is unchanged.

Tests use synthetic inputs and tracked CCD/force-field authorities, not required
ignored parent outputs. Four source-math compile controls exercise SER/THR/PRO/HIP;
they do NOT constitute staged H proposals or scientific continuation of the frozen
parents. See the [conditioning validation record](chemical-reconstruction-heavy-h-conditioning-validation.json)
for tests, replay counts and remaining gates.

Verification: **29 new regression cases**, **115 focused passes**, **32 split
boundary passes**, and **3783 full-fast passes**, three existing optional skips,
33 slow deselected. The first full-fast run failed six CLI subprocess tests because
the absolute Python invocation omitted its environment's `bin` from `PATH`; that
failure is retained, and the complete rerun with corrected `PATH` passed without
code/test/threshold changes. All **381 old-control files** across eleven arms remain
byte-identical. Ruff, agent maps and generated CLI/GUI checks pass. Scoped strict
typing passes 27 files and source-bound nonincremental strict typing passes 28;
the ambient diagnostic multiset remains exactly the same 793 legacy errors.

**Next at the first slice (implemented in the continuation below), on the SAME branch:** implement the bounded strict known-DTO research
decoder; staged heavy→H constructor and separately role-labelled raw/rounded
envelopes; constructor-free lineage/rounding/geometry/CIP/pair-coverage auditor;
synthetic corruption/failure tests. Then execute the frozen twenty-entry continuation
from saved parents without resampling, retain every H refusal, independently read
back evidence and replay old controls. H failure must NOT trigger heavy reselection,
motion, new phase/state/mask or fallback. The committed rounded heavy checkpoint
must be identical in raw-H and rounded-H envelopes. The original WP2K ASH21
negative and WP2F timeout remain open. No public publisher, `.dat` sidecar, default,
dependency, production/data-license or broad native-H qualification is enabled.

## WP2L bounded continuation — staged heavy→H and independent audit (2026-11-20)

Built on `25a2484` under the unchanged `bb292e8` protocol and `0475f76` parent
ledger. The driver hashes ALL twenty retained parent artifacts before reserving
output and forwards only typed request/parent/authority, real incoming C, pinned
amide parameters and explicit backbone types. Retrospective native-reference/RMSD
fields stay in the I/O adapter and never enter construction. The scientific run
disables heavy search/rotation and SetDihedral; eighteen saved selected parents
continue to H planning without resampling. GLY remains N/A and MET still lacks the
declared eligible parent. No H-only fallback changes those denominators.

`research_decode.py` is a private known-root DTO transport, not pickle, `.dat`,
public JSON/API migration or source validation. Trusted annotation closure and
exact field sets select only known dataclass variants; JSON cannot name/load a
class. It rejects unknown/missing fields, wrong scalar/container/literal types,
bool-as-number, duplicate keys/set items and nonfinite/overflow numbers. Bounds
are 64 MiB, two million value/container nodes and depth 64. Driver reads are bounded
as well as hashed; absent/corrupt parents fail before reservation/science. Decoder
success alone establishes neither source provenance nor scientific acceptance.

`staged_completion.py` performs one complete-H proposal against the committed
ROUNDED heavy checkpoint. `staged_contracts.py` retains separate role-labelled
heavy and H arrays and source/parent/conditioning/H-candidate digest bindings.
"Raw" means raw H after rounded heavy, not globally unrounded heavy+H. Both
envelopes keep identical committed heavy points; original source observations are
untouched. Refused/partial H stages retain available evidence but have no accepted
combined result, no heavy reselection/movement and no phase/state/mask change.
There is no PDB materializer, `.dat` sidecar or public/default route.

`staged_audit.py` independently re-audits parent topology/grid/selection/baseline
metadata and recompiles H declarations/parameters (no coordinate generation), then
remeasures raw/rounded H geometry and conditional full-target CIP, exact rounding,
ownership/masks and envelope lineage. Proposal admission declaration was extracted
unchanged into the pure `construction.proposal_admission` helper for metadata
readback. All old control outputs remain byte-identical. Explicit pair accounting
covers every eligible pair with a generated atom in the source/heavy/H union,
excluding authoritative 1-2/1-3 paths: parent covers heavy-source/heavy-heavy, H
audit covers H-source/H-heavy/H-H. Incomplete H cannot claim complete-pair coverage.

Twenty entries are retained: **17 combined geometry/CIP passes**, **one rejected
H stage**, GLY N/A, MET unavailable. The seventeen passing cases contain **80 new
heavy atoms + 141 H**, with **173280** eligible pair counts summed across cases;
raw and rounded independently undergo seventeen H geometry and seventeen conditional
CIP audits each. Known heavy atoms stay fixed in both H envelopes. These are
case-specific geometry/integration results, NOT native-H accuracy or broad chemistry.

**Retained negative: 1UBQ A/PRO19.** The parent passes its existing local heavy
profile, but its generated CD lies **0.786839 A** from the N–incoming-C–CA plane
required by the declared PRO H interface (bound **0.05 A**). H admission returns
`invalid-input / fixed-constraint-conflict`, and H proposal returns
`unsupported / request-not-admitted` WITHOUT constructing H coordinates. There
are no raw/rounded H outputs for that case and it is NOT counted as a combined
pass. This exposes an incoming-peptide-interface gap in the local heavy proposal,
not a license to weaken the plane constraint, move observed atoms, reflect a
sidechain or replace the frozen parent after seeing an H refusal. Earlier PRO
failures and the separate WP2K ASH21 native-score worsening remain intact.

Independent readback disables all proposal/internal/joint-fit/search/rotation/
SetDihedral entry points, reconstructs DTOs from saved JSON and checks original/
prepared source transport for eighteen continuations, all twenty parent hashes, committed checkpoint roles,
H declarations/envelopes, measurements and the PRO refusal. Three retained readback
reports are byte-identical. Attempting to reuse the first report path hit the
expected exclusive-write refusal; the prior report was preserved and a fresh -02
report was used, followed by an identical -03 final readback with bounded evidence
reads. All **458 files** across twelve replay arms are byte-identical
(77 staged + 381 historical). Replays/readbacks and synthetic math tests are not
independent new scientific trials.

Verification: **44 new staged/decoder cases**, **73 passes** including the prior
conditioning tests; **254 split slow-boundary passes**; **3827 full-fast passes**,
three existing optional skips and 33 slow deselected. Scoped strict typing passes
31 selected files; source-bound nonincremental strict typing passes 32. The ambient
import diagnostic multiset remains exactly the same 793 legacy errors, not a global
typing pass. Ruff passes. See the
[validation record](chemical-reconstruction-heavy-h-validation.json) for evidence,
source/report hashes, all denominators and retained negative details.

**Next:** separately freeze an actual-neighbour/PRO-ring incoming-peptide geometry
strategy before new controls. Constrain construction to the real peptide interface
rather than retroactively weakening WP2L bounds or retrying heavy poses inside the
H stage. Preserve this v1 parent/refusal and compare explicitly labelled new requests.
Independent structures/folds/complexes, native-H/name/conformational quality,
measured/partial H, polymer/link/metal/ring chemistry, public publication/DatRecord/
defaults/dependency replacement and production/data-license gates remain open.
WP2F's large combined timeout remains unresolved. No scientific qualification or
production permission follows merely from this bounded integration milestone.
