# Chemical reconstruction — active development handoff

Date: 2026-11-19. Status: WP0/WP1/WP2A/WP2B/WP2C and bounded WP2D metric geometry
checked; synthetic amide controls and final fast suite completed. General
reconstruction and scientific acceptance remain
incomplete. This document is the **starting point for
future reconstruction work**, not a claim that a universal tool has shipped.

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
