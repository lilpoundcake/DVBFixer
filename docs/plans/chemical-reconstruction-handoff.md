# Chemical reconstruction — active development handoff

Date: 2026-11-19. Status: approved development direction; implementation and
scientific acceptance are incomplete. This document is the **starting point for
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
- Current task is branch/handoff preparation, NOT implementation or another
  scientific sampling run. Subsequent implementation starts with WP0 below.
- No main merge, release, default-backend switch, dependency removal, installed
  OpenMM patch or automatic upstream publication follows from this handoff.

Documentation scope (six files): this handoff and its validation JSON; discovery
updates in `AGENTS.md`, the historical PDBFixer handoff and the existing agent
task record; one narrow `.gitignore` entry for the retained QA artifacts. User's
direct handoff request authorizes these bounded changes.
Pending links: none; proposed files/types below are not existing APIs.

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

The first implementation may use a proposed private `chemical_reconstruction`
namespace; **that namespace/API does not exist yet and its layout must be decided
in WP0**. Reuse/adapt domain/nonprotein graph concepts rather than inventing a
second ligand authority or renaming the old residue contract "universal".

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

**WP0 — first coding session; do this next, on this SAME branch.**

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

**WP1 — shared geometry primitives / first cross-chemistry slice.** Adapt tested
internal-coordinate/anchor primitives and a Class A example to the agreed graph
contract; independent validators must measure rather than call the generator.
Close source and rounded geometry before selecting conformers. Existing ideal
lengths/angles are provisional until a documented geometry authority is chosen.

**WP2 — expand protein reconstruction and H in parallel.** Whole/partial sidechains,
multi-χ and rings, then explicit backbone/terminal cases; authoritative geometry
and stereo. Geometry-preserving H construction and orientation solving are an
equal exit condition, not a promise after heavy integration. Keep empirical prior
and physical/context scoring separate from hard validity.

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

Resume with WP0. Do not restart completed PDBFixer/Rosetta/OpenMM investigations,
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
