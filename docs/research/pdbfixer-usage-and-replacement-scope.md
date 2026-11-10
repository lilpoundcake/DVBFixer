# PDBFixer in DVBFixer: actual consumers, chemistry and replacement requirements

Status: source audit and offline discovery controls; no production replacement.
Branch baseline `4a328a6`; continuation of the
[causal study](pdbfixer-reconstruction.md) and
[handoff](../plans/pdbfixer-reconstruction-handoff.md).

Implementation continuation: the [first private ALA/SER/CYS slice](protein-completion-first-slice.md)
is now executable. The requirements and unclosed acceptance gates below still
apply; no production consumer has been replaced.

## Decision summary

**Do not design a replacement for the whole PDBFixer package before deciding
which responsibility is being replaced.** Our immediate need is a safe
**heavy-coordinate completion kernel**, not new protonation, file-format,
solvation, ligand-parameterization and homology-modeling implementations.

The actual dependency is wider than five historical preparation calls:
**11 constructor sites and 9 seeded-rebuild sites across 8 materializing
modules**. Diagnostics is a separate read-only consumer. Four of the nine
rebuild sites are in diffusion-related modules, including a private comparator.
Counts refer to static source locations, not distinct commands or runtime jobs.

The [machine-readable inventory](pdbfixer-usage-inventory.json) pins file hashes,
calls and state accesses. Its AST scanner ignores comments/docstrings, resolves
the actual explicit import aliases and follows the current `fixer`/`_fixer`
receiver convention. This is not a general Python dataflow/reachability proof.
Readiness subprocess strings and indirect command dispatch are reviewed separately
below; comments alone are not treated as evidence of a capability.

After the diagnostic-only transfer to the research branch, the active source has
a separate [source-binding amendment](pdbfixer-usage-inventory-stereo-diagnostics.json).
It preserves this original inventory and template probe unchanged, records the
new `diagnose/pipeline.py` hash, and requires identical PDBFixer imports, calls,
line locations, consumer paths and counts. No template experiment was rerun and
no scientific result was replaced. Current source tests check that companion;
historical verifiers use the preserved pre-diagnostics source context. See the
[static force continuation](stereo-force-projection.md) for the migration boundary.

## 1. Actual direct consumers

Paths are relative to `src/dvbfixer/`. Line numbers refer to the audited source.

| Consumer | Concrete use | Whole missing residues? | What leaves this stage |
|---|---|---|---|
| `prepare/pipeline.py:355–956` | Parse into topology/coordinates; discover gaps; normalize selected nonstandard residues; optionally strip heterogens; inventory and rebuild heavy/terminal atoms; subsequently place H with **OpenMM**, not PDBFixer | **Yes**, retained `missingResidues` from input sequence records; requested deletions are scrubbed first | Whole prepared structure and generated-atom/variant provenance |
| `prepare/mutations.py:555–573` | `applyMutations` renames residues and deletes atoms not in target template; later rebuilding adds the target's missing atoms | Mutation is not gap construction | Mutated heavy graph, then completed by `prepare` |
| `minimize/pipeline.py:524–575` | Keep-H branch discovers missing heavy/terminal atoms; calls rebuild only when those counts are nonzero; restores disulfide pairing; adds H if needed | Not explicitly disabled, but input is an OpenMM-written temporary PDB without copied SEQRES | Whole working topology for minimization |
| `minimize/pipeline.py:667–763` | Rebuild-H/variant branch strips protein H, discovers heavy defects, invokes rebuild, restores links, then runs variant-aware H placement | Same temporary-input qualification | Whole working topology for minimization |
| `protonate.py:512–560` | Heavy-atom completion before OpenMM variant-specific H addition | **Explicitly disabled** with `missingResidues={}` | Re-protonated structure |
| `top/pipeline.py:273–477` | On demand, strips a temporary structure to protein/selected solvent ions, completes heavy atoms, places H and extracts only required ASH/GLH/HIS-family protonation H | **Explicitly disabled** | **Only selected H coordinates** are copied to topology-building records; completed heavy atoms are temporary |
| `diagnose/pipeline.py:103–118`, `diagnose/structural.py:34–98` | Missing-heavy, gap and terminal discovery only | Reports, never inserts | Findings; no coordinate repair |
| `model/diffusion/observed_completion.py:104–140` | Completes observed canonical residues **before** a new diffusion request is frozen; checks exact identities/observed coordinates and selects only requested additions | **Explicitly disabled**, gaps remain diffusion-owned | Missing-heavy records only; original source records are reinserted unchanged |
| `model/diffusion/boundary_refinement.py:160–209` | Temporary completion of unrelated incomplete sidechains so the FF system can be built; checks whether any D repair would move protected observations | **Explicitly disabled** | Final refiner returns only the declared generated/mobile atom set; temporary completions are not automatically published |
| `model/diffusion/rfdiffusion_v1.py:554–600` | Materializes sidechains on generated backbone, then selects only `request.generated_atoms` before independent refinement/validation | **Explicitly disabled**; backbone/gap already exists | Declared generated-heavy coordinates; this adapter's presence does not imply current public-profile exposure |
| `model/diffusion/atom_completion_benchmark.py:53–114` | Private seeded-template comparator; reconstructs from source and target placement, not held-out coordinates | **Explicit source-bound gap list**, overriding discovery | Only requested coordinates via experimental materializer; extra OXT is recorded but not published |

Every production `addMissingAtoms` call is centralized in
`ffutils/geometry.py:1025–1137`. Its two raw call expressions cover the
terminal/no-interior-defect fast path and seeded retry path. The nine consumers
above invoke this wrapper, never an independent coordinate builder.

### Indirect consumers and lifecycle dependency

- `zbs.py` orchestrates `renumber → model → prepare → minimize`. It does not own
  another PDBFixer algorithm. Default whole-gap modeling is Salilab MODELLER.
- Homology/GUI/jobs can reach the same preparation or diffusion input/refinement
  paths; they are consumers of those services, not extra scientific kernels.
  Do not assume every Homology run invokes PDBFixer.
- `--backend tleap-reduce` bypasses PDBFixer **inside that preparation backend**.
  A subsequent `minimize` may still instantiate PDBFixer. Selecting this prep
  backend does **not** establish a PDBFixer-free end-to-end ZBS workflow.
- `doctor.py`, `profile_readiness.py:16–37`, installation profiles and package
  dependencies check PDBFixer availability/resources. Readiness executes a Python
  string to load `soft.xml` and check canonical templates; the AST call counter
  does not count that subprocess as a scientific constructor call.
- `prepare/glycan.py` and `.dat` construction consume/mutate the returned
  `fixer.topology`/`positions` object as a carrier. Replacing the package requires
  adapting that carrier contract, not just changing one import.

### What is used versus not used

Used methods: constructor; `findMissingResidues`, `findMissingAtoms`,
`findNonstandardResidues`, `replaceNonstandardResidues`, `removeHeterogens`,
`applyMutations`, and wrapped `addMissingAtoms`.

`diagnose` conditionally calls `findMissingTerminals` only if that attribute
exists. PDBFixer 1.12 populates `missingTerminals` in `findMissingAtoms`; the
conditional call is a compatibility branch, not a separate 1.12 API requirement.

**Not called by audited application code:** PDBFixer `addMissingHydrogens`,
`addSolvent`, `addMembrane`, `removeChains`, `registerTemplate` or
`downloadTemplate`. The constructor is used with local PDB files/file-like
objects, not the `pdbid`, URL or mmCIF APIs. CIF normalization is already owned by
`structure_input.py` using Gemmi/Open Babel; do not implement another reader.

There is one important indirect possibility: upstream `applyMutations` can call
`downloadTemplate` for an unknown target. Our parser maps known protonation
variants to parents but does **not** whitelist only the 20 canonical targets.
The lack of a direct `downloadTemplate` call is therefore **not proof of complete
network independence for every current mutation request**. A new service needs
an explicit offline authority/refusal rule, not an inherited hidden download.

## 2. PDBFixer state that current code depends on

| State | Meaning | Requirement for replacement/adaptation |
|---|---|---|
| `topology`, `positions` | OpenMM objects; later replaced after H addition and used for bonds, export and `.dat` | Preserve compatible topology or add an explicit adapter; keep original identity and coordinates |
| `missingAtoms` | Residue-object keys; values are template Atom objects, not just strings | New core should use stable identity plus expected names/elements; adapt diagnostics instead of imitating object-identity hazards |
| `missingTerminals` | Residue-object keys with terminal atom names, normally OXT in this path | Separate real terminus from a temporary fragment end, gap, cap or covalent link |
| `missingResidues` | `(chain index, insertion position) → residue-name list`; mutated by callers | Keep whole gaps separate from partial-atom requests; preserve explicit user deletion and target ownership |
| `nonstandardResidues` | Selected normalization targets; GLYCAM protein anchors filtered out | Do not silently convert retained linked/PTM chemistry to a parent amino acid |
| `platform` | Selected Reference/CPU/platform for rebuilding in some callers | Backend-independent result plus explicit optimizer/platform metadata; do not force production Reference |
| Imported keep sets | `proteinResidues`, `dnaResidues`, `rnaResidues` used in cap-aware heterogen stripping | Molecular retention policy is not equivalent to coordinate-template coverage |

`templates`, `sequences`, `modifiedResidues` and the CCD cache are not directly
manipulated by the audited DVBFixer consumers, but they affect upstream methods.
Repeated residue numbers, insertion codes and repeated physical chain IDs cannot
be safely represented by a naive `(chain, resnum)` replacement dictionary.

## 3. What molecular types does PDBFixer actually handle here?

Three different claims must remain separate:

1. **Read/retain:** a component can survive a structure-processing pass.
2. **Discover/reconstruct:** its intended atom set and new coordinates are known.
3. **Prepare for MD:** its protonation, links and FF parameters are supported.

None implies the next. In particular, **no missing-atom finding for an unknown
residue means “not assessed”, not “complete”.**

### Offline probe of the installed 1.12 stock templates

The source SHA matches the earlier pinned upstream source. No network or
coordinate generation was allowed. The probe finds **30 templates**:

- All 20 canonical protein amino acids.
- RNA: **A, C, G, U**.
- DNA: **DA, DC, DG, DT**.
- Peptide caps: **ACE, NME**.

One actual heavy atom was withheld from each template and discovery was rerun:
**30/30 withheld atoms were detected**. This tests name/template inventory only,
not reconstruction geometry, nucleotide stereochemistry or end-to-end commands.

The keep-classification lists also contain RNA `I` and DNA `DI`, but **neither
has a stock template** in this installation. `NHE`, `HOH`, `UNK`, `MSE`, `NAG`,
`DAN`, `HEM`, `ZN` and GLYCAM anchors `NLN/OLS/OLT` likewise have no stock
coordinate templates. Do not confuse a built-in substitution or retention rule
with reconstruction of that original component.

### Chemistry matrix

| Component class | Upstream stock coordinate role | DVBFixer route / own-code obligation |
|---|---|---|
| Canonical proteins, peptides, antibody chains | Heavy-atom and whole-residue templates | Primary own-kernel scope: sidechain/terminal completion; all Cα and sidechain stereocentres must be checked |
| Protein protonation variants / disulfides | OpenMM loading and DVBFixer parent/variant normalization plus canonical residue templates | Keep variant decisions separate; preserve CYS SG partner and atom count; no forced canonicalization of chemistry |
| ACE/NME caps | Stock templates, but stock heterogen stripping does not retain them automatically | Preserve existing cap-aware DVBFixer logic and link geometry; use `terminal_caps.py` for intentional cap creation |
| MSE and other modified residues with a parent substitution | Can normalize selected entries via table/MODRES, then complete the **parent** | This changes chemistry; not a claim to reconstruct native MSE/PTMs. Retain explicit normalization policy and protected GLYCAM anchors |
| Canonical DNA/RNA | Stock atom templates and nucleotide-specific terminal omission rules | An all-protein replacement cannot silently take over. Keep legacy routing initially or implement/qualify nucleotide sugar/base/phosphate geometry separately |
| Ligands, free carbohydrates and arbitrary HETATM | No general intended missing-atom inventory without a registered template | Existing RDKit/Open Babel/authoritative-SMILES H handling and ligand parameters are separate. Preserve components/context; do not invent missing heavy atoms |
| Glycan trees / glycoprotein anchors | No stock GLYCAM/NAG materialization here; can be retained as context | Existing GLYCAM/link logic owns bonds/H. Completing a protein atom near a glycan must preserve the attachment and steric context; rebuilding the sugar is another capability |
| Covalent ligands / PTMs | May be retained or normalized, not generally linkage-aware reconstructed | Require link-aware authority or explicit refusal; never treat as an isolated molecule by residue name alone |
| Cofactors, metal centres, ions | Can be retained, excluded, or used as steric context; not a coordination/redox reconstruction solver | Preserve exact atoms and links; unknown metal coordination is not a missing-protein-atom problem |
| Water | Retention controlled separately; H/solvation handled elsewhere | No stock heavy reconstruction need comparable to a sidechain; preserve/remove only by explicit policy |
| User-registered / downloaded CCD components | Upstream API can expand template coverage | No direct registration/download in current heavy-repair callers; arbitrary nonstandard mutation targets are the indirect exception described above |

The code route for a complex being present is **not scientific acceptance for
all of these chemistries**. In particular the existing canonical L-Cα detector
does not qualify D-ribose, anomers, ligand stereocentres or metal geometry.

### Existing own nonprotein reconstruction must not be duplicated

There is already `nonprotein_reconstruction.py`,
`domain/nonprotein_reconstruction.py`, `ccd.py` and
`pdb_component_reconstruction.py`/`reconstruct-component`.
The implemented opt-in Class A route reconstructs **known isolated components**
from an authoritative local graph/ideal geometry with fixed observations and
independent checks. Linked glycan/PTM/covalent components and metal/redox classes
are explicitly out of that slice. It is not an automatic fallback in `prepare`.
See the [nonprotein plan](../plans/nonprotein-reconstruction-parameterization.md).
Extend/reuse its authority and identity contracts where appropriate, rather
than designing a second generic ligand builder under “replace PDBFixer”.

## 4. Concrete own-code functional scope

### Option A — recommended first: replace only unsafe coordinate generation

Leave PDBFixer discovery/normalization behind an adapter initially. Implement an
opt-in canonical-protein heavy-atom completion service, with no public default
change. This **does not remove PDBFixer as a dependency**.

Required functionality:

- **F1 — Exact request:** source digest; full residue/atom identity, including
  model/physical occurrence where needed; canonical parent and separate variant;
  expected graph; explicit missing-atom mask; fixed observations; true/temporary
  termini; explicit links; candidate/resource policy. Reuse canonical atom tables
  and variant/identity services, not copied chemistry lists.
- **F2 — Partial sidechain completion:** one missing atom through an entire
  sidechain; explicitly L Cβ construction; internal coordinates and rotamers;
  rings/PRO closures and THR/ILE sidechain stereo when admitted. Existing distal
  atoms constrain the solution just as strongly as backbone anchors.
- **F3 — Terminals/caps:** OXT and supported cap defects only when chemistry and
  chain topology require them; no fictitious OXT beside an internal gap or cap.
  A minimal sidechain-only milestone may refuse these, but cannot then claim
  drop-in equivalence with the current wrapper.
- **F4 — Missing backbone:** N/CA/C/O defects need their own anchored geometry
  and peptide-closure support. Initial v1 should refuse them explicitly rather
  than silently pretend sidechain support covers them. Complete replacement
  requires either this capability or explicit routing to an approved engine.
- **F5 — Fixed environment:** evaluate neighbouring chains, retained ligands,
  sugars, ions, disulfides and covalent anchors; never delete the context to get
  a pass. Unparameterized environmental atoms still need geometric exclusion.
- **F6 — Constrained search/ranking:** sample actual torsional freedom, reject
  stereo/graph/geometry violations before ranking, handle coupled defects and
  retain ambiguity/candidates. No reference coordinates in candidate generation.
- **F7 — Generated-only refinement:** preserve observations; explicit guards for
  every admitted stereocentre and final degeneracy rejection. No repair that
  reflects an observed sidechain under the guise of adding one missing atom.
- **F8 — Result/publication:** only requested atom identities/coordinates;
  exact generated provenance, deterministic tie-breaking, validation findings
  and refusal reasons; caller adapters own OpenMM topology and canonical
  `DatRecord` publication. Temporary FF-only atoms must be labelled and excluded
  from scientific outputs and `.dat`.
- **F9 — Routing:** distinguish protein completion, whole-gap modeling,
  nucleotide completion and nonprotein Class A reconstruction. Unsupported
  chemistry must be routed explicitly or refused, never silently dropped.

**Not part of this kernel:** PROPKA/pH/tautomer selection; H placement; CIF
normalization; MD parameter generation; solvent/membrane construction; GUI;
sequence alignment; publication storage; arbitrary whole-loop prediction.
Those responsibilities already have owners or require separate capability work.

### Option B — remove the entire dependency later

In addition to F1–F9, replace/adapt:

1. Missing-heavy and true-terminal discovery and diagnostic schemas, including
   unsupported/not-assessed states rather than false clean results.
2. Gap inventory and exact sequence placement; reuse `sequence_alignment.py` and
   existing MODELLER/diffusion ownership. Direct `prepare` currently can insert
   SEQRES gaps: changing that behavior needs an explicit compatibility decision.
3. Nonstandard-residue/MODRES policy and mutation atom-graph edits. Existing
   deletion/link cleanup in `prepare/mutations.py` remains reusable. A new precise
   insertion-code mutation syntax would be a separate user-visible enhancement,
   not an implicit replacement requirement.
4. Heterogen/water/cap retention, nucleotide fallback/support and topology/bond
   reconstruction. OpenMM can remain the parser/topology layer; replacing
   PDBFixer does not imply replacing OpenMM or Salilab MODELLER.
5. All nine scientific caller adapters, read-only diagnostics, state-carrier
   consumers, Doctor/readiness, packaging and installation profiles.

Switching one preparation call, or passing six Cβ controls, closes none of these
dependency-removal gates by itself.

## 5. Required tests and what exists now

The following is an acceptance specification, **not a claim that the own
builder or all these tests exist**. Existing files are regression assets to
retain; passing them with the legacy backend does not qualify a replacement.

| ID / requirement | Concrete new-builder tests and pass condition | Existing assets / remaining gap |
|---|---|---|
| T1 / F1 identity, inventory | Each of 20 AAs: delete each nonterminal sidechain atom individually, terminal subtrees and whole sidechain; exact missing-set equality; no unexpected atom; GLY adds no CB | `test_diffusion_atom_completion.py`, canonical tables; exhaustive own-builder matrix **needed** |
| T2 / F2 stereochemistry | Every generated chiral centre L/correct dictionary stereo; THR/ILE Cβ stereo, branches, PRO/rings; planar/undefined/nonfinite centres reject; rotation/translation and atom-order controls | `test_ffutils_chirality.py`, causal study, `test_pdbfixer_research.py`; full non-Cα own-builder checks **needed** |
| T3 / F2 retained distal atoms | Remove an internal sidechain atom while retaining distal atoms; satisfy both constraints without moving any observation; impossible geometry refuses | Partial atom inventory/materialization tests exist; actual multi-anchor coordinate solver tests **needed** |
| T4 / F3 ends and caps | OXT-only missing, true N/C ends, cap-linked ends, internal gap fragments, ACE/NME; no terminal atom added inside a peptide; correct cap/peptide bonds | `test_terminal_caps.py`, `test_prepare_broken_geom.py`; own completion parity **needed** |
| T5 / F4 backbone | Explicit N/CA/C/O deletion cases near observed neighbours; admitted cases satisfy closure/amide geometry; initial unsupported cases fail before output, not fallback to sidechain solver | Current kernel capability **not implemented**; own admission/solver tests **needed** |
| T6 / F5 mixed context | Protein near ligand/glycan/partner, disulfide pair and crowded CYS cluster; original partner atoms/links unchanged, SG degree ≤1, no severe generated/fixed clash | `test_conect_false_peptide.py`, `test_convert_glycam_roundtrip.py`, `test_zbs_e2e.py`, `test_diagnose_steric.py`; own-builder paired context controls **needed** |
| T7 / F7 preservation | Same coordinates and byte-identical retained records where the caller promises this; case-distinct chains, insertion-code siblings, MODEL/altloc/occurrence, ANISOU, links and CONECT endpoint identity | `test_diffusion_observed_completion.py`, `test_diffusion_atom_completion_validation.py`, `test_altlocs.py`, `test_conect.py`; adapter-specific assertions required |
| T8 / F6 quality and ambiguity | Bond lengths/angles, peptide omega/planarity, chi/rotamer plausibility, ring closure, all relevant stereo, clashes; rank without held-out fit; report multiple plausible outcomes and failures | `test_diagnose_protein_geometry.py`, benchmark service and withheld fixtures; source-disjoint own-builder qualification **needed** |
| T9 / F8 ownership | Exact generated-only `.dat`, no known atoms marked new; FF-only OXT/completion never leaks; exceptions/cancellation retain evidence and publish no accepted output | `test_diffusion_atom_completion_publication.py`, `test_diffusion_atom_completion_validation.py`, `test_prepare_input_file_safety.py`; own result adapter parity **needed** |
| T10 / F9 other chemistry | DNA/RNA/caps use declared supported route; unknown ligand or I/DI is not labelled complete; PTM/link/metal class refuses or uses its explicit authority; no hidden network | Offline 30-template discovery probe and `test_nonprotein_reconstruction.py`, `test_ccd_component_reconstruction.py`; own routing tests **needed** |
| T11 / all callers | Replace/mock the backend separately in prepare, both minimize branches, protonate, top H extraction, observed completion, FF-only completion and RF sidechain materialization; retain private comparator as baseline | `test_seeded_atom_rebuild.py`, `test_diffusion_workflows.py`, `test_diffusion_boundary_refinement.py`, `test_diffusion_rfdiffusion_v1.py`, `test_diffusion_observed_completion.py`; shared adapter contract suite **needed** |
| T12 / mutations/normalization if Option B | Canonical substitution grows/shrinks sidechain; GLY/PRO changes; variants preserve chemistry; deletion is not reinserted; MSE normalization explicit; glycan/disulfide-dependent cleanup retained; unknown targets do not download silently | Existing mutation application code is not an exhaustive regression suite; dedicated own-edit tests **needed** |
| T13 / reproducibility/resources | Same declared platform/seed produces reproducible candidate policy; CPU/Reference/GPU geometry checked separately; missing resource/timeout/resource limit is failure, not pass; no broad skip hiding damage | `test_seeded_atom_rebuild.py`, `test_profile_readiness.py`, `test_doctor.py`; prospective own-runtime controls **needed** |
| T14 / hydrogen separation | Heavy-kernel replacement still calls the selected H engine only after heavy completion; test variants, sibling-H, pathological SER and invariant heavy coordinates | `test_prepare_propka_integration.py`, `test_prepare_icode_variants.py`, `test_prepare_broken_geom.py`; replacing PDBFixer alone does not repair OpenMM H placement |

### Benchmark protocol before implementation acceptance

- Freeze molecule classes, missing masks, candidate budget, scoring/geometry
  gates, source-disjoint development/test split and seeds **before** generation.
- Include artificial deletions from complete held-out structures and authentic
  incomplete inputs (e.g. tracked 5J13/8CZ8); authentic absent atoms have no known
  reference coordinate, so report geometry/uncertainty instead of invented RMSD.
- Compare sidechain completion against seeded PDBFixer and applicable `tleap`;
  compare whole gaps against **full legacy MODELLER**, not just PDBFixer's arc.
- Report refusals, attempted failures, final geometry passes, fixed-atom drift,
  deposited-coordinate errors and wall time separately. Chemical name-symmetry
  scoring is separate from exact publication identity; no best-seed selection by
  reference and no pooled score that hides a stereochemical failure.

## 6. Audited conclusion and next action

**Minimal useful development target:** a private canonical-protein
sidechain/terminal heavy-atom service with exact masks, stereochemical internal
coordinates, constrained rotamer search, immutable observations, retained complex
context and independent rejection. PDBFixer can initially remain for discovery
and explicit legacy routes. DNA/RNA and linked nonprotein rebuilding are separate
capabilities, not tasks to solve by pretending all HETATM are protein.

The user-requested priority is therefore **boundary and test specification
first, implementation second**. The already executed Cβ primitive is a useful
building block, not justification to expand scope or change defaults.

### Verification of this audit

- Static inventory: 11 constructors, 9 rebuild call sites, 27 method expressions
  including the conditional terminal-discovery compatibility call.
- Offline probe: PDBFixer 1.12.0, same pinned source SHA; 30 templates, 30 withheld
  atom detections; no downloads and no new coordinates generated.
- The earlier temporary CPU runtime no longer starts (Python encoding/bootstrap
  failure). It was not repaired. This inventory-only probe used the existing
  general `dvbfixer` interpreter; it does **not** certify its `soft.xml`, optimizer
  or current scientific readiness. Earlier causal measurements remain historical
  retained evidence, not freshly re-executed results.
- New inventory regressions are stdlib-only and do not import all command
  parsers or require the scientific stack. Exact test/check results are recorded
  in the live handoff after execution.
