# Reconstruction and modeling backends: research plan

Status: research, 2026-09-22. This note records current boundaries, candidate
architectures, and proposed benchmarks. It does not add a backend, dependency,
CLI option, or delivery commitment. PDBFixer and Salilab MODELLER remain
supported production baselines; nothing here deprecates or removes either one.

## Scope and terminology

Three different problems must not be collapsed into one feature:

1. **Heavy-atom completion** restores absent atoms within a known molecular
   component whose chemistry is sufficiently defined.
2. **Hydrogen placement and protonation** choose chemical variants and place
   hydrogens after the heavy-atom graph is complete.
3. **Residue, loop, or comparative modeling** creates whole missing residues or
   target-sequence regions from one or more structural templates.

Geometry regularization is also distinct from MD parameterization. A component
can have valid dictionary geometry without an MD-ready charge/parameter set,
and a force-field template does not by itself prove that an incomplete
component has been reconstructed correctly. The related
[whole-complex relaxation note](whole-complex-relaxation.md) defines that
boundary in more detail.

“Modeller” is ambiguous in this area:

- **Salilab MODELLER** supplies comparative modeling, LoopModel, `align2d`,
  optional SALIGN, refinement, and `molpdf` scoring.
- **OpenMM `app.Modeller`** edits topologies and places hydrogens. Replacing a
  Salilab MODELLER capability would not replace OpenMM `app.Modeller`.

## Current production boundaries

### PDBFixer use

The `legacy` preparation path uses PDBFixer for topology-oriented structure
repair: loading PDB input, finding missing residues, finding and replacing
nonstandard residues, removing heterogens when requested, finding missing and
terminal atoms, applying mutations, and materializing missing heavy atoms.
Read-only diagnostics also use its missing-component discovery.

DVBfixer does not accept direct unseeded `addMissingAtoms()` output. All current
materializing callers route through
`ffutils.geometry.rebuild_missing_atoms_with_retry`, which snapshots topology
and coordinates, retries explicit seeds, identifies additions by stable atom
identity, and rejects candidates with rebuilt D-C-alpha centers or immediate
close inter-residue contacts. Legacy preparation also preserves a strict call
order because PDBFixer topology-changing operations invalidate residue-object
keys stored by `findMissingAtoms`.

Hydrogen placement is not owned by PDBFixer. The legacy path uses protonation
evidence plus variant-aware OpenMM `Modeller.addHydrogens`, followed by DVBfixer
hydrogen-geometry repair. The opt-in `tleap-reduce` backend uses AmberTools
`tleap` for pure-protein heavy atoms and MolProbity Reduce for hydrogens. It is a
useful deterministic baseline, not a general ligand, glycan, PTM, or cofactor
solution.

### Salilab MODELLER use

`dvbfixer model` currently uses MODELLER to build sequence/SEQRES gaps, generate
and refine loop candidates, preserve selected heterogens through BLK entries,
and rank candidates by `molpdf`. DVBfixer independently owns sequence placement,
input pinning policy, identity/connectivity restoration, deterministic residue
numbering, candidate publication, and `.dat` provenance.

`dvbfixer homology` uses MODELLER for comparative construction and optional loop
refinement after target/template alignment. The GUI template-plan path first
fits selected chains into a shared frame and creates one authoritative mosaic
PDB; preserving that mosaic contract is more important than preserving a
particular modeling engine.

`dvbfixer salign` already demonstrates a narrower backend split: sequence-guided
Biopython C-alpha superposition is the default, while MODELLER SALIGN remains an
explicit optional engine.

## Design principles for experimental alternatives

Any prototype should be behind a backend-neutral request/result contract rather
than embedded directly into `prepare`, `model`, or `homology`.

- Preserve case-sensitive chain IDs, insertion codes, alternate locations,
  MODEL identity, explicit covalent links, and observed-versus-generated atom
  provenance.
- Separate discovery, chemical-authority resolution, coordinate generation,
  local regularization, and validation.
- Preserve observed coordinates by policy and measure any permitted movement.
- Generate no partial public output on failure. Publication should be atomic and
  include backend, version, dictionary/model-weight version, parameters, seed,
  score provenance, warnings, and validation results.
- Treat unsupported or ambiguous chemistry as an explicit result. Successful
  execution is not evidence of chemical correctness.
- Keep PDBFixer and MODELLER available as baseline backends throughout research
  and any later opt-in evaluation period.

## Chemistry-authority-first atom reconstruction

A universal “infer the intended molecule from incomplete coordinates” engine is
not scientifically defensible. Coordinates alone do not uniquely specify bond
orders, aromaticity, formal charge, protonation, redox state, stereochemistry,
metal coordination, or intended covalent attachments. Reconstruction should
therefore begin by resolving an authoritative component graph.

Proposed authority order:

1. User-supplied mapped SMILES, SDF, MOL2, or monomer CIF.
2. A versioned local wwPDB Chemical Component Dictionary entry.
3. A supported polymer, force-field, or monomer-library template whose exact
   component and links are known.
4. Otherwise an explicit `unsupported` or `ambiguous` result.

The existing SMILES-guided ligand path is a useful precedent: supplied chemistry
is authoritative, observed heavy-atom names and coordinates are retained, graph
matching must be unique, and chemically incompatible or externally bonded
residues fail rather than guess.

A proposed reconstruction flow is:

```text
normalized structure
  -> component inventory and explicit links
  -> chemical-authority resolution
  -> unique observed-atom-to-graph mapping
  -> missing-coordinate materialization
  -> restrained local geometry regularization
  -> graph, stereo, geometry, clash, and identity validation
  -> atomic publication with provenance
```

### Candidate materializers by component class

| Component class | Candidate route | Required limits |
|---|---|---|
| Canonical protein residue | Deterministic residue templates plus internal coordinates and rotamer sampling; compare against `tleap` and PDBFixer baselines | Preserve backbone anchors and peptide links; zero detectable D-C-alpha; deterministic seeded ranking |
| Protein hydrogens | Existing protonation evidence plus MolProbity Reduce or variant-aware OpenMM placement | Run only after heavy atoms; keep protonation decision separate from coordinates |
| Known isolated ligand | CCD/Gemmi graph and ideal geometry plus RDKit constrained embedding for missing coordinates | Require unique mapping, known microspecies, and no unexplained external bond |
| Glycan, PTM, covalent ligand | Component and link dictionaries with linkage-specific restraints | Proximity is not chemical authority; preserve attachment stereochemistry and link identity |
| Cofactor or metal system | Validated component/link/redox/coordination model supplied by user or curated library | No generic isolated-ligand fallback for unsupported states |
| Unknown incomplete heterogen | None | Return unsupported rather than invent chemistry |

Gemmi is a candidate adapter for CCD/monomer-library graph and restraint data;
RDKit is a candidate for constrained small-molecule conformer generation. Neither
adapter removes the need for identity mapping, chemistry-state selection, or
independent validation. Phenix/cctbx-style dictionary regularization and a
custom OpenMM geometry-restraint objective are possible regularizers, but their
coverage, deployment, licenses, weights, and objective scaling require separate
evaluation.

### Proposed reconstruction request/result

A request should contain normalized coordinates, stable atom/residue/component
identity, observed masks, resolved chemical graph, explicit links, movable
selection, coordinate tolerances, requested candidate count, and deterministic
seed. A result should contain complete candidate coordinates, atom mapping,
generated-atom provenance, backend evidence, validation findings, and an
explicit unsupported/ambiguous/failure classification.

This contract deliberately does not imply that `PDBFixer` discovery and
materialization must be replaced together. Read-only missing-component
discovery, canonical protein materialization, and heterogen materialization can
be evaluated or migrated independently.

## Experimental loop-modeling backends

A loop backend request should include the observed structure, complete target
sequence, shared affine sequence placement, explicit gap masks, fixed/movable
atom policy, retained heterogens and covalent links, candidate count, and seed.
Its result should include candidates, generated residue/atom provenance, score
and score provenance, warnings, validation, and stable identity mapping.

### Deterministic fragment/closure backend

A narrow custom backend could combine a fragment library, cyclic-coordinate
descent or kinematic closure, side-chain rotamer packing, and local
geometry/OpenMM refinement. Ranking should combine closure error, covalent
geometry, severe clashes, Ramachandran quality, rotamer quality, anchor drift,
and independent DVBfixer diagnostics. This route offers maximum control and an
offline permissive implementation, but fragment selection, closure robustness,
long-loop sampling, terminal gaps, and multi-chain context would require
substantial scientific engineering.

### OpenStructure/ProMod3

OpenStructure/ProMod3 is a candidate open-source comparative-modeling and
loop-building stack. Evaluation must cover installability, Python API stability,
license compatibility, offline deployment, multi-chain behavior, chain and
insertion-code preservation, heterogen/link retention, candidate scoring, and
whether fixed template coordinates can be held to DVBfixer's tolerance.

### Rosetta comparator

Rosetta applications or PyRosetta are scientifically capable loop-modeling
comparators, but source, binary, and deployment licensing must be reviewed
before any integration claim. They should not become an automatic core
dependency. A benchmark adapter could still establish an external quality
reference if its use and license are explicit.

### Independent template-constrained diffusion backend

A diffusion proof of concept should use only the general concept exemplified by
[PATCHR](https://github.com/DeepFoldProtein/patchr): PATCHR must not be imported,
vendored, wrapped, copied, or treated as an implementation dependency. DVBfixer
would define its own backend-neutral interface and independently implement and
evaluate the method.

Proposed algorithm:

```text
complete target sequence/topology
  -> fixed observed atoms versus generated missing-region mask
  -> all-atom diffusion sampling
  -> reinsert fixed experimental coordinates at every denoising step
  -> weighted Kabsch synchronization of generated and template frames
  -> localized second diffusion/refinement pass around gap junctions
  -> independent DVBfixer graph/stereo/geometry/clash validation
  -> publish only a passing candidate
```

Ordinary template conditioning in embeddings is not enough: a generative model
may still move supposedly fixed deposited atoms. Per-step coordinate reinjection
makes fixed atoms an explicit constraint, while weighted Kabsch synchronization
prevents global frame drift. A localized boundary pass should target peptide
continuity, bond/angle geometry, chirality, and clashes near the splice without
regenerating the entire observed structure.

Candidate all-atom diffusion engines, including Boltz/Protenix-class models,
should be reviewed as replaceable samplers rather than adopted through a
model-specific API. Source-code licensing and model-weight licensing must be
reviewed separately. The backend must record model and weight versions, random
seeds, hardware, precision, and any nondeterministic kernels.

Required diffusion cases include internal two-anchor gaps, terminal one-anchor
gaps, short and long gaps, multi-chain interfaces, and retained covalent
heterogens. Main risks are GPU/weight requirements, nondeterminism, movement of
fixed atoms, hallucinated chemistry, confidence miscalibration, length limits,
and weak handling of rare PTMs, glycans, cofactors, or metals.

The minimum ablation compares:

1. ordinary template conditioning;
2. conditioning plus per-step fixed-coordinate reinjection;
3. reinjection plus localized boundary refinement.

## Experimental homology-modeling backends

A homology backend must preserve target/template chain mapping, the
`template-plan` mosaic, shared reference frame, zero-based half-open masks,
multi-chain identity, and the rule that covered template coordinates remain the
authoritative starting geometry outside generated regions.

### OpenStructure/ProMod3

Evaluate the same deployment and identity questions as for loops, plus
multi-template comparative construction, substitutions, insertions, deletions,
multi-chain assemblies, and ranking. A successful single-chain protein example
is insufficient for the GUI template-plan contract.

### Custom mosaic-first builder

The existing `selected_template_mosaic.pdb` can be treated as the authoritative
known geometry. A custom builder would copy covered coordinates, materialize
mutated residues through residue/rotamer backends, send insertions and deletions
to a loop backend, regularize template junctions, then rank complete candidates.
This decomposes comparative modeling into bounded capabilities and avoids asking
independent templates to define incompatible global frames.

OpenMM can support local refinement but is not by itself a comparative-model
builder. A diffusion sampler could implement insertion/substitution generation
inside the same mosaic-first contract only after it demonstrates fixed-coordinate
adherence and linkage preservation.

## Proposed benchmark corpus

The benchmark should be versioned and stratified rather than summarized by one
success percentage.

- Canonical protein side-chain truncations with withheld known coordinates.
- Backbone/terminal atom omissions and multiple simultaneous missing residues.
- Ligands with authoritative SMILES/SDF and CCD entries, including aromatic,
  charged, and stereogenic examples.
- Glycan residues, PTMs, covalent ligands, and explicit disulfide/link records.
- Complete validated cofactors and intentionally incomplete/unsupported
  cofactors such as the FAD evidence documented in the relaxation note.
- Antibody insertion codes and case-distinct chain IDs.
- Short, medium, long, terminal, and multi-chain-interface sequence gaps with
  experimentally observed coordinates withheld for evaluation.
- Single-template, multi-template, and GUI mosaic homology cases.

For every case, store input provenance, reference provenance, exact masks,
expected supported/unsupported classification, and tool/dictionary/model
versions. Avoid evaluating only structures used to tune implementation rules.

## Acceptance gates

### Shared hard gates

- No unexplained atom loss, duplication, identity change, graph change, or link
  loss.
- No new stereochemical inversion; zero detectable D-C-alpha for supported
  L-protein output.
- Explicit failure with no partial publication for ambiguous or unsupported
  chemistry.
- Stable mapping back to `(model, chain, resid, icode, altloc, atom)` identity.
- Reproducible output under the declared deterministic mode, or measured and
  documented repeat-run variability when the backend is stochastic.

### Atom reconstruction metrics

- Heavy-atom recovery rate by component class.
- RMSD of withheld atoms after alignment on observed anchors.
- Bond-length, angle, planarity, chirality, and rotamer outliers.
- Severe intra-component and inter-component contacts.
- Maximum and RMS displacement of observed atoms.
- Correct unsupported/ambiguous classification, especially for covalent,
  cofactor, redox, and metal cases.
- Runtime, peak memory, external-process reliability, and deployment footprint.

### Loop and homology metrics

- Complete target sequence and correct chain/gap placement.
- Peptide closure and junction connectivity pass rate.
- Ramachandran, rotamer, chirality, and severe-clash quality.
- Flank/anchor and fixed-template atom displacement.
- Preservation of heterogens, explicit links, and template-plan masks/frame.
- Candidate diversity, repeat-run variability, and ranking enrichment against
  withheld experimental coordinates.
- Runtime and memory by gap length and structure size.
- Score provenance; raw scores from different objectives must not be compared as
  though they share a physical scale.

Diffusion evaluation additionally records fixed-atom RMSD, the fraction of
candidates recovered by boundary refinement, seed-to-seed diversity, and the
three-way conditioning/reinjection/refinement ablation.

Numerical pass thresholds must be fixed before running the comparison and
reported per chemistry/gap stratum. A backend should become production-default
only after it beats or complements the PDBFixer/MODELLER baselines on declared
coverage and quality, not merely because it completes more inputs.

## Staged research and integration sequence

1. Freeze benchmark fixtures, masks, metrics, and publication rules. Add no
   production backend yet.
2. Extract backend-neutral request/result prototypes outside public command
   behavior, preserving existing PDBFixer and MODELLER baselines.
3. Benchmark heavy-atom materializers independently from hydrogen placement and
   geometry regularization.
4. Compare three loop directions: OpenStructure/ProMod3, a narrow custom
   fragment/closure engine, and an independently implemented
   template-constrained diffusion proof of concept. Rosetta may serve as an
   explicit external comparator.
5. Evaluate homology construction separately, with the template mosaic as a
   non-negotiable contract.
6. Only after acceptance gates pass, expose an explicit experimental backend
   selector with provenance and no automatic fallback that could hide a failed
   chemistry/modeling decision.
7. Retain PDBFixer and MODELLER as supported compatibility/baseline backends
   until a separate, evidence-backed product decision changes that policy.

No production code or dependency change is implied by this research sequence.
