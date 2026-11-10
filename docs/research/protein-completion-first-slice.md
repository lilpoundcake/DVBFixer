# Private protein completion: first implementation slice

Status: **private executable prototype**, not a registered preparation backend.
User authorized starting implementation after the
[usage/requirements audit](pdbfixer-usage-and-replacement-scope.md).
Legacy remains default and all existing PDBFixer consumers remain unchanged.

Continuation: an [optional offline single-chi adapter](protein-completion-rotamer-adapter.md)
is now implemented. The measurements and source hashes below describe the
original staggered control; real Dunbrack artifact acquisition/qualification is
still open. Do not relabel these earlier results as empirical-library output.

## Implemented, not merely proposed

`src/dvbfixer/protein_completion/` contains:

- `contracts.py`: frozen observation/request/result types, exact domain atom
  identity, source SHA, observation/context/link-bound request digest, explicit
  missing mask and generated-only candidates. Limits and duplicate checks apply
  before construction.
- `geometry.py`: own L-Cβ construction from N/CA/C projection constraints and
  internal-coordinate extension with explicit bond length, angle and signed chi.
  No PDBFixer template coordinates, random initialization or Cartesian minimizer.
- `service.py`: deterministic proposals for **ALA, SER and free CYS sidechains**;
  retains every attempted candidate and chooses only an independently passing one.
- `validation.py`: measurement-only readback of exact generated ownership,
  expected graph, supported local bond/angle geometry, positive nonplanar Cα
  volume, declared torsion and generated/context sterics. It does not regenerate
  a supposedly correct proposal and compare it with itself.

This service imports without NumPy, Biopython, OpenMM, PDBFixer or RDKit, verified
in a separate `python -S` process. It reuses `ComponentInstanceRef`/`ExactAtomRef`,
the existing canonical heavy-atom catalog and pure diagnostic measurements/radii.
It introduces no CLI flag, registry entry, dependency change or output writer.

### Exact first capability

| Request | Outcome |
|---|---|
| ALA missing Cβ; intact N/CA/C/O | Build explicitly L Cβ, validate against fixed source/context |
| SER/CYS missing terminal OG/SG, observed Cβ | Enumerate three signed chi proposals, preserving observed Cβ |
| SER/CYS missing Cβ and terminal side atom | Build L Cβ then enumerate terminal proposals |
| Complete supported residue, no requested additions | `unchanged` only after scoped validation; malformed existing geometry rejects |
| Other residues, variants, missing backbone or terminal/cap creation | Explicit `unsupported`; no automatic legacy fallback |
| Missing Cβ with retained distal OG/SG | Explicit closure-required refusal; never move the distal observation |
| Disulfide, glycosylation or other nonpeptide link touching the target | Link-aware solver required, so refuse |
| Explicit ordinary peptide links | Preserve supplied topology, validate target C–N link length and partner uniqueness; no bond inference |
| Partial occupancy, unresolved altlocs, multiple models, unknown element/radius or nearby metal | Explicit refusal; no normalization, discarded context or invented coordination |
| Impossible anchors, D/nearplanar Cα, inconsistent mask or unacceptable candidates | No selected result; malformed contracts raise, unsupported scope refuses, failed geometry rejects |

Only new coordinates are returned. Original coordinates are not an optimizer
variable and cannot be supplied by the candidate as a replacement. Generated
atoms must exactly equal the requested mask; extra known atoms, wrong elements,
duplicates, nonfinite coordinates and stale request digests fail independent audit.
Chain case, insertion codes, model identity and physical occurrence remain exact.

### Sampling and ranking are intentionally limited

SER/CYS use **−60°, +60°, 180°**, a declared staggered torsion grid. This is
**not a backbone-dependent empirical rotamer library**. No Dunbrack artifact was
downloaded, distributed or assigned made-up probabilities. Fixed v1 illustrative
construction geometry and validation tolerances are explicit in code; they are
not a calibrated residue/backbone-dependent restraint dictionary or MD force field.

First reject incorrect geometry/stereo/clashes. Among passing proposals, rank by
worst nonbonded clearance capped at zero, with stable grid order for ties. There
is no purported physical energy, H-bond score, entropy or likelihood, and no
held-out reference coordinate reaches the generator or ranking function.

## Executed controls

`scripts/benchmark_protein_completion.py` declares eight masks before generation
on the same tracked **1UBQ chain-A 1–72 protein-only projection**:

| Mask | Final scoped geometry | Selected deposited-coordinate RMSD |
|---|---|---:|
| ALA28 / CB | pass | 0.040 Å |
| ALA46 / CB | pass | 0.082 Å |
| SER20 / OG | pass | 1.966 Å |
| SER20 / CB + OG | pass | 1.443 Å |
| SER57 / OG | pass | 0.023 Å |
| SER57 / CB + OG | pass | 0.047 Å |
| SER65 / OG | pass | 2.459 Å |
| SER65 / CB + OG | pass | 1.599 Å |

The service sees **all retained projection atoms**, minus the declared withheld
ones, plus the explicitly declared peptide neighbours. This is not intact
deposited-complex acceptance: the projection excludes residues 73–76, waters and
other excluded source records as stated, without rewriting their occupancies.
PDB parsing in this script is a bounded fixture adapter, not a new production
format boundary. CYS currently has synthetic contract/geometry tests only, **not
a real-source scientific cohort**.

**20 proposals:** 16 pass and four SER65 alternatives fail sterics. All eight
masks have a selected candidate with positive Cα volume. All rejected attempts
and separately computed held-out scores remain retained. Zero observation motion
is established **by generated-only representation**, not by a claim of auditing
a published PDB: there is no serializer or scientific publication in this slice.

The large SER errors are a meaningful negative result: scoped geometric validity
does not establish native conformation, effective ranking or superiority to the
existing engines. No PDBFixer/`tleap` comparator was rerun for these eight masks.
These cases reuse a previously studied protein and are not source-disjoint
qualification or an independent holdout discovery.

Portable [evidence](protein-completion-first-slice-evidence.json) binds original,
projection and masked-input digests; exact executed module/script hashes;
generated coordinates, every audit, selected indices and all reference scores.
Full immutable request snapshots and readback are local in
`tmp/pdbfixer-research/own-kernel-controls-03.json`; earlier-run evidence is retained
separately and not counted as independent scientific repetitions.

## Tests and execution

New `tests/test_protein_completion.py` covers:

- Exact ownership, known-atom immutability and source/context digest binding.
- Independent rejection of forged observations-as-additions, duplicates,
  wrong elements, missing atoms, mirror/planar/nonfinite geometry and false chi.
- Whole supported sidechains and terminal-only sidechain defects.
- Preserved distal-atom closure refusal, backbone defects and linked-SG refusal.
- Peptide links, distant fabricated bonds, metal and unsupported context refusal.
- Context-dependent rejection/ranking and all-failed attempt retention.
- Source D centres are rejected, never repaired by reflecting an observed atom.
- Atom order, exact identity, rigid motion and **100 oblique rotations per
  supported residue**; generated bond angles and signed chi are measured.
- Readback against withheld real fixture atoms without using them for generation.
- Import without scientific dependencies.

Executed narrow suite: **131 passed in 7.72 s**, including existing identity,
diagnostic geometry and experimental completion contract/validation regressions.
Targeted Ruff passed; mypy passed for all five package files with
`--follow-imports=silent` (no claim that imported legacy modules were rechecked).
Agent documentation maps passed. No full scientific-suite/production acceptance
is inferred from these focused tests.

Initial broader test attempts in the base interpreter failed on existing tests'
OpenMM imports, not new-kernel assertions. The complete 131-test focused run then
passed in the existing general scientific interpreter using project-local pytest.
No skip was added and no external scientific environment was repaired. Type tools
were installed only under `tmp/pdbfixer-research/type-tools/` in this worktree.

```bash
# Stdlib-only new kernel + inventory tests (normal pytest installation required)
PYTHONPATH=src python -m pytest tests/test_protein_completion.py tests/test_pdbfixer_usage_inventory.py -q

# New exclusive destinations; raw requests and portable summary kept separate
PYTHONPATH=src python scripts/benchmark_protein_completion.py --output tmp/pdbfixer-research/new-controls.json --summary tmp/pdbfixer-research/new-summary.json
```

## Remaining gates before replacing even one production caller

1. **Authority/data:** pin appropriate residue-specific geometry and empirical
   backbone-dependent rotamer data with attribution/licensing and versioned
   parsing, rather than describing the current three-angle grid as a library.
2. **Chemical coverage:** expand beyond ALA/SER/free CYS; add branched and aromatic
   residues, secondary stereocentres, PRO/ring closure and observed-distal-atom
   constraints. Retained links need dedicated scope, not new exceptions to clashes.
3. **Geometry/refinement:** calibrated stereo/steric/rotamer tests, explicit
   candidate budgets and generated-only constrained refinement with independent
   postchecks. Keep failures visible; do not tune thresholds on these controls.
4. **Adapter/integrity:** authoritative full-source inventory, complete-context
   and link extraction, source SHA verification, global chirality/geometry gates,
   byte-preserving insertion and canonical `.dat` ownership. The present typed
   kernel trusts caller-provided observations/context and **does not verify
   source bytes or global chemistry outside the selected residue**.
5. **Qualification:** prospective source-disjoint natural and deletion cohorts,
   matched PDBFixer/`tleap` comparators, failure counts and runtime. Hydrogen
   engines and whole-gap MODELLER ownership remain separate.
6. **Integration:** only after these gates, explicitly adapt and test individual
   consumers. Do not patch the central retry helper to divert all nine call sites
   to a three-residue capability or remove legacy resource/dependency checks.

This is the first executable implementation milestone, not completion of the
replacement plan or authorization to use the prototype as a production backend.
