# Parameterization

Status: partial

Verified on: 2026-09-18

Verified at commit: `425f290eb85760246766f1d1500e51672c640b2d`

## Purpose

This context describes how retained nonstandard components are routed toward a
force-field treatment. The domain vocabulary prevents a residue name from being
treated as proof that generic parameters are scientifically valid.

It does not select the protein force-field family, repair missing cofactor
atoms, resolve metal coordination or redox state, or perform geometry-only
refinement.

## Scope

`src/dvbfixer/domain/parameterization.py` owns a pure five-route vocabulary and
the precedence rules in `classify_parameterization`. The function consumes
caller-supplied evidence flags; it does not inspect structures or prove that
the evidence is valid. `domain/__init__.py` re-exports the enum and classifier.

Integration is deliberately partial. `lig_params.py` calls the policy only to
block named complex cofactors from generic GAFF. The other route concepts have
working CLI or adapter behavior in places, but those paths do not consume a
single domain decision object.

## Capabilities

| Policy route | Status | Meaning | Active behavior |
|---|---|---|---|
| `EXPLICIT_STRIP` | vocabulary only | Caller explicitly excludes heterogens from parameterized minimization | `minimize --strip-heterogens` implements a separate strip-and-splice path; it does not call the classifier |
| `USER_TEMPLATE` | vocabulary only | Caller supplies a validated component template | repeatable `minimize`/`zbs --extra-ff FILE` appends OpenMM XMLs; no production caller passes `user_template=True` |
| `EXACT_NATIVE_TEMPLATE` | vocabulary only | Caller has established an exact compatible native template | runtime screening starts from template names and OpenMM later checks atoms, bonds, and external bonds; no production caller passes `exact_native_match=True` |
| `UNSUPPORTED_COMPLEX_COFACTOR` | implemented guard | Generic isolated-ligand GAFF is unsafe | `lig_params.build_ligand_generator` rejects `HEM`, `HEME`, `FAD`, `FMN`, `NAD`, and `NAP`, including non-strict auto mode |
| `GAFF_CANDIDATE` | vocabulary default | Remaining component may be attempted as an isolated organic ligand | minimize's adapter attempts GAFF2 for unknown residue names, but does not ask the classifier for this route or prove isolation, completeness, organic chemistry, or metal absence |

Classifier precedence is `strip` before `user_template`, then
`exact_native_match`, then the complex-cofactor denylist, with
`GAFF_CANDIDATE` as the default. Matching the denylist is case-insensitive.

## Entry Points

| Concern | Start symbol |
|---|---|
| Change route vocabulary or precedence | `src/dvbfixer/domain/parameterization.py::classify_parameterization` |
| Change retained-ligand GAFF screening | `src/dvbfixer/lig_params.py::build_ligand_generator` |
| Change minimize routing or fallback | `src/dvbfixer/minimize/pipeline.py::minimize` |
| Change standalone GROMACS ligand generation | `src/dvbfixer/parametrize.py::main` |
| Change `top --acpype` dispatch | `src/dvbfixer/top/acpype.py::run_acpype_mode` |
| Change AMBER+GLYCAM export | `src/dvbfixer/acpype_export.py::export_gromacs` |

## Contracts

- Route classification is pure and returns one `ParameterizationRoute`.
- Boolean route evidence is authoritative input to the classifier; deriving and
  validating that evidence remains an application-service responsibility.
- A native-template name is only a screening hint. OpenMM template matching is
  the compatibility check for the actual atom and bond graph.
- Generic GAFF is for isolated, single-residue small molecules. It does not
  parameterize cross-residue ligand bonds, covalent attachments, or metal
  coordination.
- `--extra-ff` means the user owns validation of the supplied OpenMM XML; files
  are loaded after the selected base force field.
- Explicit `--parametrize-ligands` uses `strict=True`. Automatic minimize
  routing uses `strict=False`, except the complex-cofactor guard raises before
  optional dependency checks in both modes.
- `--strip-heterogens` preserves heterogens in final output by splicing their
  input coordinates back; it is not whole-complex minimization.

## Invariants

- Never infer exact compatibility from residue name alone.
- Never send known complex cofactors through generic isolated-ligand GAFF.
- Never silently weaken an explicit request for ligand parameterization.
- Preserve `PDBFile.writeFile(..., keepIds=True)` in ligand extraction.
- Build ligand SDFs from a heavy-atom-only Open Babel submolecule, perceive its
  bonds, then reattach existing hydrogens with bond order one.
- Reuse `ffutils.ligand_valence` for ionizable-group and known-double-bond
  corrections; do not duplicate this chemistry in an adapter.
- Do not treat universal-force-field refinement as MD-ready parameterization.
- `acpype_export.add_glycam_bonds` must keep its direct imports and must not
  hide failures behind a broad exception handler.

## Callers

- `minimize.pipeline.minimize` scans retained heterogens against base OpenMM
  template names and invokes `build_ligand_generator` when retained heterogens
  are present. It uses `strict=True` only for explicit
  `--parametrize-ligands`.
- `minimize.pipeline.main` validates and appends every `--extra-ff` path before
  minimization.
- `zbs.py` forwards `--parametrize-ligands`, repeatable `--extra-ff`, and
  `--strip-heterogens` to minimize; it owns no parameterization policy.
- No production caller currently asks `classify_parameterization` for
  `EXPLICIT_STRIP`, `USER_TEMPLATE`, or `EXACT_NATIVE_TEMPLATE`.

## Adapters

- `lig_params.py` is the OpenMM in-process GAFF adapter: topology to a temporary
  PDB, heavy-atom Open Babel submolecule to SDF, OpenFF `Molecule`, then
  `GAFFTemplateGenerator(forcefield="gaff-2.11")` with AM1-BCC and a JSON cache.
- `parametrize.py` is a separate single-molecule CLI adapter: antechamber,
  parmchk2, tleap, and ParmEd produce GROMACS `.itp`, `.gro`, and restraints;
  RESP can use Gaussian, PSI4, or PySCF. It rejects multi-residue and
  metal-containing PDB input but does not call the domain classifier.
- `top --acpype` bypasses the RTP topology builder and exports a fixed
  AMBER14+GLYCAM OpenMM system through ParmEd and ACPYPE. It is not the
  `GAFF_CANDIDATE` implementation and does not add arbitrary ligand templates.
- OpenMM performs final residue-template graph matching. Open Babel, OpenFF,
  AmberTools, ParmEd, and ACPYPE are external-tool adapters, not policy owners.
- `xtb` and `obminimize` are post-OpenMM geometry refiners. They cannot replace
  a successful preceding OpenMM parameterization or emit production MD
  parameters.

## Side Effects

- `lig_params` writes a temporary PDB, temporary SDFs, and a persistent GAFF
  cache under `$DVBFIXER_LIG_CACHE` or `~/.cache/dvbfixer/lig_params/`.
- It may prepend the active environment's `bin` directory to `PATH` and
  register `AmberToolsToolkitWrapper` globally.
- Non-strict dependency, extraction, or generator failures warn and may return
  `None`; complex cofactors raise `LigandParamError` before dependency checks.
- Minimize template failure can trigger strip-and-splice fallback, leaving the
  protein-ligand interface unrelaxed and warning that coordinates may be
  strained.
- Standalone `parametrize` and ACPYPE export create external-tool
  intermediates and final topology artifacts; their cleanup policies are local
  to those adapters.

## Known Divergences

- The five routes are not five implemented backend branches. Only
  `UNSUPPORTED_COMPLEX_COFACTOR` is directly consumed from this policy by a
  production path.
- Native-template screening compares residue names before OpenMM's later graph
  match; incomplete or differently bonded residues can still fail.
- The complex-cofactor list is a small name denylist, not connected-component,
  metal, attachment, redox, or completeness analysis.
- Unknown residues are deduplicated by residue name, and ligand insertion codes
  are unsupported by the Open Babel extraction adapter.
- Automatic GAFF failure can degrade to strip-and-splice; explicit GAFF failure
  raises instead.
- Standalone `parametrize` has PDB-only structural preflight. Other accepted
  formats do not receive the same residue-count and metal checks there.
- `top --acpype` has no `--extra-ff` or GAFF-candidate routing and assumes its
  fixed AMBER14+GLYCAM templates can describe the complete input.
- Policy tests exercise HEM/FAD rejection, a GAFF default, and user-template
  precedence, but do not cover every precedence branch or case normalization.
- Direct unit tests for `build_ligand_generator` strictness, cache behavior,
  and cofactor failure through minimize are absent. ACPYPE tests cover only its
  pure option parsers, not the external export pipeline.

## Proposed Work

Introduce an application service that inventories connected molecular
components, records evidence for exact native matches, and returns typed route
decisions consumed consistently by minimize, standalone `parametrize`, and
topology exporters. It must validate component completeness, covalent
attachment, metals, and charge/redox assumptions before invoking adapters.

Do not present geometry regularization as a shipped parameterization backend;
the whole-complex design remains research described in
[`../../research/whole-complex-relaxation.md`](../../research/whole-complex-relaxation.md).

## Focused Verification

```bash
pytest -q tests/test_scientific_domain.py tests/test_parametrize_preflight.py \
  tests/test_top_acpype.py
pytest -q tests/test_zbs_e2e.py -k ligand_no_severe_clash_without_parametrize_flag
python scripts/check_agent_docs.py
```

The ZBS regression is slow and requires the scientific and AmberTools stack.
