# Topology Generation

Status: partial

Verified on: 2026-09-18

Verified at commit: `425f290eb85760246766f1d1500e51672c640b2d`

## Purpose

This context owns GROMACS topology generation from a prepared PDB or GRO
structure. It covers the default RTP builder and the opt-in ACPYPE exporter,
including their different force-field, naming, and dependency contracts.

It does not own structure preparation, protonation-state selection, arbitrary
ligand chemistry, or the standalone GAFF2 `parametrize` command.

## Scope

`src/dvbfixer/top/pipeline.py::main` orchestrates input normalization, chain
classification, topology construction, molecule counting, and output. The
in-memory builders live in `top/topology_builder.py`; file emission lives in
`top/writers.py`; glycan graph handling lives in `top/glycan.py`.

The default route parses bundled or user-supplied GROMACS force-field files.
`--acpype` bypasses that route and delegates the complete build to
`top/acpype.py::run_acpype_mode` and `acpype_export.py::export_gromacs`.

## Capabilities

| Capability | Status | Owner | Evidence |
|---|---|---|---|
| RTP AMBER/CHARMM protein topology | implemented | `TopologyBuilder.build_chain` | `tests/test_top_pipeline_baseline.py` |
| CHARMM glycan and glycolipid topology | implemented | `TopologyBuilder.build_glycan_chain`, `build_glycolipid_chain` | glycoprotein baseline |
| Modular `.top`, `.itp`, restraints, and topology-matched PDB output | implemented | `top/writers.py` | protein and glycoprotein baselines |
| Chain merge preserving coordinates and identity | implemented | `top/pipeline.py::_merge_chains` | `tests/test_top_merge_chains.py` |
| AMBER14+GLYCAM ACPYPE export with mixed 1-4 scaling | implemented | `acpype_export.py::export_gromacs` | parser unit tests; external integration coverage is limited |
| Arbitrary unknown-ligand topology inside `top` | missing | none | unknown chains warn and are dropped |
| GAFF2 templates for OpenMM minimize | implemented, separate | `lig_params.py::build_ligand_generator` | not a `top` backend |

## Entry Points

- Public command: `src/dvbfixer/top/pipeline.py::main`.
- CLI and bundled FF resolution: `src/dvbfixer/top/cli.py::parse_args` and
  `bundled_ff_root`.
- RTP construction: `src/dvbfixer/top/topology_builder.py::TopologyBuilder`.
- GROMACS parser boundary: `src/dvbfixer/rtp_parser.py::parse_rtp`, `parse_r2b`,
  `parse_arn`, `parse_tdb`, and `parse_atomtypes`.
- ACPYPE route: `src/dvbfixer/top/acpype.py::run_acpype_mode`.
- Shared exporter: `src/dvbfixer/acpype_export.py::export_gromacs`.

## Contracts

The RTP route consumes a residue/atom graph whose names can be resolved to a
loaded RTP building block. R2B chooses main and terminal blocks, ARN translates
atom spellings, TDB supplies CHARMM terminal patches, and ATP supplies masses.
`TopologyBuilder` then derives bonds, angles, proper dihedrals, 1-4 pairs,
impropers, and CMAP entries into `ChainTopology`.

Atom matching is name-driven. Supported ARN, explicit aliases, and hydrogen
number shifts bridge common PDB/OpenMM and GROMACS conventions; matching does
not prove that an arbitrary residue is chemically equivalent to a template.
The output PDB must use the same selected atoms and names as the emitted
topology, while preserving coordinates and original chain/residue identity.

The ACPYPE route is fixed to OpenMM AMBER14 + GLYCAM + TIP3P parameters. It
creates an unconstrained OpenMM `System`, converts it through ParmEd to AMBER
`prmtop`/`inpcrd`, then asks ACPYPE for GROMACS files. ACPYPE's `[ pairs_nb ]`
output is required because AMBER and GLYCAM use different 1-4 scaling factors.

GAFF2 is not selected by `top --acpype`. `lig_params.py` generates OpenMM GAFF
templates only for isolated unknown organic residues used by minimize/ZBS.
For a GROMACS ligand, parameterize it separately and include its validated ITP;
do not make the RTP classifier pretend that a name is a native template.

## Invariants

- Keep `PDBFile.writeFile(..., keepIds=True)` in topology-adjacent OpenMM paths.
- Chain IDs are case-sensitive; `D` and `d` are distinct.
- RTP atom indices and every bonded term must remain contiguous and consistent
  after terminal patches or chain merging.
- Unknown chains must emit a visible warning before being omitted; they must
  never disappear silently from `[ molecules ]`.
- `add_glycam_bonds` must directly import and use `nanometer` and the shared
  sugar detector; import failures must not be swallowed.
- ACPYPE preprocessing must preserve explicit variant intent across OpenMM's
  residue-name normalization and repair hydrogens/chirality after addition.
- Glycosidic linkage hydrogens are removed with matching charge/type changes
  unless `--keep-all-hydrogens` explicitly opts into the risky raw geometry.
- Interchain interaction includes remain after `[ molecules ]` in `topol.top`.

## Callers

- `dvbfixer top` dispatches directly to `top.pipeline.main`.
- `structure_input.py` normalizes CIF input to temporary PDB before dispatch.
- `top.pipeline.main` calls CONECT materialization before either topology route
  unless `--no-infer-conect` is set.
- `transplant --gromacs` also uses the shared ACPYPE exporter.
- External GROMACS tools consume `topol.top`, included ITP files, and the
  topology-matched PDB or ACPYPE GRO output.

## Adapters

- `rtp_parser.py` is the adapter for GROMACS RTP/R2B/ARN/TDB/ATP syntax.
- Bundled `FF/*.ff` trees are build-time data; RTP output is self-contained and
  does not require a target-machine force-field installation.
- MDAnalysis converts GRO input to a temporary PDB.
- OpenMM supplies template matching, hydrogen placement, and AMBER+GLYCAM
  parameterization on the ACPYPE route.
- ParmEd bridges OpenMM systems to AMBER files; ACPYPE bridges AMBER to GROMACS.
- Open Babel, OpenFF, AmberTools, and `GAFFTemplateGenerator` belong to the
  separate `lig_params.py` adapter chain, not RTP topology generation.

## Side Effects

The RTP route writes `topol.top`, `ffparams.itp`, one ITP and position-restraint
file per built molecule type, `water.itp`, `ions.itp`, optional
`interchain_ss.itp`, and a topology-matched PDB. Outputs are written directly,
not transactionally; a late failure can leave partial artifacts.

The ACPYPE route writes `_gmx_temp.pdb`, temporary AMBER files, and an
`*.amb2gmx` directory beside the input, temporarily changes the process working
directory, copies final files into the output directory, then removes known
temporaries. External tool exceptions propagate.

`lig_params.py` separately creates temporary PDB/SDF files and a persistent
GAFF cache under `~/.cache/dvbfixer/lig_params/` or `$DVBFIXER_LIG_CACHE`.

## Known Divergences

- RTP mode warns and drops a wholly unrecognized chain, then succeeds if other
  recognized molecules remain. It exits when nothing recognized can be built.
- A recognized chain can still fail atom/template resolution; builders warn and
  return `None`, and orchestration may continue with other chains.
- Residue maps in topology code commonly key by `(chain, resseq)` and therefore
  do not consistently preserve insertion-code identity.
- Missing RTP atoms may be skipped and their charge moved to a bonded retained
  atom; this permissive behavior is not a general completeness validator.
- ACPYPE ignores RTP `--ff`, `--water`, `--ignh`, and `--merge`; selecting
  CHARMM with `--acpype` warns and still uses AMBER14+GLYCAM.
- Malformed or unknown ACPYPE protonation states can be silently omitted after
  parsing, unlike the stricter validation in the RTP orchestration path.
- Full ACPYPE export depends on external scientific packages and has less
  automated end-to-end coverage than the RTP route.

## Proposed Work

Add structured topology-build results that report every retained, omitted, and
failed component before committing output. Route unknown components through the
shared parameterization policy without silently treating GAFF candidates as RTP
templates, and preserve `(chain, resseq, icode)` throughout topology maps.

Add direct parser tests for RTP/R2B/ARN/TDB and an external-dependency-gated
ACPYPE export test that validates `[ pairs_nb ]`, cleanup, and final artifacts.

## Focused Verification

```bash
pytest -q tests/test_top_pipeline_baseline.py tests/test_top_merge_chains.py \
  tests/test_top_acpype.py tests/test_cli_stabilization.py \
  tests/test_prepare_propka_integration.py
python scripts/check_agent_docs.py
git diff --check -- docs/agent/contexts/topology-generation.md
```

See [`../../force-fields.md`](../../force-fields.md) for the two `--ff`
namespaces and [`../../commands/top.md`](../../commands/top.md) for user-facing
operation and output details.
