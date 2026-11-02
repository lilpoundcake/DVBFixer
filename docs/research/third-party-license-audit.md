# Third-party components and licensing audit

Initial audit: **2026-10-27**. Repository/metadata refresh: **2026-10-31**,
baseline **`4a328a6` / 0.10.1**. Status: **preliminary inventory, not release
clearance or legal advice**. The initial audit did not select a project license;
the owner subsequently chose BSD-3-Clause and confirmed Stanislav Skripkin.
Original-code LICENSE/metadata were implemented on 2026-11-02; third-party
clearance remains open.

The [2026-11-02 exact-source continuation](dependency-license-review-2026-11-02.md)
supersedes generic observations where actual selected artifacts differ: the
Open Babel wheel recipe pins a fork, ACPYPE includes a native Amber payload,
Reduce 4.15 has custom notices, and the CHARMM content reference now has a
byte-level comparison. The owner reports CHARMM acquisition from GROMACS;
its exact package/release is still needed to establish that acquisition chain.

[Implementation plan](../plans/licensing-and-attribution.md) ·
[Execution handoff](../plans/licensing-audit-handoff.md) ·
[Existing notices](../../THIRD_PARTY_NOTICES/README.md) ·
[Dependency/environment audit](../plans/environment-dependency-audit.md)

## Scope and evidence limits

The inventory covers declared Python/GUI dependencies, Python imports and external
executables, deployment manifests, model identities, bundled force fields and
test-data provenance. Evidence comes from `pyproject.toml`, `environment.yml`,
`gui/package.json`, `gui/package-lock.json`, `src/dvbfixer/`, `deploy/`, and
upstream license files/model records inspected during the audit.

The 2026-10-31 refresh additionally checked the canonical
`src/dvbfixer/data/installation-profiles.toml`, runtime provisioning, packaging
declarations and current GUI lockfile. Original upstream findings below are not
all reverified by this refresh; new PyPI observations are explicitly metadata
only. No full environment export, binary inspection or transitive SBOM was made.

The subsequent implementation collected eight installed Apple metadata/legal
inventories and inspected actual wheel/sdist snapshots. See the
[continuation findings and owner gate](licensing-decisions-2026-10-31.md) and
[portable evidence index](licensing-evidence-2026-10-31.json). Earlier evidence
limits below describe initial discovery unless explicitly superseded there;
cross-platform closure and complete compatibility review remain open.

The [generated declaration inventory](license-inventory.json) is the first WP1
implementation slice. Run `python3 scripts/gen_license_inventory.py --check`
to detect changes to declarations, selected tracked data and review-input hashes.
It currently contains 804 records (uses/assets/lock entries, not unique software
projects): 20 Python runtime requirements, five build/dev requirements, 40
catalog entries, five external integrations, 32 direct GUI declarations, 558
npm lock entries, 120 declared packaged files and 24 pending review inputs.
Untracked material is excluded from this register, not guaranteed absent from a
local build: setuptools globs and actual archives still need inspection. All
compatibility decisions remain `unknown`; npm license fields are `declared`
metadata only. Existing upstream observations in this prose have not been
promoted automatically into exact-version, artifact-verified evidence records.

An upstream license observation is not proof that every permitted dependency
version, binary wheel, container or auxiliary dataset has identical terms.
Unpinned packages need release-specific verification. This is not a complete
transitive Python/system SBOM, a copied-code similarity scan or a compatibility
opinion. Distinguish source imports, subprocess tools, bundled material and
operator-provided artifacts; they have different distribution boundaries.

## Priority findings

1. At the audited baseline, no root project license or Python declaration existed.
   This original-code declaration gap was addressed on 2026-11-02 with the
   owner-approved BSD-3-Clause grant. It does not settle third-party rights.
2. Protpardelle-1c source is **MIT**, but its official weights archive is
   **CC-BY-4.0**. Do not label both MIT.
3. Open Babel and ACPYPE are directly imported GPL components, not exclusively
   independent CLI programs. Review the combined distribution before selecting
   the project's license. Merely making an import optional does not establish
   compatibility.
4. MODELLER has its own restrictive terms. Academic access does not establish
   permission to redistribute it or operate every commercial/service deployment.
5. Existing third-party notices cover Amber19SB/GROMACS only, while packaging
   also distributes two other force-field trees with incomplete license evidence.

## Implemented installation boundaries (2026-10-31)

Standard and Desktop are implemented, not future dependency-slimming proposals.
Standard includes PySCF/RESP, MDAnalysis/Plotly, MODELLER, AmberTools, Reduce,
MAFFT and ANARCI/HMMER. Desktop inherits Standard and adds Node/npm and GUI.
`openff-toolkit-base` is the conda route; openmmforcefields, ACPYPE and
openbabel-wheel are pip packages. Inspect their actual installed distributions,
not the license of a similarly named conda metapackage.

The `msa-alternatives` add-on adds MUSCLE 5 and Clustal Omega. Probe 2.18 is an
optional add-on (managed Linux provisioning; native Apple can be operator
provided). Its experimental 2.26 candidate is not the production dependency.
`ambertools-compat` changes Python/NumPy/Biopython bounds and adds `pdb2pqr`;
it needs its own resolved inventory. PSI4/psiresp is separate; Gaussian is
operator-provided. xTB is removed from supported profiles.

The helper provisions pinned, locally patched Protenix v1 and cc94 CUDA/MPS
sources and their runtime environments. Review fetched source/data terms and
modification notices as well as any future redistributed image or archive.
Installation or checkpoint-backed scientific acceptance is not license clearance.
Historical cc89 and research Boltz/RFdiffusion paths remain separate audit scopes.

The source/wheel, resolved installed environments, browser assets, standalone
server, helper-fetched artifacts, research/build tools and operator-provided
software must not be collapsed into one undifferentiated dependency list.

## Model source, weights and auxiliary artifacts

| Model | Repository use | Source terms | Weight terms and evidence |
|---|---|---|---|
| [Protpardelle-1c](https://github.com/ProteinDesignLab/protpardelle-1c), cc94 | New explicit production requests use `cc94_epoch3100.pth` | [MIT at pinned commit](https://github.com/ProteinDesignLab/protpardelle-1c/blob/ee378400f25b801fa481028000f9060183d7fb4c/LICENSE) | **CC-BY-4.0** in the [official Zenodo model record](https://zenodo.org/records/16817230), also exposed by its [record API](https://zenodo.org/api/records/16817230) |
| Protpardelle-1c, cc89 | Frozen historical-request replay, `cc89_epoch415.pth` | MIT | CC-BY-4.0 for weights obtained from the same official archive |
| [Protenix v1](https://github.com/bytedance/Protenix) | Optional integration, `protenix_base_default_v1.0.0.pt` | [Apache-2.0 at pinned commit](https://github.com/bytedance/Protenix/blob/85767b811c40ed46e73a9b39519cf6bfca8701ba/LICENSE) | The [pinned README](https://github.com/bytedance/Protenix/blob/85767b811c40ed46e73a9b39519cf6bfca8701ba/README.md#license) explicitly includes model parameters under Apache-2.0 |
| [Boltz-2](https://github.com/jwohlwend/boltz) | Research integration, not a public DVBFixer backend; gap smoke uses `boltz2_conf.ckpt` | [MIT at pinned commit](https://github.com/jwohlwend/boltz/blob/b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc/LICENSE) | The [pinned README](https://github.com/jwohlwend/boltz/blob/b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc/README.md) explicitly includes weights; the [pinned Hugging Face card](https://huggingface.co/boltz-community/boltz-2/blob/6fdef46d763fee7fbb83ca5501ccceff43b85607/README.md) also declares MIT |
| [RFdiffusion v1](https://github.com/RosettaCommons/RFdiffusion) | Internal research adapter, `Base_ckpt.pt` | [BSD-3-Clause at pinned commit](https://github.com/RosettaCommons/RFdiffusion/blob/bf42b54c20a99dd7350456c85985ed4d83b95d48/LICENSE) | [Current upstream LICENSE](https://github.com/RosettaCommons/RFdiffusion/blob/main/LICENSE) explicitly includes README-linked weights; that clarification postdates the pinned source. Preserve an immutable copy/reference before release approval |

### Identities and repository evidence

| Artifact | SHA-256 | Local authority |
|---|---|---|
| `cc94_epoch3100.pth` | `a58c009ead22b6f79ef61d44e258b432b0fada1dfd8d44cb7a137d7459038a15` | [cc94 deployment guide](../../deploy/protpardelle-1c/cc94.md) |
| `cc89_epoch415.pth` | `dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483` | [Protpardelle deployment guide](../../deploy/protpardelle-1c/README.md) |
| `protenix_base_default_v1.0.0.pt` | `2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04` | [Protenix profile lock](../../deploy/protenix-v1/profile-lock.json) |
| `boltz2_conf.ckpt` | `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1` | [Boltz research guide](../../deploy/boltz-2/README.md) |
| `Base_ckpt.pt` | `0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca` | [RFdiffusion deployment guide](../../deploy/rfdiffusion-v1/README.md) |

- No tracked `.pt`, `.pth`, `.ckpt`, `.onnx` or `.safetensors` weights were found.
  Operator-provided weights still need documented terms; embedding them in an
  image is a separate distribution decision.
- The initial inspected Protpardelle paths omitted checkpoint-license provenance.
  The 2026-10-31 continuation corrects the shared adapter to populate MIT source
  and CC-BY-4.0 checkpoint fields separately, without changing sampling. Existing
  frozen requests/results are not rewritten.
- The pinned Protenix README explicitly excludes **v2 weights** from open-source
  licensing. Do not generalize the v1 finding to v2 or future weights.
- Boltz's deployment guide also records `boltz2_aff.ckpt` and `mols.tar`.
  A recorded cache artifact is not proof that every smoke uses it. Audit the
  molecular archive separately rather than inheriting a checkpoint license.
- Protenix's lock records CCD, derived RDKit molecules, cluster and release-date
  tables. Verify each auxiliary artifact's source and terms.
- ProteinMPNN/LigandMPNN sampling is disabled in the inspected Protpardelle
  adapter (`skip_mpnn_proportion=1.0`, full-MPNN flags false); ESMFold execution
  was not established in that path. Upstream installation suggestions alone
  do not prove these weights are DVBFixer runtime dependencies.
- Local callback/Apple portability patches must retain relevant upstream notices
  and document modifications if patched engine source is distributed.

## Scientific Python dependencies

Links identify upstream projects or package metadata. Exact release notices,
bundled code and binary-library obligations remain part of release validation.

| Component | Role | Observed license / caveat |
|---|---|---|
| [NumPy](https://github.com/numpy/numpy) | Numerical arrays | BSD-3-Clause core; binary/vendored components have additional notices; keep the project's `<2.5` constraint |
| [SciPy](https://github.com/scipy/scipy) | Geometry, optimization, spatial search | BSD-3-Clause core plus bundled-component notices |
| [Biopython](https://github.com/biopython/biopython/blob/biopython-188/LICENSE.rst) | Sequence and structure operations | Biopython License Agreement; some files additionally BSD-3-Clause, not a blanket dual-license claim |
| [Gemmi](https://github.com/project-gemmi/gemmi/blob/v0.7.5/LICENSE.txt) | CIF/mmCIF and CCD | MPL-2.0, checked at `0.7.5` |
| [OpenMM](https://github.com/openmm/openmm/blob/master/docs-source/usersguide/library/01_introduction.rst) | Molecular mechanics, hydrogens, refinement | MIT core/CPU/application layer; GPU platforms LGPL; inspect release license inventory |
| [PDBFixer](https://github.com/openmm/pdbfixer) | Structure repair | MIT |
| [PROPKA](https://github.com/jensengroup/propka/blob/master/LICENSE) | Protonation | LGPL-2.1 |
| [MDAnalysis](https://github.com/MDAnalysis/mdanalysis/blob/release-2.10.0/LICENSE) | Structure analysis and clustering | LGPL-3.0-or-later for checked `2.10.0`; do not apply this retrospectively to every older version |
| [openmmforcefields](https://github.com/openmm/openmmforcefields) | Parameter generators and force-field integration | MIT code; parameter/data terms need separate review |
| [Open FF Toolkit](https://github.com/openforcefield/openff-toolkit) | Ligand parameterization | MIT code; separately inventory parameter datasets |
| [RDKit](https://github.com/rdkit/rdkit) | Molecular graphs and hydrogens | BSD-3-Clause core; check packaged notices |
| [Open Babel](https://github.com/openbabel/openbabel/blob/openbabel-3-1-1/README.md) / [openbabel-wheel](https://pypi.org/project/openbabel-wheel/) | Chemistry, conversion and refinement | GPL v2; inspect exact build and file-level terms before compatibility approval |
| [ParmEd](https://github.com/ParmEd/ParmEd) | Parameter/topology export | Package metadata says LGPL; exact edition remains to be verified for the selected release |
| [ACPYPE](https://github.com/alanwilter/acpype) | GROMACS export | Current package metadata declares GPL-3.0-or-later; verify the release actually distributed |
| [Plotly.py](https://github.com/plotly/plotly.py) | Plots | MIT |
| [ANARCI](https://github.com/oxpig/ANARCI/blob/master/LICENCE) | Antibody numbering | BSD-3-Clause code; inventory HMM profiles and their underlying data separately |
| [PySCF](https://github.com/pyscf/pyscf) | Standard, required Python dependency for RESP QM | Apache-2.0 core; PyPI `2.14.0` metadata rechecked on 2026-10-31; inspect bundled/native libraries separately |
| [TinyDB](https://pypi.org/pypi/tinydb/4.9.0/json) | Force-field template cache | MIT in checked PyPI metadata |
| [lxml](https://pypi.org/pypi/lxml/6.1.3/json) | Force-field XML | BSD-3-Clause in checked PyPI metadata; inspect libxml2/libxslt and wheel notices separately |
| [PyYAML](https://pypi.org/pypi/pyyaml/6.0.3/json) | Force-field configuration | MIT in checked PyPI metadata; inspect bundled LibYAML if present |
| [validators](https://pypi.org/pypi/validators/0.36.0/json) | Force-field source validation | MIT in checked PyPI metadata |
| [threadpoolctl](https://pypi.org/pypi/threadpoolctl/3.7.0/json) | BLAS runtime verification | BSD-3-Clause in checked PyPI metadata; does not determine the loaded BLAS license |

These tables now name all **20** direct requirements in `pyproject.toml`, plus
selected catalog/external components; this is enumeration, not 20 completed
license reviews. The added versioned PyPI records were read on 2026-10-31 and
are not installed-version evidence or archived/hash-bound legal texts. ParmEd
`4.3.1` PyPI metadata still reports only `LGPL`; the subsequently captured conda
build declares LGPL-2.1-or-later and includes LGPL 2.1 text. That exact-build
finding does not clear every version/file or the combined distribution.

The supported catalog additionally installs Python, pip, OpenBLAS-selected
BLAS/CBLAS/LAPACK and (Desktop) Node/npm. Review the selected builds, native
runtime/compiler libraries and bundled dependencies. The
[`pdb2pqr` 3.8.0 metadata](https://pypi.org/pypi/pdb2pqr/3.8.0/json) contains
BSD-style three-clause terms, but this does not clear all data or transitive
dependencies of the compatibility add-on. Build requirements (setuptools),
Python dev extras and the test stack also require a separate build/source
inventory; they are not automatically shipped runtime components.

Open Babel imports occur in chemistry/preparation/parameterization paths;
`src/dvbfixer/acpype_export.py` imports `acpype.topol.MolTopol`. Do not claim that
all GPL dependencies are isolated subprocess tools. Likewise, a GPL dependency
does not by itself settle which license can be applied to all original source;
review the actual combined work and distribution, including GPL-version issues.

## External scientific programs

| Tool | Role | Observed terms / remaining check |
|---|---|---|
| [MODELLER](https://salilab.org/modeller/registration.html) | Default gap/comparative modeling, optional SALIGN engine | Public EULA captured in continuation: internal non-commercial use, non-transferable grant, written permission required for distribution; owner-specific agreement/deployment rights remain open |
| [AmberTools](https://ambermd.org/AmberTools.php) | `tleap`, `antechamber`, `parmchk2`, associated calculations | Mixed suite, predominantly GPL; identify the actual release's component licenses, including indirect executables |
| [Reduce](https://github.com/rlabduke/reduce) | Hydrogen placement and tautomer/orientation decisions | Selected bioconda 4.15 build has UCSF/Word/ReliaSolve terms; metadata says BSD-4-Clause-UC. Earlier generic Apache-2.0 observation must not be applied to this artifact; see exact-source continuation |
| [Probe](https://github.com/rlabduke/probe/blob/master/LICENSE.txt) | Contact diagnostics | Apache-2.0 in inspected upstream |
| [xTB](https://github.com/grimme-lab/xtb) | Historical research only; removed from current installation profiles | LGPL-3.0 family; audit only if a historical/research artifact actually includes it |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | MSA | Standard source without extensions is BSD; Windows/extensions packages can have different terms |
| [MUSCLE 5](https://github.com/rcedgar/muscle) | MSA | GPL-3.0 family; use v5, not v3 |
| [Clustal Omega](http://www.clustal.org/omega/) | MSA | GPL-2.0-or-later in package recipe |
| [HMMER](https://github.com/EddyRivasLab/hmmer) | ANARCI support | BSD; include Easel/component notices where distributed |
| [Kalign 3.6.0](https://github.com/TimoLassmann/kalign/blob/v3.6.0/COPYING) | Protenix environment | Apache-2.0 |
| [Psi4](https://github.com/psi4/psi4) / [psiresp](https://pypi.org/project/psiresp/) | Alternative isolated RESP backend | LGPL-3.0 / LGPLv3 metadata; inspect actual releases |
| [Gaussian](https://gaussian.com/) | Input preparation and supplied-log consumption | Proprietary external program, not distributed by this path |

`doctor` detecting `gmx` does not establish that native GROMACS is required by
the standalone `top` generator. Generating a topology and distributing GROMACS
itself are different boundaries.

## Bundled force fields and structural data

`pyproject.toml` explicitly installs three FF trees and `THIRD_PARTY_NOTICES/*`:

| Material | Existing evidence | Missing work |
|---|---|---|
| `FF/amber19sb.ff/` | [Pinned GROMACS source/provenance and checksums](../provenance/amber19sb-gromacs.md); LGPL-2.1 text and attribution already shipped | Preserve and validate notices in built artifacts |
| `FF/amber99sb-ildn-lipid21.ff/` | `forcefield.doc` cites AMBER99SB-ILDN/Lipid21; `forcefield.itp` credits the Sorin ffamber port | Recover exact upstream archive/version, Lipid21 port origin, modifications and applicable permissions for all portions |
| `FF/charmm36_ljpme-jul2022.ff/` | Owner reports GROMACS acquisition; [official-archive content comparison](../provenance/charmm36-jul2022.md) finds 59 identical files, three BUF-modified files and omitted `watermodels.dat` | Identify actual GROMACS package/release, BUF modification provenance and parameter/port redistribution terms |
| `src/dvbfixer/data/charmm-terminal-caps.xml` | Explicit-residue equivalents of CHARMM36 ACE/CT3 patches | Document parameter source and applicable attribution/terms |

Do not infer that a CHARMM parameter archive is LGPL merely because it uses a
GROMACS format, or that a converter's license applies to all converted data.
Do not assign one license to every file under `FF/`: additional local FF trees
were untracked at audit time and were not approved distribution inputs.

[Force-field selection documentation](../force-fields.md) also identifies the
OpenMM AMBER ff19SB/ff14SB, GLYCAM, Lipid17, DNA OL15/RNA OL3, CHARMM36, water and
ion parameter sets. GAFF/GAFF2 comes through the parameterization stack. Record
parameter provenance independently from Python wrapper licenses.

`src/dvbfixer/ccd.py` downloads CCD records from RCSB. The [RCSB usage policy](https://www.rcsb.org/pages/policies)
states PDB archive data are CC0 and recommends structure-author citations.
The [fixture inventory](../../tests/fixtures/README.md) contains deposited
structures, user-provided derivatives, CHARMM-GUI material and an AmberClassic
example. Do not label every fixture CC0 without checking its actual origin and
any additional contributions. Preserve accession IDs, source URLs and hashes.

## GUI, build tools and model runtimes

The inspected `gui/package-lock.json` contained **558 non-root package entries**,
including development dependencies and platform variants. All declared a
license: 499 MIT, 21 Apache-2.0, 21 ISC, six BSD-2-Clause, six BSD-3-Clause, and
one each BlueOak-1.0.0, Python-2.0, CC-BY-4.0, 0BSD and `(MIT OR CC0-1.0)`.
This is a lockfile observation, not a count of packages shipped in the browser.
The count and license histogram were recomputed unchanged on 2026-10-31;
the root GUI package still has no declared license. `private: true` prevents
accidental npm publication but does not grant usage or redistribution rights.

Direct components include [React](https://github.com/facebook/react),
[Mol*](https://github.com/molstar/molstar), [MUI](https://github.com/mui/material-ui),
Emotion, Zustand, TypeBox, openapi-fetch, Connect, serve-static, Sass and Rollup
WASM (MIT), and [FlexLayout](https://github.com/caplin/FlexLayout) (ISC).
Inspect embedded assets and component notices as well as npm metadata.
`build:server` uses esbuild with `--packages=external`: inspect installed server
dependencies as well as `dist-server/server.js`. Browser and server inventories
must be separate, and source maps/assets must be included when distributed.
`caniuse-lite` is CC-BY-4.0; assess actual build/distribution inclusion rather
than assuming every development dependency ships to users.

Deployment environments also introduce PyTorch, Lightning, Biotite, Hydra,
OmegaConf, Transformers, Hugging Face Hub, e3nn, DGL, ProDy and other packages.
They are not all core dependencies or all fully license-audited here. An image
SBOM must include its exact Python/system libraries, BLAS, SE(3)-Transformer,
NVIDIA dependencies and CUDA/cuDNN redistribution terms. Development tooling
(pytest, ruff, mypy, Vite, TypeScript, ESLint, Vitest, etc.) belongs in the
source/build inventory, not automatically in every runtime notice bundle.

## Attribution, citations and unresolved release gates

Legal notices and scientific citations serve different purposes. A project URL
or paper citation is not a substitute for a required license/copyright notice,
Apache NOTICE, change notice, source availability or copyleft obligation.

Bibliography work should cover the actually used tools/models/parameter sets.
In particular, upstream Protpardelle requests both
[Chu et al.](https://doi.org/10.1073/pnas.2311500121) and
[Lu & Shuai et al.](https://doi.org/10.1101/2025.08.18.670959), Boltz-2 has its
[technical report](https://doi.org/10.1101/2025.06.14.659707), and Mol* has its
[viewer paper](https://doi.org/10.1093/nar/gkab314). Preserve the parameter
citations already present in `forcefield.doc`; complete and validate the
bibliography rather than treating these examples as exhaustive.

Release gates remain open for combined-distribution compatibility, MODELLER terms,
the two older bundled FF trees, auxiliary model/HMM data, contributed fixture
rights and exact transitive/binary notices. The linked implementation plan
tracks deliverables and owner decisions; this report records findings only.
