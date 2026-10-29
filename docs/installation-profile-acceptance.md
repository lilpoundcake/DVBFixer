# Installation profile acceptance

## Readiness contract

Readiness is scoped to the selected profile and add-ons:

- **ready**: all selected checks executed successfully;
- **action-required**: a selected licensed/operator component needs configuration
  (MODELLER license, checkpoints, Gaussian, or a separate runtime);
- **error**: installation/process failure, missing required dependency/data,
  incompatible native ABI, or unsupported platform;
- **not-tested**: acceptance evidence has not been executed on that manager/host.
  This is an evidence state, never equivalent to ready.

A successful resolver or import is not full scientific acceptance. The catalog's
module/executable checks are prerequisites; WP5 also requires real workflows.
Optional Probe absence uses the existing Python diagnostic fallback. Diffusion
runtime readiness and scientific/hardware acceptance are separate gates.

## Baseline inspection (2026-10-28, Linux x86_64)

- Checkout: `feature/installation-profiles`, initial HEAD `cc0716c`.
- Existing environment left unchanged. Its default editable import resolves to
  another worktree (`diffusion-cli-integration`, version 0.9.0); development
  checks explicitly set `PYTHONPATH` to this checkout's `src` (version 0.10.1).
- Existing ACPYPE 2026.9.4 PyPI metadata requires
  `openbabel-wheel>=3.1.1.23`, `rich>=14.2.0`, `typer>=0.20.0`.
- `micromamba repoquery depends acpype -c conda-forge -c bioconda --json`
  reports conda-forge ACPYPE 2026.10.15 `pyhd8ed1ab_0`, SHA256
  `2f84f602eb3eadac11701da6c92f70d03a90b52a1ccad67e33f3175c80677573`.
  Conda requirements are AmberTools, Open Babel, Python >=3.12, Rich, Typer.
- The conda-forge recipe uses the PyPI sdist with no dependency-metadata patch;
  PyPI 2026.10.15 still declares `openbabel-wheel`. Therefore conda solver
  metadata alone does **not** yet prove a compatible single-provider route.
  Installed conda ACPYPE was inspected: its Python metadata still requires the
  wheel. That route was rejected. The selected route uses pip ACPYPE and
  `openbabel-wheel` together, with no conda Open Babel and no `--no-deps` bypass.
- The first actual solve also installed NAGL/PyTorch through the full OpenFF
  metapackage. That experimental prefix was rejected and retained for evidence.
  `openff-toolkit-base` provides the used Molecule/AmberTools APIs without NAGL.
  PyPI `openmmforcefields==0.15.1` avoids reintroducing the conda metapackage;
  its undeclared TinyDB/lxml/PyYAML/validators requirements are explicit in the
  catalog and Python metadata. GAFF template generation was actually exercised.

## Composition

The [dependency matrix](installation-dependency-matrix.md) and files under
`environments/` are generated using `scripts/gen_installation_profiles.py`.
The TOML catalog is the one composition source; manager commands belong in
adapters. After Linux replacement acceptance, `environment.yml` is the generated
Standard alias; `environments/legacy-full.yml` preserves Desktop plus alternative
MSA as an explicit compatibility input. Native Apple acceptance remains open.

Diffusion provisioning instructions and phase boundaries are in
[separate runtime provisioning](installation-runtime-provisioning.md).

## Implemented scope and next gate

- **WP1 complete:** composition/readiness contract, dependency matrix, actual
  ACPYPE metadata and single-provider investigation.
- **WP2 complete:** generated Standard/Desktop, MSA alternatives, optional Probe
  and separate PSI4 inputs; Python metadata and bundled catalog; real Python
  3.11/3.13 installs and native solves. Existing diffusion deployment inputs are
  referenced, not redefined as scientific profiles.
- **WP3 implemented (native MPS acceptance not-tested):** interactive/noninteractive helper, three manager adapters,
  new-prefix ownership records, resume/reconciliation, protected MODELLER setup
  and Desktop lockfile build/standalone smoke work. The helper now provisions
  pinned diffusion sources and isolated runtimes through `runtime_install.py`,
  stages registration using `install_diffusion_backend.py`, and publishes only
  after canonical runner preflight. The lower-level script remains registration-only.
  The helper now runs the real scientific smoke script after installation and
  persists a unique evidence bundle below the target prefix. Configure-only
  completion runs the licensed modeling/Homology steps without rebuilding GUI.
- **WP4 implemented:** selected-profile doctor and scientific-command preflights
  distinguish unavailable dependencies, calculation failures and explicit
  disabling; ANARCI/HMMER/data failures do not become negative domain findings;
  MUSCLE 5 and topology/trajectory/TPR checks are early and actionable. Included
  PySCF and separate PSI4 actual RESP paths were exercised.
- **WP5 Linux acceptance executed:** final main-profile and CUDA/QM manager
  routes, slow structural regressions, extended format/ligand/topology/analysis
  checks and native GROMACS preprocessing passed. Native macOS Standard/Desktop,
  optional native dependencies and cc94 MPS remain **not-tested**.
- **WP6 Linux candidate migrated:** default/manual/generated specs, profile CI,
  GUI/MSA separation, user/agent docs and the protected-license CI readiness
  path are synchronized. Whole-plan completion and native platform claims await
  the [Apple acceptance handoff](plans/installation-profiles-apple-acceptance.md).

## Executed Linux matrix

All target prefixes were newly created below `/tmp/opencode/`; existing named
environments were preserved. Mamba and conda were bootstrapped in a separate
tools prefix. No manager was installed into or substituted for the user's
existing scientific environment.

| Manager | Standard scientific smoke | Desktop science | Desktop build + standalone | MODELLER import/model/Homology |
|---|---|---|---|---|
| micromamba 2.9.0 | 12/12 passed | 12/12 passed | passed | passed in both profiles |
| mamba 2.10.0 | 12/12 passed | 12/12 passed | passed | passed in both profiles |
| conda 26.9.4 | 12/12 passed | 12/12 passed | passed | passed in both profiles |

The 12 scientific steps are RDKit/Open Babel SDF/PDB + OpenFF/GAFF2 system
creation with finite Reference energy, actual HF/6-31G* two-stage PySCF RESP
fit on water, ANARCI IMGT numbering using HMMER/database, MDAnalysis DCD
alignment and Plotly HTML, legacy preparation, opt-in tleap-reduce preparation,
the wheel's obminimize UFF executable, protonation, diagnose, a 10-iteration
Reference `minimize --no-solvent`, protein topology generation, and MAFFT
alignment. Newer smoke reports also assert nonempty prepared/minimized PDBs
and the all-L invariant. MODELLER smoke deletes a real internal residue from
tracked 1UBQ, reconstructs it using FASTA, and runs a separate Homology model.

Representative resolved versions/builds (full exports are retained in each
local `report.json`):

| Component | Python 3.13 Linux route | Python 3.11 baseline route |
|---|---|---|
| Python | 3.13.18 `hf47f18c_100_cp313` | 3.11.18 `h5f976f7_0_cpython` |
| NumPy | 2.4.6 | 2.4.6 |
| SciPy | 1.18.1 | 1.17.1 |
| OpenMM | 8.6.1 | 8.6.1 |
| OpenFF Toolkit Base | 0.19.0 | 0.18.0 `pyhd8ed1ab_6` |
| ACPYPE | 2026.10.15 | 2023.10.27 (pip's compatible Python selection) |
| openbabel-wheel | 3.1.1.23 (native OB reports 3.1.0) | 3.1.1.23 |
| PySCF | 2.14.0 | 2.14.0 |
| MODELLER | 10.8 | installed, licensed modeling not retested on 3.11 |
| AmberTools | 26.0 | 26.0 `cuda_None_nompi_py311h232f31e_100` |
| libblas | 3.11.0 `11_h4a7cf45_openblas` | same wrapper build |
| libcblas | 3.11.0 `11_h0358290_openblas` | same wrapper build |
| liblapack | 3.11.0 `11_h47877c9_openblas` | same wrapper build |

Both real NumPy/SciPy solves passed. `threadpoolctl` inspection **after** the
solves reports the environment's loaded `libopenblasp-r0.3.34.so`. AmberTools
still brings an MKL runtime dependency; that is not evidence that NumPy/SciPy
are executing MKL. Python 3.11 also passed the same 12 scientific smoke steps.

### Bugs exposed by actual installs

- The Open Babel wheel bundles `obminimize` but publishes only an `obabel`
  console entrypoint. DVBFixer now supplies a thin `obminimize` entrypoint that
  imports the wheel to initialize its plugin/data paths, then execs its bundled
  native binary. Actual UFF minimization on ethanol passed.
- `conda run` default capture disconnects stdin; protected MODELLER setup failed
  in the first conda attempt. The adapter uses `--no-capture-output`; the helper
  still captures/suppresses configuration/import output. Actual import, gap
  rebuilding and Homology then passed in Standard and Desktop.
- An interrupted native transaction can leave `conda-meta/history` before all
  packages are linked. Resume reconciles the native spec before pip; it must not
  silently replace missing native packages with wheels.
- **Whitespace in environment prefixes is rejected early.** A real Desktop
  install in a spaced prefix built and passed 11/12 scientific steps, but
  MAFFT's bundled shell script failed on its baked-in path. The failed prefix
  is retained, not counted as accepted. Choose a prefix without whitespace;
  adapter argv preserves spaces in source/spec/input paths.
- `pip check` is **not green** for the upstream AmberTools 26 auxiliary Python
  distributions: edgembar/ndfes/fetkutils declare NumPy <2, proprep declares
  NumPy <2 and Biopython <1.86, and unrelated pdb2pqr requirements are missing.
  These are recorded upstream metadata conflicts, not hidden with `--no-deps`
  or fixed by downgrading the tested scientific stack without evidence.

### Local evidence and replay commands

Logs, per-step commands, resolved versions/builds and reports live under the
ignored `.artifacts/installation-profiles/` directory. Important bundles:

- `micromamba-standard-licensed-smoke/` — 14/14 including model/Homology;
- `mamba-standard-smoke/`, `conda-standard-smoke/` — scientific 12/12;
- `micromamba-desktop-science-smoke/`, `mamba-desktop-science-smoke/`,
  `conda-desktop-science-smoke/` — scientific 12/12;
- `mamba-standard-modeling/`, `conda-standard-modeling-v2/`,
  `micromamba-desktop-modeling/`, `mamba-desktop-modeling/`,
  `conda-desktop-modeling-v2/` — actual licensed modeling/Homology;
- `py311-science-smoke/` — Python 3.11 scientific 12/12;
- `*-install.log`, `*-configure*.log`, `desktop-resume.log` — manager commands,
  GUI builds and deferred/finished license states;
- `clean-scientific-regressions.log` — **47 passed** in the new Standard env
  (SMILES, GLYCAM roundtrip, CIF/structure input, ACPYPE and Homology contracts);
- `final-fast-tests.log` — **1681 passed, 7 skipped, 31 deselected**;
  focused helper/composition/readiness/doctor/environment checks: **36 passed**.
  The seven skips are one unavailable Biotite and six unavailable GROMACS
  checks; they are not accepted scientific capabilities.
- `helper-automatic-smoke.log` — the helper itself ran all **14/14** licensed
  scientific/modeling steps and finished with profile prerequisites `ready`.
  `helper-configure-automatic-smoke.log` proves the configure-only conda path
  runs model/Homology and finishes with `ready`.
- `distribution-build.log` — wheel/source build passed; a new dependency-free
  wheel environment passed `scripts/smoke_distribution.py`, including catalog
  inheritance. The source archive was checked to contain the TOML catalog and
  exclude the local credential file.

```bash
# Use a NEW prefix, never an existing environment. No activation is required.
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --prefix /path/to/new-standard --non-interactive --skip-modeller

# Finish setup later; KEY_FILE is a user-owned 0600 file. Its VALUE is not argv.
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --prefix /path/to/new-standard --non-interactive --resume --configure-only \
  --modeller-key-fd 3 3< "$KEY_FILE"

micromamba run -p /path/to/new-standard python scripts/smoke_installation_profile.py \
  --with-modeller --output /path/to/new-evidence-directory
```

MODELLER setup was first tested as deferred (`action-required`, exit 2), then
completed using the **operator-provided protected file**, never copied from an
existing environment. The target config is atomically replaced with mode 0600
to break conda package-cache hardlinks. The local key filename is Git-ignored
and explicitly excluded from source distributions. No key value belongs in
commands, logs, manifests or this evidence document.

GUI source verification: **403 tests passed**, typecheck passed. Desktop
builds from fresh `npm ci` use the lockfile and include standalone HTTP smoke.
Source tests used the existing checkout's GUI dependencies; that is distinct
from the new Desktop lockfile-build evidence.

Final Ruff, strict-target mypy, versions, agent/user documentation checks,
Python-3.11 CLI/GUI generators, profile generator and `git diff --check` passed.
CI now covers the stdlib-only composition/helper/readiness contracts and
generated-spec drift. Full CI profile-lane migration remains a later gate.

## WP3 diffusion-runtime increment (after `2ffc9fc`)

The earlier Standard/Desktop implementation was committed as `2ffc9fc` before
starting this increment. New provisioning is in `src/dvbfixer/runtime_install.py`,
reachable from the same `scripts/install_dvbfixer.py --runtime ...` helper.
It clones and patches an owned engine checkout, creates the isolated native/pip
environments through the selected manager, verifies operator assets, runs
actual CUDA/Reference kernels, captures package/build inventories, stages
registration through the existing script, and invokes the existing bounded
runner preflight before atomic launcher publication.

| Linux manager | cc94 sampling + separate Reference refiner | Protenix frozen backend (sampling + refinement) |
|---|---|---|
| micromamba 2.9.0 | clean install, import/kernels, canonical preflight passed | clean install, inference-entry import/kernels, canonical preflight passed |
| mamba 2.10.0 | clean install, import/kernels, canonical preflight passed | clean install, inference-entry import/kernels, canonical preflight passed |
| conda 26.9.4 | clean install, import/kernels, canonical preflight passed | clean install, inference-entry import/kernels, canonical preflight passed |

Both new micromamba backends completed a checkpoint-backed public `model`
execution on the tracked 1UBQ internal-residue deletion, including refinement
and native-filesystem bundle publication. This is an installation/execution
smoke outside those profiles' qualified length/case evidence; it is **not**
general gap acceptance or a transfer of historical cc89 claims.

- cc94 uses the maintained sampling environment and callback patch, cc94
  checkpoint/config digests, and a **new separate** Python 3.13/OpenMM 8.6.1
  Reference refiner. No CUDA refiner or cc89 rotary intervention is introduced.
- Protenix preserves Python 3.13.15, PyTorch 2.13.0/CUDA 12.9 and all existing
  critical `profile-lock.json` identities. The new environment is self-contained,
  unlike the historical backend venv that inherited a pre-existing Standard env.
  Its frozen runner continues to refine in its sampling process; this increment
  does not change that scientific execution layout.
- A native MPS route and separate OpenMM 8.5.2 CPU-refiner input are implemented,
  but installation/execution on macOS arm64 remains **not-tested**. Linux rejects
  that route before creating directories. cc94 MPS scientific acceptance is
  still pending.
- Missing weights/data/Kalign were actually exercised: the helper returned
  `action-required`, retained successful kernel/source/env evidence and did not
  register a launcher. Existing local operator-owned checkpoint/data files were
  then supplied explicitly; no weights or datasets were downloaded or altered.
- Recipe/manager identities are saved for resume. A changed recipe requires a
  new root; source edits/untracked modules are rejected without reset/clean.
  Existing different launchers require explicit replacement; identical private
  launchers resume idempotently. All test registrations used contained runner
  directories under ignored evidence, preserving existing user registrations.

Real sampling caught two issues not exposed by solve/import/preflight alone:

1. Protenix's inference entrypoint additionally needs **Click/Requests**. These
   are explicit in its recipe, and the smoke imports `runner.batch_inference`
   rather than only the network module. The superseded incomplete prefix was
   retained; a new clean prefix was created for the corrected recipe.
2. The current staged refiner received a relative output path from the Protenix
   runner, then tried to relativize it against an absolute workspace. The shared
   refiner now normalizes the output path at entry and rejects paths outside the
   workspace before writing. Relative-output and no-outside-write regressions
   cover this orchestration fix; no scientific algorithm/default was changed.

The first attempted bundle publication on the checkout filesystem also failed
its existing atomic no-replace check. The succeeding public smokes used native
`/tmp/opencode` storage; the rejected attempt and its diagnostics were retained.

Local command/results evidence:

- `cc94-provision-without-checkpoint-v3.log`, `cc94-provision-with-checkpoint.log`;
- `protenix-provision-without-assets.log`, `protenix-provision-v2.log`;
- `cc94-runtime-{mamba,conda}-install.log`;
- `protenix-runtime-{mamba,conda}-install.log`;
- `cc94-runtime-public-model-native-smoke.log`;
- `protenix-runtime-public-model-native-smoke-v3.log`;
- per-step native build/package snapshots and immutable recipe records below
  each newly created runtime's `evidence/` and `runtime-plan.json`.

Final increment verification: **1699 passed, 7 skipped, 31 deselected** in
`runtime-final-fast-tests-v2.log`; Ruff, strict-target mypy, generated profiles,
CLI/GUI schema checks, versions and agent/user documentation checks passed.
`cc94-runtime-idempotent-resume.log` also verifies actual resume with an already
published identical launcher, including its recorded SHA256. Interrupted auxiliary
copies are staged and hash-checked before no-clobber publication, so a partial
copy does not become a resume-visible input.

## WP4–WP6 Linux completion and Apple gate

### Functional results

- Final candidate Standard and Desktop were freshly installed below new
  `/tmp/opencode/dvbfixer-linux-final-*` roots using micromamba, mamba and conda.
  All six finished **ready**, with real licensed model/Homology, 14 scientific
  smokes, and fresh Desktop lockfile build/standalone checks.
- Final installed Desktop tests exposed Node 26's changed
  `fs.cpSync(errorOnExist)` directory behavior. Native evidence copying now merges
  directories explicitly with exclusive file creation, while retaining cache,
  symlink/hardlink and no-overwrite rules. All **405 tests passed** in source
  and in the actual fresh-lockfile installed Desktop for **each** Linux manager.
  Source typecheck and all installed build/standalone checks also passed.
- The first parallel matrix hit a mixed-manager package-cache `.partial` rename
  race and two tool timeouts. Evidence was retained. Sequential owned-prefix
  `--resume` completed the affected routes; no existing environment was removed.
- `linux-slow-acceptance.log`: **31 passed** across the tracked slow input-class
  and ZBS regressions, in the new Standard environment with `--no-solvent` on
  minimize/ZBS. This includes GLYCAM, disulfides, ligand chemistry/SMILES,
  missing atoms, insertion codes, input safety and pathological modeled inputs.
- `linux-psi4-integration-test.log`: the additional real optional two-stage
  separate-runtime RESP test passed with its selected prefix. The main profile
  remains free of PSI4 and PyTorch.
- `linux-extended-acceptance.log`: actual PDBx/mmCIF and small-molecule CIF,
  native Open Babel, GAFF/OpenFF, BCC/PySCF/PSI4 complete ethanol parameterization,
  ANARCI renumbering, DCD/GRO/NetCDF, structural fit, cluster/Plotly, ACPYPE and
  CHARMM paths passed. An attempted Amber19SB `--ff-dir` combination was correctly
  rejected by existing policy; it was corrected to the supported Amber custom-FF
  route rather than relaxing that policy.
- `linux-topology-acceptance.log`: supported custom FF, bundled Amber19SB,
  native `editconf` and strict `grompp -maxwarn 0` passed. The GROMACS
  preprocessing fixture runs zero MD steps; no solvent minimization was added.
- `linux-glycam-acpype.log` and `linux-glycam-grompp.log`: the tracked hinge/CH3
  glycoprotein exported **4117 atoms** through AMBER/GLYCAM, ParmEd and ACPYPE,
  including nonuniform 1–4 scales, then passed strict native GROMACS preprocessing.
- `linux-gromacs-forcefield-tests.log`: **56 passed**, including all six
  previously unexecuted GROMACS water/ion preprocessing cases.
- `linux-tpr-reader.log`: the application input loader read the GROMACS 2026.3
  TPR and matching GRO (**1231 atoms**, one frame). Unsupported/corrupt TPR and
  mismatched atom-count failure paths are covered separately.
- A NumPy ABI RuntimeWarning appeared once in the long slow suite. Native imports
  and the implicated lipid regression passed again with RuntimeWarnings as errors
  (`linux-abi-warning-check.log`); real DCD/GRO/NetCDF roundtrips also passed. It
  was not hidden by suppressing warnings or changing the tested main NumPy stack.
- Both newly provisioned CUDA public-model bundles were consumed through the
  digest-verified publication API: `status == success`, one selected candidate
  each (`linux-runtime-bundle-validation.log`). This remains case-specific
  installation evidence, not new general or MPS scientific qualification.

### Optional components and defects found by real acceptance

- `linux-standard-addons-resume.log`: a new licensed Standard prefix with
  `probe` and `msa-alternatives` finished ready. MUSCLE 5.3 and Clustal Omega
  1.2.4 are actually **bioconda** packages; explicit conda-forge qualifiers were
  corrected after a failed clean solve, and the resumed install passed. Probe
  2.18's actual engine/version path was checked; unselected Probe stays optional.
- MSA now rejects MUSCLE 3 before launching alignment, and the recipe selects
  `>=5,<6`. MAFFT auto preference and explicit Clustal Omega remain unchanged.
- Legacy preparation reports `disabled`, `unavailable`, `failed` or `completed`
  separately for PROPKA/Reduce while retaining documented fallbacks. Standalone
  protonation requires its enabled dependencies before calculations. Tleap's
  explicit backend checks both tools/input chemistry early, honors `--no-propka`
  and rejects `--no-protassign` because Reduce is required for H placement.
- ANARCI Python API, hmmscan and all pressed HMM files are checked together.
  Engine/data failures now abort explicit antibody numbering before output;
  a successful no-domain result remains a legitimate negative.
- `cluster` checks requested plotting before analysis (retaining the documented
  numeric-output fallback), validates atom counts/first frame and gives explicit
  TPR version/matching-GRO/PDB advice. No topology is silently substituted.

### PSI4 is separate; PySCF remains included

The actual optional PSI4 path exposed three defects: psiresp 0.4.2 cannot import
with Pydantic 2, unbounded Libint 2.13.1 segfaulted in the native PSI4 kernel,
and the old worker called a nonexistent `Molecule.from_psi4` API.

The generated separate recipe pins **PSI4 1.9.1 / psiresp 0.4.2 / Pydantic <2 /
Libint 2.9.0**. Its Python 3.9 is isolated and does not change the main
`>=3.11,<3.14` contract. `--runtime psi4` provisions this independent environment
and checks an actual HF native kernel through all three Linux managers; receipts
are `linux-psi4-managed-{install,mamba,conda}.log`.

The worker uses the public QCSchema/AtomicResult bridge and psiresp's own grid and
two-stage fitting API. Input connectivity/atom order is carried across the process
boundary, and ESP is computed on the Job-owned orientation (Pydantic copies nested
models). No new fit implementation or ligand-perception rules were introduced.
Water RESP and complete ethanol GAFF2/RESP export passed; the optional integration
test runs only with an explicitly supplied `DVBFIXER_TEST_PSI4_PREFIX`.

Gaussian input generation also passed (`linux-gaussian-input-smoke.log`);
licensed external Gaussian execution was not requested or claimed. This optional
operator component is not an Apple platform blocker for Standard/Desktop.

### Migration and next host

`environment.yml` is generated Standard; Node is Desktop-only, alternative MSA
is explicit, and `environments/legacy-full.yml` retains the old feature composition.
PySCF/RESP, MDAnalysis/Plotly, ACPYPE/ParmEd and normal chemistry remain in Standard.
No xTB route or ML metapackage is introduced there.

CI now uses explicit Standard/Desktop specs, validates generated aliases, and
executes the isolated PSI4 integration test. MODELLER is configured only from a
protected CI secret when provided; otherwise readiness records action-required
and licensed slow acceptance is not claimed. Environment caching is disabled for
the licensed target; download caching does not receive the private replacement
config. Hosted workflow execution has not been claimed from local checks alone.

Further required platform acceptance now needs a **native Apple Silicon host**.
See [the handoff](plans/installation-profiles-apple-acceptance.md); all native Mac
and MPS rows remain not-tested until executed there. The full plan is not declared
complete at the Linux checkpoint.

Final checkpoint checks are recorded in `linux-handoff-fast-tests.log`:
**1728 passed, 3 skipped, 32 deselected** in Standard with the native GROMACS
test executable available. The three skips are two optional Torch checks and
one independent Biotite comparison, consistent with excluding ML from Standard.
The Torch cases passed separately (**2 passed**, `linux-optional-torch-tests.log`)
and the independent Biotite aggregation agreed exactly in the isolated sampling
runtime (`linux-optional-biotite-check.log`). They are not presented as Standard
dependencies or as empty passing skips.

Source and each installed Desktop manager route passed **405 GUI tests**; final
typecheck/build/standalone checks passed. `linux-final-cluster-plot.log` confirms
actual nonempty HTML output after clustering. The optional PSI4 integration
passed separately with its supplied runtime prefix. Final wheel/source builds
and the dependency-free installed-wheel smoke passed, including stdlib-only
catalog/runtime bootstrap and bundled force fields. The source archive contains
the catalog and excludes the local credential and ignored evidence tree.

Ruff, strict-target mypy, actionlint, version synchronization, agent/user docs,
generated profiles and Python-3.11 CLI/GUI schema checks passed. None of these
local checks is recorded as a native Mac or hosted workflow execution.
