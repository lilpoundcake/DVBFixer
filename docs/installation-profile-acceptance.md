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
- **WP3 implemented (native MPS provisioning/example execution passed; broader scientific acceptance pending):** interactive/noninteractive helper, three manager adapters,
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
- **WP5 Linux acceptance executed; native Apple increment underway:** final main-profile and CUDA/QM manager
  routes, slow structural regressions, extended format/ligand/topology/analysis
  checks and native GROMACS preprocessing passed. Native micromamba/mamba
  Standard/Desktop, MSA alternatives, separate PSI4 and cc94 MPS installation/
  example execution have passed checks below. Native conda/Probe restrictions
  and broader cc94 MPS scientific qualification remain open.
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

## Native Apple acceptance: bootstrap (2026-10-29)

Work started in a separate `feature/installation-profiles` worktree at
`6662503c69a1b101ce5c8198e7594ca1a86ecf43`, with a clean initial Git status.
Existing environments and the original checkout were preserved.

- Host: native macOS 26.6.2 (25G83), arm64; bootstrap Python 3.12.7.
  Explicit `PYTHONPATH="$PWD/src"` resolves DVBFixer 0.10.1 to this worktree.
- Installed managers: micromamba 2.0.2, mamba 1.5.9, conda 24.9.2.
  Standard helper `--dry-run` passed with micromamba and conda. Mamba was
  rejected before prefix creation because the helper requires mamba >=2.0;
  the existing manager was not upgraded or replaced.
- Native bootstrap contract tests: **34 passed** using
  `PYTHONPATH="$PWD/src" python3 -m pytest tests/test_install_helper.py
  tests/test_installation_profiles.py tests/test_runtime_install.py -q`.
  Generated installation profiles passed their `--check` with the same
  explicit source import path. These are not installed scientific smokes.
- Only approximately **11 GiB** was available on the target volume at initial
  inspection. Installation sizing must precede creation of new scientific,
  Desktop and MPS prefixes. No existing environments/caches were removed.
- The native micromamba Standard `create --dry-run --json --override-channels`
  solve completed successfully: 244 LINK packages and 179 FETCH packages
  (392,046,727 declared download bytes at solve time). Selected native builds
  include Python 3.13.18 `hb59dee6_100_cp313`, OpenMM 8.6.1
  `py313h6e8b90c_0_khronos`, MODELLER 10.8 `py313h8a0c2ca_0`, AmberTools 26.0
  `cuda_None_nompi_py313h47fe9ef_100`, plus noarch OpenFF Toolkit Base 0.19.0.
  This is resolver evidence only, not installed/import/scientific acceptance.
  Declared package sizes do not establish the unpacked environment footprint;
  FETCH requirements can change after cache cleanup.
- Local bootstrap evidence is retained under
  `.artifacts/installation-profiles/apple/` (`host.json`,
  `bootstrap-tests.log`, `helper-dry-runs.json`, manager dry-run logs).

Standard/Desktop scientific and licensed workflows, installed Desktop/browser,
add-ons, optional PSI4 and cc94 MPS remain **not-tested** on this host. Neither
bootstrap tests nor helper dry runs close those platform gates. Licensed
acceptance requires an operator-provided protected MODELLER credential source;
no key is retrieved from an existing installation.

Subsequent storage cleanup was separately authorized by the user: removal of
the large multi-checkpoint download archive while preserving verified extracted
cc94/cc89 assets, lossless verified recovery archives for historical gfnff files
and handoff essentials, and `micromamba clean --tarballs --yes`. No existing
environments or unpacked package-cache entries were removed. After tarball
cleanup, measured free space was **23.85 GiB**. Detailed local cleanup evidence
is in `artifact-audit/README.md` and `micromamba-clean-tarballs.log` below the
Apple evidence directory.

Native installation has now been started with micromamba in a new owned prefix,
initially using `--skip-modeller`. The operator subsequently supplied a protected
credential file; ownership and mode 0600 were checked without reading its value
into reports. Licensed configure-only/model/Homology acceptance is still pending.
The broader bootstrap readiness test selection returned **54 passed, 2 failed**:
the existing bootstrap Python lacks Gemmi and SciPy. These are not failures of
the new Standard prefix (installation is pending); repeat the selection there.
The separate cc94 MPS helper dry run accepted the retained operator checkpoint
and config identities; it did not provision or execute MPS. Evidence:
`bootstrap-readiness-tests.log`, `micromamba-standard-install.log`, and
`cc94-mps-plan.json` below the Apple evidence directory.

### Native micromamba stalled attempt and channel isolation fix

The initial helper installation made no prefix and emitted no native progress
for over 85 minutes. A native `sample` captured an approximately 43.7 GiB
physical footprint; symbols were unavailable, so the stack does not identify
the underlying resolver operation. The operator's rc files add a channel not
declared by the profile. Unlike the successful dry-run command, the helper's
micromamba create/update commands lacked `--override-channels` and could merge
those user channels. No credential-bearing channel URLs are recorded here.

Stopped only the owned stalled micromamba child with SIGTERM, preserving its
log and sample. The micromamba adapter now passes `--override-channels` for
both creation and resume/update, without modifying user configuration. Manager
argv and runtime/helper/composition tests pass (**34 passed**). This correction
does not by itself prove the cause of the high memory use or native readiness.
An owned-prefix resume was started with protected MODELLER fd input and a
30-minute process-group timeout; evidence is in
`micromamba-standard-install-v2.log`. Its result remains pending.

That resume reached the 30-minute limit and was terminated (exit 124), still
without creating the target prefix. At a 16-minute sample it used one CPU core
and approximately 1.5 GiB RSS. Channel isolation therefore did not resolve the
stall; it must not be reported as a proven root-cause fix. No native scientific
or licensed smoke had run at that point.

Bootstrapped native micromamba **2.9.0** in a new separate tools directory,
without replacing the existing 2.0.2 binary or installing into an existing
environment. The conda-forge osx-arm64 distribution was checked against the
release API's SHA256
`500f5074feb8d02c4296ef9921c3650ed2874171805a9fbb8fbb53896433646b`.
Its executable version check passed. A strict-priority, profile-channel-only
Standard dry run was started with a five-minute process-group limit; its logs
and elapsed-time record are `micromamba-2.9-standard-solve*`. This remains
diagnostic work, not platform acceptance. Focused helper/catalog/runtime and
readiness checks now pass **41 tests** in bootstrap Python.

Micromamba 2.9.0 completed that strict native dry run successfully in **30.76
seconds** (245 LINK packages, 182 FETCH packages). Selected native builds retain
Python 3.13.18, NumPy 2.4.6, OpenMM 8.6.1 and AmberTools 26.0; SciPy resolves to
1.18.1 `py313hccac1f9_2`. This is evidence of materially different behavior with
the newer binary on this host, not a diagnosis of the old binary's internal
failure. A real helper resume was then started with the separate 2.9.0 tools
directory prepended to PATH, protected MODELLER fd input and the same 30-minute
process-group limit. Installation and science remain pending; logs/timing are
`micromamba-standard-install-v3*`.

The 2.9.0 real installation completed its native package transaction and reached
the pip stage (including native arm64 Open Babel/ACPYPE/PySCF distributions).
Scientific and licensed readiness are still pending that attempt's completion.
At the user's explicit request, `~/bin/micromamba` was also atomically updated
from 2.0.2 to the verified 2.9.0 binary, with the old executable retained in the
separate tools directory for rollback. No existing environments or user rc files
were modified by the binary update; `user-micromamba-update.json` records it.

Native source GUI verification passed a fresh `npm ci`, typecheck, **405 tests**
in 57 files, client/server build and standalone HTTP smoke. Evidence:
`source-gui-{ci,typecheck,tests,build}.log` and `source-gui-results.json`. This is
source-worktree evidence, not installed Desktop or real-browser acceptance.

### Native Standard first successful installation

The micromamba **2.9.0** helper resume finished successfully (exit 0) in
**445.45 seconds**, including protected MODELLER configuration and actual
licensed execution. All **14/14** scientific smoke steps passed: chemistry/
OpenFF/GAFF2, real PySCF HF/6-31G* RESP, ANARCI/HMMER/data, DCD/GRO/NetCDF/
structural fitting/Plotly, both preparation backends, native obminimize,
protonation, diagnose, Reference minimize `--no-solvent`, topology, MAFFT,
MODELLER gap modeling and Homology. Repaired/minimized outputs also passed
the smoke's nonempty/all-L assertions. Native linear algebra used OpenBLAS;
`doctor --profile standard` reported `ready`. Probe and alternative MSA are
not selected by this baseline; their absence is not counted as a passed add-on.

Evidence: `micromamba-standard-install-v3.log`, its timing JSON and copied
`standard-licensed-report.json`. The immutable prefix-owned original report is
under `share/dvbfixer/installation-smoke/87ee687ee9734faaa8d1c7e677af6ddf/`.
Installed native/Python versions include Python 3.13.18, OpenMM 8.6.1,
NumPy 2.4.6, SciPy 1.18.1, RDKit 2026.3.1, MDAnalysis 2.11.0, PySCF 2.14.0,
ACPYPE 2026.10.15, Open Babel wheel 3.1.1.23 and MODELLER 10.8.

| Native Apple route | Current executed result |
|---|---|
| micromamba 2.9.0 Standard | **passed**, licensed 14/14 smoke; doctor ready |
| micromamba 2.9.0 Desktop | **passed**, licensed 14/14 smoke, fresh GUI build/HTTP smoke, 405 installed tests |
| Source GUI | **passed**, fresh npm ci/typecheck/405 tests/build/HTTP smoke |
| Installed Desktop/browser | **passed**, real Chrome workspace import and scientific diagnose job |
| mamba 2.10.0 Standard/Desktop | **passed**, licensed 14/14 each and 405 installed Desktop tests |
| conda 26.9.4 profiles | **failed** Standard attempts; Desktop not-tested pending resolution |
| MSA alternatives | **passed**, separate licensed prefix and 30/30 extended checks including explicit MUSCLE/Clustal Omega |
| Probe add-on | **failed** native solve: no bioconda Probe 2.18 osx-arm64 package |
| Optional PSI4 | **passed**, separate native kernel, integration RESP and full ethanol export |
| cc94 MPS | **passed** corrected provisioning/preflight and one public installation case; broader research qualification pending |

Extended Standard acceptance was started separately with native GROMACS
**2026.0**, whose version reports ARM_NEON_ASIMD and CPU-only execution.
This existing operator installation is external to Standard; its presence does
not make GROMACS a profile dependency. The extended run passed **27/27** steps
(`standard-extended/report.json`), including real BCC/PySCF full ethanol
parameterization, native CIF conversions, antibody numbering, trajectory formats,
cluster/Plotly, ACPYPE, custom Amber/Amber19SB/CHARMM and strict GROMACS
preprocessing. The application also read the generated 2026.0 TPR with its
matching GRO (**1,231 atoms, one frame**).

Native Standard fast tests: **1,728 passed, 3 skipped, 32 deselected**
(`standard-fast-tests.log`). Skips are two optional Torch tests and one Biotite
comparison; these remain unexecuted here, not accepted Standard capabilities.
The run recorded two upstream SWIG deprecation warnings and four unraisable
`batch._Mirror.flush` exceptions when a captured target stream was already
closed. Tests passed, but those warnings are retained, not suppressed or called
warning-free acceptance. Ruff, strict-target mypy and **41** focused source
helper/readiness tests passed. Version, agent-doc, profile and GUI generators
passed. CLI-reference checking under Python 3.13 reports argparse alias-format
drift; its detailed diff is retained without rewriting generated docs. The
bootstrap Python 3.12 cannot complete that check because it lacks OpenMM.
A canonical Python-3.11 generator check remains outstanding.

Native micromamba Desktop installed in **135.04 seconds**, with real licensed
science, fresh lockfile install/typecheck/client/server build/HTTP smoke and
**405 installed GUI tests** in 57 files (`micromamba-desktop-install*`,
`installed-desktop-tests.log`). Source npm ci reported zero vulnerabilities;
the build retained Vite's large-chunk warning.

Installed Desktop real-browser acceptance used existing native Google Chrome
**154.0.8037.98** through Playwright Core in headless mode, with its own fresh
browser context and workspace-data directory. The installed Node server bound
only to loopback and invoked the real installed Desktop `dvbfixer`, not the
version-only CLI stub used by the build's HTTP smoke. Through the UI it created
a workspace, imported tracked 1UBQ, rendered the structure, selected Analysis /
Diagnose and ran a real managed scientific job. The persisted job finished
`succeeded`, exit 0, and the UI displayed its report. No page errors occurred.
Evidence: `browser-smoke-v3.log` and `browser-1793286991875/` (screenshot, server
log, report and owned workspace/job files). Two earlier harness attempts failed
on ambiguous/wrong selectors and are retained separately; they are not counted
as application failures or passing browser runs.

Separate native tools with mamba **2.10.0** and conda **26.9.4** were installed
in a new tools prefix using a distinct package cache. The original mamba/conda
executables and environments were not upgraded. Their four Standard/Desktop
helper installations and installed Desktop tests were queued sequentially;
results remain pending. Verified mamba 2.10 help supports channel overrides,
so the same create/update channel-isolation fix now covers both libmamba
adapters; conda's argv is unchanged. cc94 MPS provisioning was also started in
its own new runtime root using the operator's verified retained checkpoint and
a separate cache, with fallback disabled and a separate CPU refiner. Runtime
preflight and actual MPS sampling remain unaccepted until their execution.

### Additional native results and defects

- Mamba **2.10.0** Standard/Desktop both finished ready, with actual licensed
  14/14 science/model/Homology in each profile (120.56 / 141.52 seconds) and
  **405 installed Desktop tests**. Evidence: `additional-manager-matrix.json`,
  `mamba-{standard,desktop}-install.log`, `mamba-desktop-installed-tests.log`.
- Conda **26.9.4** first failed at pip because inherited caller activation state
  left tools-Python ahead of the target prefix after `conda run`. Reproduced
  the wrong interpreter, then cleared only child `CONDA_PREFIX`, nested-prefix,
  `CONDA_SHLVL`, default-env and prompt state: the same real conda invocation
  then selected the intended target Python. `child_environment` now implements
  this without modifying the user's shell; regression tests preserve the
  selected package cache and parent activation variables.
- The initial conda attempt also reported an existing corrupted blosc cache
  entry and upstream package path overlaps. A second attempt used a new prefix,
  a fresh dedicated cache, enabled safety checks and path-conflict prevention.
  It correctly stopped on gettext/libintl dylib and tk/xorg-libx11 header
  collisions. SHA256 comparison showed that the overlapping dylib and two
  headers are NOT identical (one remaining header is identical). No `--clobber`
  workaround, cache deletion or claim of conda readiness was made. Both failed
  prefixes/logs are preserved (`conda-standard-install.log`,
  `conda-standard-v2-install.log`, `conda-v2-matrix.json`). Desktop via conda
  remains not-tested. This upstream transaction issue is separate from the
  corrected child-activation bug.
- Separate native PSI4 provisioning passed in **117.84 seconds**: actual
  PSI4 1.9.1 / psiresp 0.4.2 / Pydantic <2 / Libint 2.9.0 installation and HF
  kernel. The explicitly selected-runtime integration test passed (**1 passed**,
  `psi4-integration.log`), and complete ethanol GAFF2/two-stage RESP export also
  passed (`psi4-ethanol/command.log`). The Python 3.9 runtime contains no DVBFixer
  install; main Python constraints and included PySCF are unchanged.
- The combined `msa-alternatives` + `probe` clean install failed before creating
  an environment: `probe=2.18` is unavailable from the declared channels on
  osx-arm64. No Rosetta/x86 or different Probe version was substituted. A new
  MSA-only Standard prefix was started separately; failed-prefix ownership
  markers are not repurposed with different add-ons. Evidence:
  `standard-addons-install.log` and `standard-msa-*`.
- The first cc94 MPS runtime installed and passed real PyTorch 2.6.0 MPS tensor
  math plus the separate OpenMM CPU-refiner kernel, but full runner preflight
  failed: sampling imported validation, whose quality reference arrays require
  MDAnalysis distribution data. The Apple sampling recipe lacked that dependency
  even though CUDA already declares it. Added the same MDAnalysis 2.10.0 and
  psutil 7.2.2 requirements; no scientific gate/default was relaxed. Sampling
  kernel checks now import the shared validator so this class of omission is
  caught before staging. A second attempt exposed incorrectly indented/folded
  pip entries in the interim recipe; that was corrected and tests now assert
  both dependencies are separate correctly indented pip list items. Both failed
  runtimes are retained; the corrected third clean root is provisioning. No
  launcher was published from either failed runtime and no MPS sampling or
  candidate acceptance is claimed. Direct preflight stderr and all kernel logs
  are retained. Focused source helper/policy tests passed **47 tests** before
  the final indentation regression was tightened. A correctly parsed final
  YAML was also checked with PyYAML, confirming separate dependency entries.
- Native slow scientific regressions completed successfully in **538.34 seconds**;
  the exact test summary is retained in `standard-slow-tests.log`. The new
  environment used source-worktree imports, two CPU threads and the tracked
  tests' `--no-solvent` development policy. Exact summary: **31 passed, 1
  skipped, 1,731 deselected**. The skipped optional PSI4 integration subsequently
  passed separately with its explicitly supplied new native runtime.

### Native MSA, generator baseline and MPS execution checkpoint

- A separate licensed Standard prefix with **only `msa-alternatives`** finished
  ready, then passed **30/30 extended steps**, including MAFFT auto, actual
  explicit MUSCLE 5/Clustal Omega alignments, optional PSI4 full parameterization
  and GROMACS. The unprovisionable Probe add-on was not silently removed from a
  requested composition: this is a distinct new prefix and selection.
  Evidence: `standard-msa-results.json`, `standard-msa-extended/report.json`.
- A new managed **Python 3.11** generator-only environment completed both
  `gen_cli_reference.py --check` and `gen_gui_spec.py --check`. This confirms
  the Python-3.13 drift was argparse rendering, not a changed CLI surface.
  Generated reference files were not rewritten. Evidence:
  `generator-py311-results.json`, `gen_cli_reference.py-py311.log`,
  `gen_gui_spec.py-py311.log`. This environment is not represented as a full
  scientific Standard profile acceptance route.
- The corrected third cc94 MPS runtime provisioned in **132.57 seconds**,
  passed native sampling/CPU-refinement kernels and the canonical runner
  handshake, then published its contained launcher. Earlier failed roots were
  retained, not rewritten after changing their immutable recipe identity.
- From the new Standard CLI, the public 1UBQ internal-residue deletion ran in
  **38.11 seconds**. `read_published_bundle(...).selected()` succeeded only
  after its independent digest/geometry validation: bundle `status=success`,
  one seed-7 candidate, **500** denoising updates on actual **MPS**, fallback
  disabled, correct cc94 checkpoint/source identities and source-physical/native
  relative features. Fixed-heavy-atom RMSD/max displacement were both **0.0 Å**
  and hard-gate failures were empty. The candidate passed at **raw** stage;
  refinement mode was `none`, so this case does NOT establish an actual
  scientific minimization through the separate CPU refiner (whose native
  kernel and preflight did pass). Profile evidence labels remain experimental/
  MPS-acceptance-pending. Evidence: `cc94-public-model*`,
  `cc94-public-bundle-validation.json` and the complete `cc94-public-bundle/`.
- An attempted case-bounded five-seed replay (7, 19, 42, 101, 202) used the
  digest-verified archived single-gap 9CQD request, changing only profile,
  native chain phase and seed while checking 58 generated / 3,441 fixed atoms,
  1,399 conditioning atoms and 189 tokens. **All five requests were rejected
  before sampling** by current public admission:
  `sequence-placement-gap-coverage-mismatch`. Thus **zero scientific seeds
  executed**; these are unsupported admissions, not a 0/5 MPS geometry result
  and certainly not passing scientific evidence. The old CUDA research harness
  uses a research-only selected-gap path and is not a portable public admission
  proof. Do not bypass scope/coverage gates or retry altered favorable inputs to
  close the MPS gate. Requests/results remain in `cc94-9cqd-mps-n0/` and its log.
  A separately specified valid representative model/Homology/multi-context
  evaluation remains necessary before broader cc94 MPS scientific claims.

Final source-worktree fast tests after the adapter/dependency regressions:
**1,730 passed, 3 skipped, 32 deselected** (`final-fast-tests.log`). The same
six SWIG/closed-capture warnings were retained. Final Ruff and strict-target
mypy passed; Ruff additionally reports the pre-existing removed UP038 ignore.
Wheel and source distributions built successfully. The source archive contains
the TOML catalog and excludes `.artifacts/` and the credential filename; a new
dependency-free installed-wheel environment passed `smoke_distribution.py`
(`distribution-build.log`, `wheel-smoke.log`). No existing environments or
artifacts were deleted to obtain these results. The native conda/Probe and
broader MPS gates above remain unresolved; the whole plan is not complete.
