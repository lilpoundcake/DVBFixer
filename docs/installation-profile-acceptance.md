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
| Probe | managed add-on solve **failed**; separate pinned native source build/doctor/electron+nuclear diagnostics **passed** |
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

### Follow-up after native checkpoint commit `71f99e9`

The checkpoint was committed and pushed to `origin/feature/installation-profiles`
at the operator's request. Subsequent native Probe work did not alter the
existing source checkout, scientific prefixes, user PATH or manager configuration.

**Native Probe manual route passed.** The previously retained upstream source
has a clean tracked tree at `091d117f5bceef87389002dccfe26df78516fdda` from
`https://github.com/rlabduke/probe.git`. Cloned it with `--no-hardlinks` into a
new native tools checkout and built with Apple's compiler via
`make probe CFLAGS=-O2`. Apple clang 21.0.0 reports an arm64 target; `file`
confirms the result is a Mach-O 64-bit arm64 executable. The exact supported
`-version` result is **`probe.2.18.211005`**; binary SHA256 is
`d4a195af8604276bd8bad6d0acbd38c12161e34b4e8995fe6f713d250885e007`.
Compiler-version changes may change the binary digest; this is executed
provenance, not a cross-toolchain bit-reproducibility claim.

Placed the new binary directory on PATH only for child invocations from the
new MSA Standard prefix. Doctor with explicit `msa-alternatives` + `probe`
reported ready. Both electron and nuclear radius diagnostics used
`actual_engine=molprobity-probe`, exact supported version, `probe_status=complete`
and `fallback=false`. Each diagnostic returned 1 because actual structure
findings were present, an accepted diagnose result—not a failed Probe process
or a clean-geometry claim. The no-Probe control used
`actual_engine=dvbfixer-python-steric`, `probe_status=unavailable`, `fallback=true`
and the explicit executable-not-found reason. The engines/radii produced their
own actual findings; they were not asserted numerically equivalent.

Evidence: `native-probe-build.log`, `native-probe-build.json`,
`native-probe-doctor.log`, `native-probe-results.json`, and
`native-probe-v2-{electron,nuclear,fallback}.{json,log}`. The first local harness
assertion expected the wrong status spelling (`completed` vs `complete`);
its evidence was retained and a new run used the actual report schema.
This establishes a documented **manual source-built** native route, not a
successful conda-based `--addon probe` solve; that managed packaging gate remains.

**Conda investigation continued without bypassing conflicts.** Reverse native
dependency inspection identified the gettext collision: conda's gawk 5.3.0
requires monolithic gettext >=0.21.1,<1, while libglib 2.90.1 requires libintl >=1.
The successful micromamba solution instead uses gawk 5.4.1 and the split
libintl 0.25.1 family. This explains one graph-dependent collision but does not
fix the independent Tk/X11 ownership conflict. Reviewed all three available
8.6.13 `_3` arm64 Tk builds, verifying downloaded archive SHA256 values against
the upstream API and reading their `info/paths.json` without installing them.
All still contain the same colliding Xlib.h/Xutil.h/cursorfont.h paths; blindly
pinning an older Tk build is therefore not a supported fix. Retained package
metadata and overlap evidence in `*-native-repodata.json` and
`tk-older-build-overlap-review.json`. No compatibility pins, custom package
repacking, header deletion, `--clobber` or safety-check weakening was introduced.
Both solver outputs contain Tk 8.6.13 and xorg-libx11 1.8.14; conda's link-time
ownership checks expose shared headers not rejected by the successful libmamba
transactions. Those libmamba results establish the executed scientific routes,
not absence of every third-party file overlap or identical transaction-safety
semantics between managers. See `conda-native-collision-graph.json`.

**Managed add-on preflight now fails early and explicitly.** The catalog records
Probe's managed platform restriction and its manual guide link. The installer
checks that restriction before manager discovery or prefix reservation, including
dry runs; it does not silently omit a selected component. The generated
dependency matrix now distinguishes managed provisioning from readiness checks
for an operator-provided executable. Linux Probe and Apple MSA selections retain
their existing behavior, while manual native Probe still passes source-worktree
doctor with both explicit add-ons. A real native rejected dry run created neither
prefix nor ownership marker (`native-probe-guard-check.log`).

Follow-up focused tests: **51 passed** (`probe-followup-tests.log`), covering
the early no-side-effect rejection, unchanged supported selections, catalog/doc
platform rendering and retained operator-tool readiness. Ruff (after correcting
test import formatting), strict-target mypy, generated profile/matrix drift,
user/agent doc checks and the real source-worktree doctor passed. These follow-up
changes are separate from the already pushed `71f99e9` checkpoint.

**Public multi-seed MPS installation stress case passed 5/5.** A separate
predeclared 1UBQ case removed ATOM records A10–A13 and all HETATM records from
the tracked source, retaining the full target FASTA. Recorded source/input
digests, seeds **7, 19, 42, 101, 202**, and one attempt per seed before running
the public `model` CLI. This is not an altered retry of the unsupported partial
9CQD requests and is not presented as that study's N0/N1 comparison.

The run completed in **59.17 seconds** and publication verification accepted
all five candidates. Each records **500 MPS updates**, fallback disabled,
empty independent hard-gate failures and **0.0 Å** maximum fixed-heavy-atom
displacement. Unlike the first one-residue installation example, every candidate
required and successfully completed **generated-only refinement** in the
separate pinned OpenMM 8.5.2 CPU worker. All raw and stage-1 artifacts and
refinement reports remain in the verified bundle. The normal publication ranking
selected seed 42, but all five results were retained and assessed; no scientific
attempt was retried or discarded. Evidence:
`cc94-ubq-gap4/{protocol,timing,verified-results}.json`, `model.log` and `bundle/`.
This extends native public model/refiner installation evidence, not general
terminal, interchain-context or cc94 MPS scientific qualification.

**Closed-capture warning fix verified.** Native fast acceptance
exposed `batch._Mirror` being finalized after its borrowed pytest capture target
was closed. Its flush now skips a closed target and tolerates a target closing
during flush, while still propagating ValueError/OSError from an open target.
It never closes the borrowed target; captured diagnostics remain available.
Five small lifetime/error regressions were added. Targeted batch/installer/
composition tests passed **42 tests** (`batch-probe-followup-tests.log`), and
the five isolated stream tests also passed with unraisable warnings promoted to
errors in bootstrap Python. After installing the current source into the new
Standard test prefix, the full source fast suite passed **1,739 tests, 3 skips,
32 deselections** (`followup-fast-tests.log`). The four closed-capture warnings
are gone; only the two upstream SWIG deprecations remain in the test summary.
Ruff and strict-target mypy passed. The new pure stream-lifetime tests were
added to the lightweight CI selection without importing scientific backends;
hosted CI was not executed and local actionlint is unavailable. FD-level runtime
teeing and diagnostic summaries are unchanged.

**Public native Homology/MPS path passed.** A separate template plan masked
the same A10–A13 interval using zero-based half-open masks `[0,9)` and `[13,76)`
on the tracked 1UBQ template. Ran the actual `homology --template-plan --backend
diffusion` CLI with seed 7, without a custom scientific materializer or adapter.
It completed in **18.98 seconds** and wrote the normal
`output_homology_diffusion` bundle. Publication API verification confirmed
`status=success`, 500 MPS updates, fallback disabled, actual geometry-gated
refinement, empty hard-gate failures and fixed-heavy-atom RMSD/max displacement
of **0.0 Å**. Evidence: `cc94-homology-gap4/{template-plan,timing,verified-results}.json`,
`homology.log` and the complete bundle. This validates the native public
template-plan/materialization/runner/refiner/publication integration for one
controlled case; multi-chain/context research qualification remains pending.

**Conda safety policy is now enforced in child commands.** The helper sets
`CONDA_SAFETY_CHECKS=enabled` and `CONDA_PATH_CONFLICT=prevent` only in its child
environment, regardless of permissive caller values; the parent's configuration
and environment are untouched. A real conda 26.9.4 attempt deliberately inherited
`disabled`/`clobber` from the test launcher, yet stopped at transaction
verification on the known upstream overlaps, before executing the transaction.
This is a successful safety-guard test, not a passing conda profile installation
(`conda-safety-guard.log`, `conda-safety-guard-results.json`).

A new micromamba Standard prefix under the current helper completed in **128.29
seconds**, with **14/14 licensed scientific/model/Homology checks passed**, from
a separate fresh-native package cache despite the same deliberately permissive
caller values. Evidence: `safe-standard-install.log`,
`safe-standard-install-results.json`, `safe-standard-scientific-report.json`.
No manual cache cleanup, old-prefix removal, package repacking or file-conflict
repair was used to obtain this result; native manager transactions were retained
as executed, with the manager-safety distinction noted above.

Actual native `pip check` is **not green** (`native-standard-pip-check.log`):
AmberTools auxiliary distributions proprep/ndfes/fetkutils/edgembar declare NumPy
<2, proprep also declares Biopython <1.86, and proprep/packmol-memgen require the
uninstalled pdb2pqr package. These parallel the documented Linux auxiliary
metadata issues and are not hidden with `--no-deps`, ignored as a passing check
or repaired by arbitrarily downgrading the scientifically tested main stack.
The accepted 14-step route does not establish every auxiliary AmberTools command.

Final safety-policy source fast suite: **1,740 passed, 3 skipped, 32 deselected**,
with two upstream SWIG warnings and no unraisable capture warnings
(`safety-followup-fast-tests-v2.log`). A preceding harness invocation omitted the
prefix bin directory from PATH and failed six tests that spawn bare `dvbfixer`;
its log is retained as `safety-followup-fast-tests.log`, not counted as a code
regression or a passing run. The corrected invocation explicitly prepends the
new Standard prefix. Strict mypy, Ruff, generated profile/matrix and Python-3.11
CLI/GUI checks, versions, user/agent docs, source/wheel builds and a new
dependency-free installed-wheel smoke all passed. The installed-wheel smoke
also asserts the packaged catalog restriction and conda child safety defaults.
CI YAML parsed locally; actionlint and hosted workflow execution remain
not-tested on this host. An interpreter-exit SWIG swigvarlink deprecation is
retained in the logs in addition to the two test-summary SWIG warnings.

### 1TIM source geometry and bounded input repair

The tracked RCSB 1TIM dimer was used for a controlled internal interface deletion:
only A73–A78 (target indices `[71,77)`, 44 heavy atoms) were withheld in a separate
copy; both chains and their full 247-residue target sequences were retained. Window
selection used native interface contacts before sampling. The fixture was not
modified. This is case-specific execution evidence, not general reconstruction
accuracy or a replacement for the archived 9CQD research requests.

The initial unprepared `local` cohort published seeds 7, 19, 42 and 202, but not
101. The then-current pipeline removed its private workspace on partial success;
seed 101's precise outcome cannot be reconstructed. This is an evidence loss,
not a justified geometric-failure classification. The planned `none` control
was rejected before sampling because that profile requires its source-selected
motif; it is not a 0/5 scientific arm. All remaining initial files are retained
under `cc94-1tim-interface/`. Partial success now retains private diagnostics;
prepared-input runs retain them even when every candidate passes. Homology also
preserves its mosaic workdir on failure or when input-repair evidence is present.

New pre-sampling validation found **29 severe clashes and 150 warnings** in the
measured source itself. The user authorized automatic local minimization before
diffusion. `bounded-observed-input-repair-v1` selects only observed residues
named by exact affected-atom evidence, with a 64-residue / 25% locality ceiling.
It uses restrained, generated-mask-independent CPU minimization without solvent
or stochastic restarts. Residues outside the selected mask remain byte-identical.
Per-residue movement bounds are 0.5 Å backbone RMSD, 0.5 Å Cα displacement and
1.25 Å maximum heavy-atom displacement. Temporary force-field fragments have
standard termini at known sequence gaps rather than a fictitious C–N spring
across missing residues; those temporary atoms/termini never enter the source.
Unsupported or excessive damage fails with the source and attempted repair
retained. Frozen requests are never silently repaired or reinterpreted.

The first native CPU attempt repaired **47 observed residues** in **11.08 s**:
**zero severe findings** after independent whole-source revalidation, with
maximum per-residue backbone RMSD **0.429 Å**, Cα displacement **0.244 Å**, and
heavy-atom displacement **1.176 Å**. The **109 remaining warnings are retained
and shown**, not counted as zero findings. The independent audit recomputes
selection, movement and all outside-mask records from the original bytes.
Evidence: `cc94-1tim-auto-repair-v1/` and the dedicated real 1TIM regression test.

A separate public `model` run started from the original damaged deletion input
and invoked preparation automatically, before motif selection. All five
predeclared seeds **7, 19, 42, 101, 202** passed, one attempt each, in **182.64 s**
total. Each records **500 actual MPS updates**, fallback disabled, chain-B
conditioning and successful generated-only CPU refinement. Independent candidate
hard-gate failures are empty; maximum fixed-heavy-atom displacement is **0.0 Å
relative to the explicitly repaired source**, not the raw original. The normal
publication ranking selected seed 42; no failed/unfavorable attempt was rerun
or replaced. The input original, prepared source, movement report and final
request snapshot are digest-bound into every candidate's publication evidence;
the reader independently re-audits them, rather than trusting a repair-success
flag. Original and new cohorts, preparation directories and all private runtime
diagnostics remain retained.

Evidence: `cc94-1tim-auto-repair-public-v1/{protocol,execution,verified-results}.json`,
`model.log`, `bundle/` and `work/`. The protocol records input and implementation
SHA256 values before execution. This establishes the executed native public
preparation/motif/sampling/refinement/publication path for this case, not broad
cc94 MPS scientific qualification or a causal comparison with the unsupported
`none` arm. Upstream conda/managed-Probe/AmberTools metadata limitations remain.

The exact executed implementation and inputs were copied into
`implementation-snapshot/` and checked against the pre-run protocol hashes
(`snapshot-integrity.json`). A public frozen `--diffusion-resume` invocation
also passed (**11.16 s**) using the digest-verified final request copied from
the bundle: no new preparation/runtime directories, no coordinate changes and
no resampling (`frozen-resume.log`, `frozen-resume-check.json`). New prepared-input runs also
retain a directly replayable `request-final.json` beside the prepared source
and repair evidence; this is the post-context request, not the intermediate
pre-context snapshot.

Verification of the current source: **1,761 fast tests passed, 3 skipped,
33 deselected** (`input-repair-fast-verified.log`), plus the executed slow real-1TIM
input-repair regression and the existing boundary-refinement regression.
GUI: **407 tests passed**, typecheck passed. Ruff, strict-target mypy,
CLI/GUI generation checks, version and user/agent documentation checks passed.
The final source/wheel builds and a fresh dependency-free installed-wheel smoke
passed; local evidence and credentials are excluded from the source archive.
The 48-test diffusion CLI subset also passed after adding an explicit assertion
that frozen requests never call automatic preparation. No broad scientific
qualification or hosted-CI pass is inferred from these checks.

## Native installed-package finalization (2026-10-29)

Evidence root: `.artifacts/installation-profiles/apple/native-finalization/`
on the native Mac checkout. These local runtime assets are intentionally not
Git payloads; keep the Mac checkout, environments and evidence when continuing
on Linux. Previous package files are preserved in verified tarball backups;
the entire previous Desktop GUI/data/dependencies tree remains at
`share/dvbfixer/gui-before-native-finalization`. No environment or input fixture
was deleted. Package refreshes used built wheels with `--no-deps`, no source
overlay or dependency changes. A second Python-only refresh installed the
heavy-atom coincidence guard; the GUI did not need another replacement.

### Installed science and Desktop

- Standard: **28/28** licensed extended checks, including real PySCF/RESP,
  separate PSI4/RESP export, GAFF/ACPYPE, both prepare backends, native GROMACS,
  formats/analysis, MODELLER and Homology; Desktop: **14/14** licensed checks.
  Both profile doctors returned 0. `pip check` still returns 1 for the recorded
  upstream AmberTools auxiliary metadata conflicts; passing executed science
  does not make that package-metadata gate green.
- Installed Desktop: **407 GUI tests**, typecheck, client/server build and
  standalone smoke passed. Actual Chrome **154.0.8037.98** passed workspace
  import/diagnose, public Protenix context choices `auto/local/full`, and hidden
  Protpardelle context overrides. Browser attempts v1/v2 were harness selector
  failures and remain retained, not scientific failures.
- Browser v3 exposed a diagnostic gap: collapsed bonded CA/CB reached the input
  optimizer, which refused with NaN before sampling. The shared geometry layer
  now detects heavy-heavy coincidences (including bonded/cross-chain pairs)
  with the existing 0.5 Å coincidence tolerance, a spatial grid and a bounded
  report. Ordinary `diagnose` and diffusion both call the shared formatter;
  hydrogen repair behavior is unchanged. **Browser v4 passed** with the explicit
  `cannot safely handle coincident_atoms` refusal and `Sampling was not started`,
  without optimizer NaN. Installed ordinary diagnose independently reported
  `heavy_atom_coincidence.v1`. All browser attempts/evidence are preserved.
- Installed Protenix production template construction passed at **152/512**
  target residues using full input, and **513** using a **48-residue motif**.
  Explicit oversized `full` rejects; incomplete full input does not silently
  crop. This is input/transport evidence, **not checkpoint inference**.
  Native doctor correctly reports no maintained Protenix profile for darwin/arm64.
- After the coincidence-guard refresh, both installed scientific smokes passed
  again (**28/28** Standard, **14/14** Desktop). The current CLI also executed
  MAFFT auto, MUSCLE 5 and Clustal Omega against existing alternative-engine
  executables, plus exact-version manual Probe electron/nuclear diagnostics
  and the no-Probe fallback control (`final-msa-report.json`,
  `final-probe-report.json`). Probe reported `complete`, not a substituted
  engine; the absent-executable control explicitly reported Python fallback.

### Installed cc94 MPS cohort

The protocol in `installed-mps/protocol.json` declared seeds **7, 19, 42, 101,
202**, one attempt per case/seed, before execution. All **25/25** candidates
published and passed independent bundle reading:

| Public workflow/input | Accepted seeds | Total seconds |
|---|---:|---:|
| Model, 1TIM A73–A78 deletion with bounded input repair | 5/5 | 197.85 |
| Model, 1UBQ two internal four-residue gaps | 5/5 | 115.71 |
| Model, 1UBQ N-terminal three residues | 5/5 | 67.43 |
| Model, 1UBQ C-terminal three residues | 5/5 | 67.72 |
| Multichain Homology, 1TIM with bounded input repair | 5/5 | 189.01 |

Each actual invocation records **500 MPS updates**, fallback disabled. The
two-gap case retains separate `raw/gap-0/` and `raw/gap-1/` 500-update traces
for every seed: the merged 500 count is **per invocation**, not their sum.
Twenty candidates selected generated-only CPU refinement; five terminal
candidates selected raw coordinates. All have empty hard-gate failures and
**0.0 Å fixed-heavy displacement**, relative to the explicitly prepared source
where repair was needed. Whole-source warnings remain visible. All private
runner workspaces now remain retained even on complete, clean success (`INFO`
rather than a spurious warning); failures/partial success retain warning-level
notifications. This extends case-specific installation stress evidence, not
source-disjoint scientific qualification.

### Still open

Conda 26.9.4 again rejected the native upstream package path collisions in
**10.48 s**, before transaction execution, despite deliberately unsafe inherited
settings in the control. Child safety settings were enforced, not bypassed.
The public bioconda API still declares **zero osx-arm64 Probe packages**; the
working manual native build is not managed add-on packaging acceptance.
Actual Protenix checkpoint execution and current CUDA end-to-end acceptance
belong on Linux. Broader cc94 MPS/context scientific qualification and hosted
CI acceptance remain unclaimed.

The first final fast-suite attempt exposed a synthetic Homology ownership
fixture whose case-sensitive `D`/`d` chains had identical coordinates. The new
shared guard correctly refused that collapsed input. The ownership fixture now
places the chains separately without changing its sequence/mask assertions;
dedicated coincidence tests retain bonded/cross-chain collapse coverage. The
original failed test log remains preserved (`final-fast-tests.log`), rather
than weakening production admission. The focused final set passed **66 tests**.
The corrected complete fast run passed **1,773 tests, 3 skipped, 33 deselected**
(`final-fast-tests-v2.log`); two additional detector argument-validation cases
were then included in the passing focused set. Skips require Torch/Biotite in
the separate scientific backend environments; they are not claimed as passes.
The real restrained 1TIM repair regression, strict-target mypy, Ruff, version,
CLI/GUI/profile generators, user/agent documentation and installed-distribution
smokes passed. The source archive excludes local evidence, credentials and Node
dependencies; protein-context fixture checksums remain unchanged.

### Fresh supplied 9CQD and 5J13 requests

`old-friends/protocol.json` declared the original tracked PDB/FASTA pairs,
seeds **7, 19, 42, 101, 202**, default terminal inclusion and one attempt per
case/seed before execution. No archived request was edited, no source overlay
used, and input SHA256 values were verified unchanged afterward.

**9CQD:** the new public request was admitted with target lengths **M/W/X =
27/220/221**. It requested four gaps: W `[0,3)`, `[25,28)`, `[217,220)` and
X `[139,141)` (11 residues total). These are not the archived partial-gap
research request or the terminal-excluded W214 request. Source assessment
reported **0 errors, 7 warnings** and did not repair the input. The run took
**569.37 s**; all five seeds actually sampled, with four retained gap traces
per seed showing **500 MPS updates each**, fallback disabled. Independent
bundle reading accepted **4/5 seeds: 7, 19, 42, 202**; normal ranking selected
**42**. Seed **101** was rejected after completed sampling/refinement with
`geometry-gated-refinement-failed`, `generated-or-junction-bond-angle`, and
`anchor-repair-expanded-geometry`; its raw/refined coordinates, result and
traces remain available. Exit 0 reflects a successful accepted subset, not 5/5.

Unlike the preceding controlled 25-candidate cohort, **all four accepted 9CQD
variants selected generated-and-anchors refinement**. Published maximum
observed-heavy displacement was **1.082 Å** (seed 42); only explicitly
authorized anchor movement is allowed by independent publication validation.
Recorded per-anchor Cα displacement remained below 0.5 Å and heavy displacement
below 1.25 Å. This is not a 0.0 Å fixed-coordinate claim for these anchors.
Warnings and the experimental terminal/CYS qualification notices remain visible.
Evidence: `old-friends/9cqd_split/{bundle,work,independent-audit.json,command.log}`.

**5J13:** the unchanged supplied file was refused at public admission in
**0.24 s**: corrected Protpardelle requires canonical protein without nonwater
heterogens. The actual nonwater residues are **SO4 B301/B302/B303** (15 atoms);
they were not silently stripped. No seed sampled, so this is **not a 0/5
scientific result**. An additional ordinary installed `diagnose` run retained
the source report: 255 errors, 24 warnings and 3 informational findings,
including 20 missing-heavy-atom residue findings and 219 hydrogen-position
findings. This uses ordinary diagnose semantics, not gap-aware diffusion input
admission; it does not establish that removing sulfate alone would produce a
safe/admissible structure. A separately authorized and clearly labeled prepared
protein-only derivative is a distinct future case, not a replacement for this
record. Evidence: `old-friends/5j13_selected/` (including the original refusal,
ordinary diagnose and exact input-class audit). All original fixtures remain
unchanged. These case outcomes do not close broader scientific qualification.

## Linux continuation after native finalization (2026-10-29/30)

The initial publication failure described below is historical. The subsequent
[publication compatibility fix](#project-filesystem-publication-compatibility-fix-2026-10-30)
now permits publication on the same project filesystem and has fresh public
workflow evidence. The interrupted supplied 9CQD cohort remains incomplete.

Source commit: `3309bd81e979d2b165633cda8756af26ade862ae`, branch
`feature/installation-profiles`. Execution used a separate detached source
checkout below `.artifacts/installation-profiles/linux-continuation-3309bd8/source`.
Host: Linux x86_64, AMD EPYC 7742, NVIDIA A100-SXM4-40GB, driver **535.104.05**
(CUDA driver API **12.2**), micromamba **2.9.0**. Actual runtime kernels passed;
the driver API version is not the CUDA runtime version.

All outcome-bearing evidence is retained **inside the project**, as requested:

- `tmp/installation-profiles-3309bd8/`: frozen protocols, controlled inputs,
  per-command logs, complete private workspaces, independently checked candidate
  PDBs/traces, `independent-audit.json` per case and aggregate `summary.json`;
- `.artifacts/installation-profiles/linux-continuation-3309bd8/`: source checkout,
  new owned Standard/ML environments, helper evidence, tests and initial Protenix
  attempt. The initial Protenix candidate also has a convenient audited copy in
  `tmp/installation-profiles-3309bd8/protenix-auto-initial/`.

These are ordinary project directories, not the system `/tmp`; their contents
are local evidence and intentionally excluded from Git. An initial detached
checkout was created in `/tmp/opencode`, but no outcome-bearing run from this
increment depends on that checkout. Old environments/checkpoints were retained.
The existing `dvbfixer` environment imported an old editable worktree and was
therefore not used as installed-package acceptance evidence.

### New installations and installed science

- New Standard imports resolve to its own `site-packages`, with `PYTHONPATH` and
  `PYTHONHOME` unset: **Python 3.13.18**, **OpenMM 8.6.1**, **PySCF 2.14.0**.
  Protected direct-bootstrap license configuration passed real MODELLER model
  and Homology checks; profile doctor reported `standard: ready`.
- Installed extended scientific smoke: **25/25 passed**, including both prepare
  backends, protonation, `minimize --no-solvent`, GAFF/OpenFF, PySCF RESP and
  ethanol RESP export, formats/trajectories, cluster/Plotly, topology/ACPYPE,
  CHARMM, custom FF, MAFFT, MODELLER and Homology. Optional GROMACS, PSI4 and MSA
  alternatives were not supplied to this new prefix; this count does not replay
  their earlier separate acceptance or a new Desktop/manager matrix.
- New Protenix managed runtime passed source/asset identities, **PyTorch 2.13.0 /
  CUDA 12.9** real tensor checks, Reference kernel and canonical preflight. It
  retains sampling/refinement in its frozen combined backend process.
- New cc94 runtime `cc94-v2` passed real **Python 3.12 / PyTorch 2.6.0 / CUDA
  12.4** kernels, separate **Python 3.13 / OpenMM 8.6.1 / Reference** refiner,
  source/asset checks and canonical preflight. Checkpoint identity remains
  `a58c009ead22b6f79ef61d44e258b432b0fada1dfd8d44cb7a137d7459038a15`.

Retained setup failures are not readiness passes: the first Standard install
was interrupted by a 120-second outer tool timeout; an owned-prefix resume
installed the package but its outer micromamba bootstrap did not forward the
license FD. A direct invocation of the new Python completed protected
configuration. The first cc94 root failed cloning because the local source path
was misspelled; `cc94-v2` used the existing `protpardelle-1c-audit` clone as a
read-only Git source and succeeded. No old prefix was replaced.

### CUDA computation and independent validation

Protocols declared seeds before their runs. These are **validated private
candidates, not successfully published bundles**:

| Public workflow/control | Seeds | Independently valid private candidates |
|---|---|---:|
| cc94 model, 1UBQ two internal four-residue gaps | 101, 102, 103, 104, 105 | 5/5 |
| cc94 model, 1UBQ N-terminal three residues | 101 | 1/1 |
| cc94 model, 1UBQ C-terminal three residues | 101 | 1/1 |
| cc94 model, 1TIM author A73–A78 deletion and bounded source repair | 101 | 1/1 |
| cc94 multichain Homology, same 1TIM source/deletion | 101 | 1/1 |
| Protenix model, 36hb, auto/full input | 101 | 1/1 |
| Protenix model, 36hb, explicit local motif | 101 | 1/1 |
| Protenix Homology, 36hb, auto/full input | 101 | 1/1 |
| Protenix model, translated-copy size control, auto, 152 targets | 101 | 1/1 |
| Protenix model, translated-copy size control, auto, 512 targets | 101 | 1/1 |
| Protenix model, translated-copy size control, auto, 513 targets | 101 | 1/1 |

Total: **15/15** completed candidates passed a fresh installed
`validate_runner_result` audit, with zero detectable D-Cα, empty hard-gate
failures and **0.0 Å fixed-heavy displacement** relative to the explicit prepared
source where applicable. CUDA traces have fallback disabled. Protenix used
200 updates; cc94 used 500 updates per invocation (per gap, not summed into one
sampler count). Complete workspaces retain per-gap evidence.

The translated-copy controls follow `tests/test_protenix_auto_input.py`; they are
size/adapter controls, not biological assemblies or source-disjoint scientific
qualification. Actual checkpoint inference used full input at **152/512** target
residues and a **48-residue motif** at **513**. Explicit `full` at 513 and an
otherwise identical 152-residue input missing B/20:CB refused before sampling,
without silent chemistry-dependent cropping. This is now CUDA execution evidence,
beyond the native Mac's earlier input-construction-only checks.

Independent source-repair audits also recomputed the two 1TIM repairs: 47
observed residues selected, all movement limits met, no remaining severe input
findings; 109 warnings remain visible. All original tracked and controlled input
hashes remain unchanged. The first Homology harness plan incorrectly equated
author residue numbering with SEQRES columns and failed before sampling. Its log
and protocol remain; a separately declared v2 plan uses the shared affine
sequence placement for the same deleted atoms and seed. No scientific gate or
source mask was changed to obtain the completed v2 result.

The original tracked 5J13 PDB/FASTA was refused for nonwater heterogens in 0.4 s,
with no runner workspace or sampling. This preserves the SO4 admission refusal;
it is not a scientific seed rejection or a protein-only derivative.

### Publication failure and interrupted supplied 9CQD cohort

**Zero bundles were published.** Both the project root and its new `tmp/`
subdirectory are on the same Ceph filesystem. A direct empty-directory probe
confirmed `renameat2(..., RENAME_NOREPLACE)` returns **EINVAL (22)**. Each
completed inference reached the unchanged fail-closed publication gate and
retained diagnostics; a new subdirectory cannot change filesystem semantics.
The subprocess-recorded model/Homology/refusal controls returned **1**, as
expected. The first shell-only invocation did not retain an independently
recorded return code; no exit-0-error claim is made from that wrapper output.

Consequently public bundle consumption, batch/resume, ZBS diffusion handoff and
the frozen **36hb seed-7 repeated-bundle hardware verifier remain blocked/not
executed** for this increment. Private validator success is not a substitute for
`read_published_bundle(...).selected()`. Atomic publication was not bypassed and
no rejected/partial structure was promoted into a public output.

A fresh supplied 9CQD full-PDB/full-FASTA protocol declared seeds **101–105**.
Admission preserved targets M/W/X **27/220/221**, all four gaps and **0 input
errors / 7 warnings**. Seed 101 produced its raw result and entered Reference
refinement, but the **600-second outer tool timeout** terminated the controlling
CLI. The remaining orphan sampler/refiner process group was explicitly
terminated, with raw arrays/results and workdir retained. This is an
**infrastructure-interrupted attempt**, not a final geometry rejection; seeds
102–105 did not execute. The protocol and timeout evidence are retained at
`tmp/installation-profiles-3309bd8/cc94-supplied-9cqd/` and in `summary.json`.
The complete supplied cohort remains pending with suitable orchestration limits.

### Repository checks and remaining gates

- Focused tests: **186 passed**.
- Full fast suite on the clean source: **1771 passed, 7 skipped, 33 deselected**;
  skips are one missing Biotite and six unavailable GROMACS cases. Two SWIG
  deprecation warnings remain, plus the interpreter-exit SWIG notice.
- Ruff, strict-target mypy, versions and agent/user documentation checks passed.
- Python-3.11 CLI/GUI generators and installation-profile generator passed after
  providing OpenMM/PDBFixer to the isolated generator environment. The initial
  missing-OpenMM generator attempt was an environment setup failure, not drift.
- GUI source was not modified in this increment; no new installed Desktop or
  browser acceptance is claimed.

The current Linux publication filesystem gate, complete supplied cohort and
public workflow/repeat verification remain open, alongside the previously
recorded native conda/Probe packaging and broader scientific qualification gates.
The cross-platform installation plan is **not complete**.

## Project-filesystem publication compatibility fix (2026-10-30)

Implemented on top of `3309bd8` in `model/diffusion/bundle_commit.py` and the
pipeline's publication dispatch. Unsupported directory-rename flags
(`EINVAL`, `ENOSYS`, `ENOTSUP`, or a missing platform rename function) now select
an exclusive **manifest-last commit**. Other I/O/permission errors still fail.

The publisher atomically reserves the output name with `mkdir`, transfers only
private flat regular staged files via exclusive hard links, removes their
staging links, fsyncs all payload entries, and atomically links `bundle.json`
last without replacing an existing entry. Completed artifacts retain one link,
schema 2 and their existing digest/provenance/scientific validation rules.
Filesystems supporting the original whole-directory no-replace rename still
use that path.

This makes the readiness boundary explicit: a directory can be visible while
being populated, but **no completed bundle exists without its index**. Both the
Python reader and GUI transport reader require that index and verify linked
payloads. A crash before commit leaves an unready reservation that is refused
by readers and later publishers, never silently adopted. Handled failures roll
back only owned entries, removing the index before payloads and preserving
unrelated entries. Publication does not substitute ordinary overwrite-capable
directory rename, disable checks, or move outputs outside the project.

### Actual project CephFS evidence

Evidence root: `tmp/installation-profiles-3309bd8/publication-fix/` inside the
project. The existing new Standard prefix was refreshed with a built package,
without changing scientific dependencies; installed checks unset `PYTHONPATH`
and `PYTHONHOME`. Backend environments/checkpoints were unchanged.

1. **Recovery, without resampling:** all **15** previously validated retained
   candidates were independently validated again, published as **11 bundles**,
   then read with the real `read_published_bundle` in request/source-bound resume
   mode. This includes both source-repaired 1TIM cases and Protenix at
   152/512/513 targets. `recovered-bundles.json` records exact paths and seeds;
   `recovery.log` records outcomes. The earlier failed attempts remain retained.
2. **Fresh public checkpoint execution:** `scripts/smoke_diffusion_workflows.py`
   passed on the project filesystem with the newly installed CLI and managed
   Protenix runtime. Two independent same-input batch runs published identical
   selected PDBs; resume verified both existing bundles; `zbs --no-solvent`
   completed and preserved the modeled `.dat` handoff; mosaic Homology sampled,
   published and passed bundle validation. The source remained unchanged.
   Evidence: `workflow-smoke/workflow-smoke.json` and all command/private logs.
   The selected batch PDB SHA256 is
   `60d80ca655cdffd7c9c17eadc1fc86a563871e5832ce35d6910fc34b88ec31c7`;
   Homology is `d8e266680726294beb2aaecc963e8b6fd6b44e3d0d93442dcd0a136fd2a80ee0`.
   These current auto/full-input controls are not the historical fixed-digest
   36hb hardware-repeat fixture, whose verifier remains a separate gate.
3. **GUI transport:** the unchanged real TypeScript
   `readDiffusionPublication` reader successfully consumed the recovered and
   fresh bundles (`gui-transport-check.json`). This verifies transport
   compatibility, not a new browser/installed-Desktop acceptance run.
4. **Tests:** initial publication/workflow tests **44 passed**; the full fast
   suite passed **1786 tests, 7 skipped, 33 deselected**. New coverage includes
   unsupported rename flags, incomplete-reader refusal, empty-destination races,
   two competing publishers, transfer/fsync failure rollback, preservation of
   unrelated entries, and an actual subprocess crash before index commit.
   Ruff and strict-target mypy passed. Scientific admission/geometry thresholds
   did not change.

The project-filesystem publication blocker is resolved for these executed
paths. The interrupted supplied 9CQD cohort, frozen hardware-repeat fixture,
native package gates and broader scientific qualification remain open.

## Supplied-case continuation and hardware replay audit (2026-10-30)

Source baseline: `0405d81`. Evidence is retained inside the project at
`tmp/installation-profiles-0405d81/`; `protocol.json`, per-job `execution.json`,
command logs, complete private workspaces, published bundles, independent audits
and `summary.json` are available there. New runs use the installed Standard and
the managed runtimes from the preceding Linux acceptance, with no source overlay
in Standard and no change to scientific thresholds or the Reference refiner.

### Completion of the five declared 9CQD seeds

The previous seed-101 attempt remains untouched. Its complete private workspace
was copied to `seed-101-continuation/workspace/`; the saved raw result, four gap
sampling traces and single-seed request were digest-checked before and after
continuation. Only refinement was rerun, into a new `refined-continuation/`
directory. No second neural sample was generated for seed 101. Seeds 102–105
then executed their first public single-seed `model` requests with the original
tracked PDB and full companion FASTA. They retain targets M/W/X **27/220/221**,
four gaps and all terminal tails, without an archived partial-request bypass.

The supervisor used a **10,800-second runner limit** and **11,100-second job
limit**, detached from individual terminal-tool timeouts. GPU sampling was
sequential, with at most three live single-thread Reference refiners; temporary
scheduler pauses allowed hardware replay to use the GPU independently while
already-running CPU refiners continued. Every job and both supervisors finished;
there are no remaining acceptance jobs from this cohort.

| Seed | Workflow | Outcome | Observed-heavy max displacement | Job seconds |
|---:|---|---|---:|---:|
| 101 | Saved raw-result refinement continuation | accepted | 0.994 Å | 993.18 |
| 102 | Fresh public model | validation_failed | 1.056 Å | 1278.87 |
| 103 | Fresh public model | accepted | 0.901 Å | 737.64 |
| 104 | Fresh public model | accepted | 1.028 Å | 786.17 |
| 105 | Fresh public model | validation_failed | 1.155 Å | 806.59 |

All **five** completed bundles passed digest/request-bound consumption and fresh
independent validation of their recorded status. **3/5** candidates are accepted
(101, 103, 104); rejected bundles 102 and 105 are not selectable. Actual GUI
transport verification likewise returned no primary structure for the two
rejections (`gui-transport-check.json`). No failed seed was omitted or rerun to
obtain a favorable structure.

All five retain four per-gap CUDA traces with **500 updates per invocation**,
fallback disabled and zero detectable D-Cα. All selected the
`generated-and-anchors` stage, so this is **not a zero-observed-motion result**:
movement within the authorized anchor mask is independently audited. Every
observed atom outside the reported refinement residues stayed exactly fixed.
The original PDB/FASTA hashes and seed-101 original raw artifacts remain unchanged.

Seed 102 failed `geometry-gated-refinement-failed`,
`generated-or-junction-ramachandran` and `anchor-repair-expanded-geometry`.
Seed 105 additionally failed `anchor-repair-per-anchor-movement-limit`.
Complete raw/processed coordinates, reports and traces remain retained in their
explicitly rejected bundles/workspaces. These Linux outcomes do not replace the
different native-MPS seed cohort or establish broader scientific qualification.

### Hardware replay outcomes: current repeatability is not the old frozen hash

The original strict `scripts/check_protenix_hardware_acceptance.py` and its
expected PDB SHA256 `664e65f90a9a7ec6abdd762008b6f28fe20d2e47eaab576bdf1830a246b61ea4`
were **not changed**. Its historical frozen gate remains **not passed**:

1. The first harness incorrectly selected the early 36hb materialization
   (gap `[38,43)`). Both runs completed with identical PDBs and passed geometry,
   but this was the wrong input for the frozen hash. Logs/protocol remain as
   an explicit harness failure, not a frozen scientific comparison.
2. Corrected input selection used cohort-500-v5, gap `[39,44)`, normalized-input
   SHA256 `8aa3012f6cfd18ba142ef1086a0c857baa7248fc8045c5edfd77808b427dd9e7`.
   Two independent current-backend runs produced identical accepted PDBs with
   SHA256 `6348985007d25877e5ce3e5a0e518a7c28e9ff6715da0b66c013a58869414106`.
   Their request hashes match the archived baseline, CUDA/callback evidence is
   complete (200 updates, zero projection error), and fixed displacement/D-Cα
   counts are zero. The strict frozen verifier correctly failed the old file
   identity comparison.
3. An offline audit established that all **896 raw atom coordinates** match
   the archived raw output exactly. The historic final baseline used refinement
   **v6**, while the current runtime uses **v8**; final coordinates differ for
   48 generated atoms, up to **0.325 Å**, and current output also carries updated
   sequence metadata. This is not merely a file-header difference or a sampling
   repeatability failure. See `hardware-independent-audit.json`.
4. A separately declared isolated historical-source replay at
   `1b49bc57696e45ef4a47b06bbf35ee1fcc57b6b9` (v6) also did not establish the
   frozen gate: its first candidate passed geometry but had another PDB digest;
   the second attempt was stopped by the workspace monitor before a result.
   The v6 label alone does not prove the full historical shared-source identity.
   The isolated checkout/private launcher and all failed evidence remain; no
   current production source, launcher, model, threshold or checkpoint was
   downgraded. Further replay requires the exact frozen execution source chain,
   or an explicitly versioned current-v8 acceptance protocol; substituting an
   observed digest into the old checker is not an acceptable resolution.

### Runtime guard hardening exposed by acceptance

The aborted historical repeat reported “workspace size limit” although its
retained regular files total only **202,275 bytes**. The monitor previously
treated *every* `OSError` during a live directory scan as exceeding the byte
limit, including ordinary `FileNotFoundError` when a runner removes a temporary
file/directory between `scandir` and `stat`. The original abort did not retain
the underlying OS exception, so its precise transient cause is not claimed as
proven from the postmortem size alone.

`runner.py::_workspace_size_exceeds` now tolerates only disappearing **child**
entries/directories and continues counting the rest. Missing/unreadable roots,
permission/other I/O errors and actual oversized workspaces still fail closed.
The byte limit is unchanged. Regression tests exercise disappeared entries,
disappeared child directories, continued counting of oversized files and root/
child access failures; the existing real oversized-runner termination test
still passes.

The updated package was installed into the acceptance Standard without changing
dependencies. A fresh public current-v8 run then published, passed the real
bundle reader and reproduced both preceding current-v8 PDBs exactly
(`monitor-fixed-current-v8/verification.json`). It is diagnostic current-version
repeat evidence, not a replacement frozen-v6 baseline.

Verification: **78 focused tests passed**; full fast suite **1793 passed,
7 skipped, 33 deselected** (Biotite/GROMACS availability skips), Ruff and strict
mypy passed. The new source change affects operational workspace accounting,
not geometry, sampling, refinement or publication policy. Native macOS package
collisions/Probe packaging, broader scientific qualification and exact frozen
historical hardware replay remain open.

## Final Linux acceptance before return to Mac (2026-10-30)

Portable aggregate: [Linux finalization JSON](research/linux-installation-finalization-2026-10-30.json).
Full evidence remains under **`tmp/installation-profiles-linux-final/` inside the
project**; previous `3309bd8`/`0405d81` attempts, original inputs, environments
and failed logs are preserved. This finalization uses the preceding `0405d81`
baseline plus the workspace scan fix and acceptance-tool/fixture changes
recorded here. The next published commit carries those changes; it is not a
release-version bump or a redefinition of scientific qualification.

### New installed Desktop and executable/add-on matrix

- A **new micromamba 2.9.0 Desktop prefix** completed protected MODELLER setup,
  all **14/14** licensed scientific smokes, profile doctor, fresh GUI lockfile
  installation, typecheck, client/server build and standalone HTTP smoke.
- Its installed GUI passed **407 tests in 57 files**. A real Chrome Headless
  Shell **156.0.8078.4** session created a workspace through the UI, imported the
  tracked 1UBQ structure through the browser API, completed an actual installed
  `diagnose` job and loaded the structure into the Mol* viewport through the UI.
  No browser runtime errors occurred. Evidence:
  `browser-attempt-v4/browser-acceptance.json` and screenshot/workspace records.
- The first browser attempt could not load host shared libraries. An isolated
  project-local browser-libraries prefix supplied them without modifying system
  packages. Browser v2 passed workspace/job transport; v3 used an incorrect
  filename selector; v4 selected the actual displayed basename and verified
  viewport loading. All attempts remain retained. This is a Linux browser pass,
  not native macOS proof.
- Installed extended science passed **30/30** checks with actual optional
  executables: GROMACS **2026.3**, MUSCLE **5**, Clustal Omega and a new separate
  managed **PSI4 1.9.1 / psiresp 0.4.2** runtime (Pydantic <2, Libint 2.9.0).
  The separate PSI4 RESP integration test passed too. Normal PySCF/RESP remains
  included and executed in Desktop, as in Standard.
- GROMACS generated the boxed coordinates and TPR with `grompp`; the application
  reader loaded that **TPR with matching coordinates, 1231 atoms**, successfully.
  This route did not replace an unsupported TPR silently with another topology.
- Exact **Probe 2.18** actually executed electron and nuclear diagnostics, both
  reporting `complete`/`molprobity-probe`. A PATH-isolated absent-Probe control
  explicitly selected the Python fallback. `doctor --profile desktop --addon
  msa-alternatives --addon probe` passed using these external executables. The
  optional executable prefix is separate acceptance tooling, not a claim that
  this Desktop install automatically included unselected add-ons.

Installed imports were independently checked in `site-packages`; no checkout
`src` overlay was supplied to scientific runs. `pip check` **still fails** for
the already documented upstream AmberTools auxiliary metadata (proprep/
packmol-memgen missing pdb2pqr, auxiliary NumPy <2 and Biopython <1.86 declarations).
Executed science passing does not turn this metadata issue green or justify
loosening the project's NumPy/Python constraints.

### Remaining CUDA workflow paths executed

`scripts/smoke_diffusion_workflows.py` now accepts explicit
`--diffusion-model`, `--case` and `--seed`. Its defaults preserve the existing
Protenix/36hb/seed-7 smoke; the reported profile is read from the actual selected
sampler trace rather than hard-coded. It validates its one-chain/one-gap input
before creating the work directory. Scientific dispatch remains in the normal
public CLI.

A declared **cc94 1UBQ A10–A13 deletion / seed 101** control passed two independent
batch publications with identical PDBs, verified resume, `zbs --no-solvent` with
preserved modeled `.dat` identities, and independently validated mosaic Homology.
Evidence: `cc94-workflows/workflow-smoke.json`; batch PDB SHA256
`f0fb1724598676ae1bd25cba1a931181a0ddc95558c1e8cdf5aea22692c2a36c`,
Homology `10a60c705fcc98039c55675f3321959fc8853e91310a54e3c0a7e3befcce9e88`.
This complements the earlier completed Protenix batch/resume/ZBS/Homology run.

Protenix **model and Homology** also executed the controlled **1TIM A73–A78 /
seed 101** source-repair case. Both repaired the observed source before selecting
the **494-residue full-input axis**, completed actual CUDA inference, published
accepted bundles and passed fresh independent candidate/bundle reading.
Source-repair evidence is bound to the final request. Fixed-heavy displacement
is **0.0 Å relative to the explicitly prepared source**, not the unrepaired
original. The original and local repair reports/warnings remain available.
Evidence: `protenix-repair-results.json` and both complete private workspaces.

### Exact historical hardware source replay: retained failure, not a Mac task

Git history identified **`05c16c9`**, the actual commit introducing the strict
hardware verifier and expected PDB digest. A new detached read-only checkout and
contained launcher replayed that historical source using the frozen 36hb
cohort-500-v5 gap `[39,44)`, seed 7, unchanged engine/checkpoint/profile lock and
the frozen combined backend process layout. Only this explicitly declared
historical replay used the old source; current production runtimes remain v8.

Both historical runs finished, published and passed current independent
request-bound reading; their PDBs match each other at SHA256
`bea7eeaa0ae72f06ecab536bb8c5bbdf79ae85afb15f840a42c653be6049e501`,
but **not** the unchanged archived expected
`664e65f90a9a7ec6abdd762008b6f28fe20d2e47eaab576bdf1830a246b61ea4`.
The strict verifier therefore still **fails**. This eliminates selection of an
arbitrary v6-label source as the sole explanation; a historical source commit
and available package lock do not reconstruct every unrecorded shared/native
runtime detail. No expected hash was edited after execution, no production
algorithm was downgraded, and old evidence was not replaced. The current-v8
repeatability/geometry passes remain separate.

Exact historical byte identity remains an explicitly open **Linux forensic
gate**, not something the Mac can demonstrate or an excuse to repeat CUDA on
macOS. Its resolution requires recovering the missing historical execution
identity/evidence or a separately approved, preregistered versioned baseline;
the observed hashes here do not retroactively become such a baseline.

### Final verification and retained attempts

- Complete **installed** Python suite, including all **33 slow tests**:
  **1831 passed, 2 skipped, 12 warnings** (`installed-full-suite-v3.log` and
  JUnit XML). The two skips require Torch, intentionally absent from Standard/
  Desktop; both exact controls separately passed in the ML-capable test env.
  PSI4/RESP, Biotite reference metrics and GROMACS checks executed in this matrix.
- The first sequential complete suite hit a **2400-second external timeout**
  near the slow ZBS tail, without a summary; it is not a complete pass. A
  detached, bounded 3-worker supervisor retained per-test outcomes. Its first
  completed run had **1829 passed, 3 skipped, 1 failed**: a provenance fixture
  incorrectly assumed pytest temp paths were outside Git. Production correctly
  found the parent checkout commit. The fixture now supplies a nonexistent
  repository path; its focused recheck and the final complete run passed.
  No production provenance rule was changed to satisfy that assertion.
- Warnings retained: SWIG deprecations, deposited CONECT-to-TER notices and two
  NumPy-extension size warnings during slow structural tests. They were not
  suppressed or reclassified as errors. All structural assertions completed.
- Ruff, strict-target mypy, Python-3.11 CLI/GUI/profile generators, version/
  user/agent docs checks and `git diff --check` passed. Tracked structural
  fixture checksums remain unchanged. Wheel and sdist built successfully and
  were inspected to exclude local evidence, credentials, node_modules and
  temporary/environment trees; the new publication helper is included.

The available Linux installation/science/workflow regression matrix is now
executed. **The whole cross-platform plan is not complete**: the old strict
hardware digest remains failed, upstream auxiliary pip metadata remains
inconsistent, and native conda/Probe packaging plus broader MPS qualification
remain open. Continue on the Mac using the current handoff rather than treating
historical bootstrap “not-tested” paragraphs as the latest outcome.

## Native Mac return after Linux finalization (2026-10-30)

Pulled `feature/installation-profiles` fast-forward to **`8f0cd33`** on native
arm64 / macOS **26.6.2**, initially clean. Evidence is project-local at
`tmp/installation-profiles-mac-return-8f0cd33/`; Linux runtime trees were not
assumed available. Micromamba **2.9.0** and conda **26.9.4** were rechecked.
The shell's mamba resolves to an older **1.5.9** executable; this continuation
does not claim a new mamba 2.10.0 acceptance run.

### Publication/runner regression and existing-prefix findings

The five focused runner/bundle-commit/pipeline/workflow/provenance files passed
**87 tests** against source and again against wheel-installed Python modules.
These native checks include normal Darwin exclusive publication, forced
unsupported-rename manifest-last commit, competing publishers, interruption,
incomplete reservations and quota-scan failure semantics. Version, generated
CLI/GUI/profile specs, user/agent docs and Ruff checks passed.

Existing Standard/Desktop Python packages were refreshed from the new wheel
without dependencies; verified prior-package tarballs are retained. Their
scientific rechecks exposed missing environment data, **not a successful
installation acceptance**: each old prefix lacks **3,030 conda-owned files in
45 packages**, including ANARCI's pressed HMM indices, `hmmpress`, and
`pdbfixer/soft.xml`. The cause is unestablished. A proposed native `hmmpress`
repair stopped before execution because that executable itself was absent;
the original HMM was backed up and was not modified. No broad repair, deletion,
clobber policy or silent environment substitution was applied.

The old-prefix fast attempt retained **1,787 passed / 10 failed / 3 skipped**
with 33 slow tests deselected. Missing package data caused scientific failures;
the outside-checkout harness working directory also broke a relative-path test.
Keep this failure record separate from the new-prefix run, rather than weaken
tests or count it as a pass. Existing installed Desktop independently passed
typecheck, tests, build/standalone smoke and the retained Chrome import/diagnose/
visible-failure protocol. These GUI passes do not establish old scientific
environment integrity.

### Fresh native Standard and MPS workflows

A new owned micromamba Standard prefix, sibling
`dvbfixer-apple-return-standard-8f0cd33`, passed conda-file presence audit
(**zero missing owned files**), protected-fd MODELLER configuration and the
installer's licensed smoke. Its independent extended science run passed
**28/28** checks, including PySCF/RESP, separate PSI4/RESP, GROMACS,
MODELLER/Homology and both prepare backends. Profile doctor reports **ready**;
actual cc94 MPS preflight passed with Torch **2.6.0**, fallback disabled and
separate OpenMM CPU refinement. Installed distribution smoke and the real
restrained 1TIM repair regression passed. Upstream auxiliary `pip check`
conflicts still return 1 and are not fixed by this clean installation.

`mps-protocol.json` froze the existing 1UBQ C-terminal-three-residue case and
seed **7**, retaining historical originals and their digests. The first MPS
workflow attempt passed two identical independent batch outputs and resume,
then failed ZBS on the old prefix's missing `soft.xml`; Homology did not execute
in that attempt. The distinct new-prefix **v2 passed public batch, resume,
`zbs --no-solvent` and mosaic Homology**, with preserved modeled-atom `.dat`
handoff. All three bundles were independently read; each trace records the
actual `protpardelle-1c-cc94-mps-v1` profile, **500 MPS updates**, fallback
disabled. Original input, both workflow attempts, private diagnostics and
selected results remain retained. This is a native workflow regression, not
source-disjoint qualification or a replay of the failed historical CUDA hash.

### Upstream gates rechecked

Conda's first attempt hit missing files in its historical cache and is not a
clean upstream retest. A second attempt with a **new project-local package
cache and new prefix** reached transaction verification and failed closed on
the existing gettext/libintl dylib and Tk/X11 header ownership collisions.
No transaction execution or `--clobber` bypass is claimed. The public bioconda
API still declares **zero osx-arm64 Probe files**. These packaging gates remain
open; historical manual native Probe evidence remains a separate alternative.

### Completed fresh-prefix continuation

The new installed Standard fast suite passed **1,797 tests, 3 skipped,
33 deselected** (`new-installed-fast.log`, **184.56 s**). The missing Torch and
Biotite modules are confined to the explicitly recorded skips, not passes;
the separate MPS runtime executed real inference. This run used the checkout
as working directory with `PYTHONPATH`/`PYTHONHOME` unset and installed DVBFixer
imports, resolving the earlier harness-relative-path issue without changing a
test. The real 1TIM bounded-input-repair slow regression also passed; the other
32 slow cases were not rerun in this native increment.

A distinct new micromamba **Desktop** prefix,
`dvbfixer-apple-return-desktop-8f0cd33`, also passed the protected-license
installer smoke and conda-file presence audit (**zero missing files**). Its
installed GUI passed **407 tests in 57 files**, typecheck, client/server build
and standalone smoke. A fresh **Chrome 154.0.8037.98** run against this new
Desktop passed real workspace creation/import/diagnose, Protenix context-choice
visibility, hidden Protpardelle context override, and the explicit collapsed-input
refusal before sampling; no page errors occurred. Profile doctor returned 0.
The earlier existing-Desktop browser/build checks and all failed science attempts
remain separate evidence, not overwritten by this result.

The existing pinned manual **Probe 2.18** binary was also exercised from the
new Standard: electron and nuclear modes reported `molprobity-probe`, complete,
without fallback. The absent-executable control explicitly reported the Python
steric fallback. This does not close the missing native managed-package gate.

Portable summary:
[native return results](research/mac-installation-return-2026-10-30.json).
All evidence and original inputs remain retained locally; no credentials or
runtime environments are included in Git. The historical Linux PDB-hash gate
remains failed, auxiliary pip metadata remains inconsistent, native conda and
managed Probe remain blocked, and wider source-disjoint MPS qualification is
still pending. **Whole-plan completion is not claimed.**
