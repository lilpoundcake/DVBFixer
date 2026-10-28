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
adapters. The old `environment.yml` remains the compatibility manual route
until replacement acceptance is demonstrated.

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
- **WP4 partial groundwork:** profile-aware doctor imports required modules in
  isolated children, checks executable presence, ANARCI pressed database,
  Open Babel/ML isolation and native linear algebra. Full scientific-command
  preflights (PROPKA/Reduce state classification, antibody error vs no domain,
  topology/trajectory/TPR compatibility) remain pending.
- **WP5 partial:** the Linux matrix below is executed evidence, not the complete
  promised structural/format/force-field matrix. Optional runtimes and native
  macOS/MPS remain not-tested.
- **WP6 pending:** CI migration, final manual-path replacement and GUI/MSA moves
  wait for the remaining gates. `environment.yml` is still the compatibility
  input and has not been silently slimmed.

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
