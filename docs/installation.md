# Installation

[User guide](README.md) · [Quick start](quick-start.md) · [Troubleshooting](known-issues.md)

## 1. Get the source

Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
and Git first. The commands below assume a terminal on macOS or Linux; on
Windows use WSL. Check backend availability with `doctor` after installation.

```bash
git clone https://github.com/lilpoundcake/DVBFixer.git
cd DVBFixer
```

## 2. Create the environment

Run these commands from the repository root:

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
```

The default `environment.yml` is generated from **Standard**: the scientific
Python stack and normal CLI external tools, including PySCF/RESP and analysis.
Keep Python in the supported `>=3.11,<3.14` range. Editable installation uses
the checkout directly; keep it in place while using this environment.

## 3. Configure MODELLER

Gap rebuilding and comparative modeling need a MODELLER license. Complete the
[license instructions](#modeller-installation-and-license) before those workflows.
Other tools can be checked without it.

## 4. Verify the installation

```bash
dvbfixer --version
# dvbfixer 0.10.1
dvbfixer doctor
dvbfixer msa --list-engines
```

Read unavailable-backend messages in the capability report. A missing optional
tool does not block every workflow; consult the command you plan to run.

Without activating the environment, prefix a command with
`micromamba run -n dvbfixer`, for example:

```bash
micromamba run -n dvbfixer dvbfixer doctor
```

Continue to the [quick start](quick-start.md) or [GUI guide](gui.md).
The remaining sections cover alternative installations and specific dependencies.

## Standard/Desktop source-checkout helper

The helper and the two main profiles have Linux clean-install evidence with micromamba 2.9.0, mamba
2.10.0 and conda 26.9.4. See the [dependency matrix](installation-dependency-matrix.md)
and [acceptance record](installation-profile-acceptance.md) for implemented scope,
resolved builds and remaining Apple acceptance gates. The manual `environment.yml`
route above is the Standard alias. `environments/legacy-full.yml` retains the
former Desktop-plus-alternative-MSA composition as an explicit compatibility input.

On native Apple Silicon, use micromamba **2.9.0** or mamba **2.10.0**: both
passed licensed Standard/Desktop workflows and installed GUI checks. Micromamba
2.0.2 stalled in native solves on the acceptance host. Conda 26.9.4 currently
has unresolved upstream macOS package-file collisions; the optional Probe 2.18
conda package is unavailable for osx-arm64. Do not bypass conflicts with
`--clobber` or substitute x86 packages. See the
[native installation limitations](known-issues.md#native-apple-installation-profiles).

- **Standard** includes scientific CLI features, MODELLER, ligand/topology and
  antibody tools, MDAnalysis/Plotly, MAFFT and **PySCF/RESP**.
- **Desktop** inherits Standard and adds Node/npm plus a GUI built from its lockfile.
- `--addon msa-alternatives` adds MUSCLE 5 and Clustal Omega;
  `--addon probe` selects the optional supported MolProbity Probe.
- Diffusion, PSI4 and external Gaussian remain separate runtimes/operator
  integrations. [Diffusion-runtime provisioning](installation-runtime-provisioning.md)
  prepares isolated frozen environments and performs checked registration.
  The lower-level launcher-registration script remains usable independently.

Use Python **3.11+** to bootstrap from the repository root. An interactive
invocation is `python3 scripts/install_dvbfixer.py`; for automation:

```bash
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --prefix /path/to/new-standard --non-interactive --skip-modeller

python3 scripts/install_dvbfixer.py --profile desktop --manager micromamba \
  --prefix /path/to/new-desktop --non-interactive --skip-modeller
```

Choose a **new prefix without whitespace**: MAFFT's upstream shell scripts
failed the real spaced-prefix smoke. Source/input paths may contain spaces.
The helper refuses existing environments; `--resume` is limited to installations
it previously recorded. Resume reconciles native dependencies before pip.
`--dry-run` prints composition and the selected manager without creating files.

MODELLER supports hidden interactive key entry, or an open protected file
descriptor for automation. Configure a deferred license without reinstalling:

```bash
# KEY_FILE is a user-owned 0600 file containing only the key, not its value in argv.
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --prefix /path/to/new-standard --non-interactive --resume --configure-only \
  --modeller-key-fd 3 3< "$KEY_FILE"
```

The helper runs real functional smokes and records their commands, results and
versions below `<prefix>/share/dvbfixer/installation-smoke/`. Finish with
`dvbfixer doctor --profile standard` (or `desktop`). Exit codes are **0** for
passed installation prerequisites, **2** for required manual action such as a
deferred MODELLER license, and **1** for failure. Scientific/hardware acceptance
remains scoped separately in the acceptance record. Run Desktop without
activation using:

```bash
micromamba run -p /path/to/new-desktop npm --prefix /path/to/new-desktop/share/dvbfixer/gui start
```

For manual profile installation without the helper, generate/check the specs
with `PYTHONPATH=src python scripts/gen_installation_profiles.py --check`, create a new prefix
from `environments/standard.yml` or `environments/desktop.yml`, then install
DVBFixer with `python -m pip install -e .`. For manual Desktop, run `npm ci` and
`npm run build` in the checkout's `gui/` and start it with `npm start` while that
environment is active. `doctor --profile desktop` recognizes this editable
source build as well as the helper's installed build. These specs use only
**openbabel-wheel** for Open Babel. Do not overlay conda `openbabel`; ACPYPE's
Python metadata requires the wheel. The `obminimize` console entrypoint launches
the wheel's bundled native tool with its plugin paths initialized.

## Wheel installation and bundled force fields

For a non-editable installation, build or install the wheel normally. Release
wheels include the AMBER99SB-ILDN/Lipid21, official GROMACS `v2026.3`
Amber19SB, and CHARMM36 GROMACS force-field trees used by `dvbfixer top`; no
repository checkout or separate GROMACS force-field installation is required. At runtime they resolve under the active
environment's `share/dvbfixer/FF/` directory. Editable installs continue to use
the checkout's top-level `FF/` directory.

The Amber19SB tree is redistributed unmodified under the GROMACS LGPL 2.1
terms. Wheels install its complete byte-verified inventory and corresponding
notice under `share/dvbfixer/`. Source provenance, hashes, citations, and the
license assessment are recorded in
[`provenance/amber19sb-gromacs.md`](provenance/amber19sb-gromacs.md).

## CIF conversion backends

PDBx/mmCIF input uses Gemmi; small-molecule crystallographic CIF uses Open
Babel. Both are installed by `environment.yml`, and pip installations declare
`gemmi` plus the prebuilt `openbabel-wheel` explicitly. Verify the active
environment with `python -c "import gemmi; from openbabel import openbabel"`.
If only PDB input is used, the conversion layer is not invoked.

## AmberTools and MolProbity Reduce

The `tleap-reduce` preparation backend requires both `tleap` (AmberTools) and
the `reduce` executable (MolProbity). `environment.yml` declares `ambertools`
and `reduce` separately; an existing environment containing `tleap` is not proof
that Reduce is installed. Repair an incomplete environment with:

```bash
micromamba install -n dvbfixer -c conda-forge -c bioconda ambertools reduce
micromamba run -n dvbfixer dvbfixer doctor
```

Missing Reduce is fatal for the explicitly selected `tleap-reduce` backend.
Legacy preparation reports the missing optional Reduce step and falls back to
its documented input-H/tautomer policy; that warning does not establish that the
external tool ran successfully.

## Multiple-sequence alignment executables

`dvbfixer msa`, `dvbfixer homology --template-plan`, and the GUI Homology
workspace run an external alignment program. Standard and Desktop install
MAFFT. Select `--addon msa-alternatives` for MUSCLE 5 and Clustal Omega, or use
the explicit compatibility input `environments/legacy-full.yml` for all three:

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
dvbfixer msa --list-engines
```

For an existing conda/micromamba environment, install one or more engines
from conda-forge/bioconda. MAFFT is the default and is sufficient by itself:

```bash
micromamba install -n dvbfixer -c conda-forge -c bioconda mafft
# Optional alternatives:
micromamba install -n dvbfixer -c conda-forge -c bioconda "muscle>=5,<6" clustalo
```

The commands must be available on `PATH` under these exact names:

| Engine | Required executable | Check |
|---|---|---|
| MAFFT | `mafft` | `mafft --version` |
| MUSCLE 5 | `muscle` | `muscle -version` |
| Clustal Omega | `clustalo` | `clustalo --version` |

On macOS, Homebrew users can install the native packages with
`brew install mafft muscle clustal-omega`. On Debian/Ubuntu, distribution
packages are commonly available as `mafft`, `muscle`, and `clustalo`, but
verify that the `muscle` package is version 5: MUSCLE 3 uses a different CLI
and is not supported. Native Windows users should use WSL with the Linux
instructions, or download upstream binaries and put their directory on
`PATH`. When installing a standalone MUSCLE binary, rename it to `muscle`
(`muscle.exe` on Windows); rename a standalone Clustal Omega binary to
`clustalo` (`clustalo.exe` on Windows).

After any standalone installation, restart the GUI/dev server so it inherits
the updated `PATH`, then confirm detection with:

```bash
dvbfixer msa --list-engines
```

`--engine auto` selects the first installed engine in this order: MAFFT,
MUSCLE 5, Clustal Omega. See the official download pages for
[MAFFT](https://mafft.cbrc.jp/alignment/software/),
[MUSCLE 5](https://drive5.com/muscle5/manual/install.html), and
[Clustal Omega](http://www.clustal.org/omega/).

## Modeller installation and license

Salilab MODELLER is currently required for `dvbfixer model` when gaps need to be
built and for final `dvbfixer homology` comparative-model construction. It is
also required when `dvbfixer salign --engine modeller` is selected; the default
Biopython SALIGN engine does not need MODELLER. Fast `model` inputs with no
sequence gaps can take the Python-only copy path without invoking it.

MODELLER needs a free academic license key. Register at
<https://salilab.org/modeller/registration.html>, then set the key in
`<env>/lib/modeller-10.8/modlib/modeller/config.py`.

MODELLER remains the default modeling backend. Experimental runners below are
separate installations, not replacements installed by the main environment.

## Experimental diffusion runner environments

`dvbfixer model --backend diffusion` is an explicit, narrow experimental path;
MODELLER remains the default. The core DVBFixer environment intentionally does
not install PyTorch, diffusion-specific CUDA/MPS stacks, Protenix, Protpardelle,
or checkpoints. OpenMM may still bring its own scientific GPU dependencies.
Create the selected runner environment separately and register it once with its
deployment README. The registration verifies the local checkpoint and installs
a private launcher; normal `model --backend diffusion` calls do not accept or
require runner, profile, checkpoint, or digest paths. DVBFixer never downloads
checkpoints and never falls back to another profile or MODELLER.

- Linux/NVIDIA `protenix-v1-cuda`: follow
  [`deploy/protenix-v1/README.md`](../deploy/protenix-v1/README.md). The portable
  wrapper and one public-CLI reconstruction passed on the pinned NVIDIA host;
  the frozen 231-case production-wrapper replay remains pending.
- Protpardelle now uses multichain **cc94_epoch3100**, with
  `protpardelle-1c-cc94-cuda-v1` on Linux/NVIDIA and
  `protpardelle-1c-cc94-mps-v1` on native Apple Silicon. Follow the
  [cc94 deployment guide](../deploy/protpardelle-1c/cc94.md) for exact weights,
  hashes, environments and registration. Both use 500 ODE updates with native
  relative chain encoding, fixed sequence and source-only local context.
  Linux refinement uses OpenMM Reference; Apple uses a separate CPU interpreter
  with fallback disabled. cc94 has case-specific CUDA evidence on 9CQD; native
  MPS acceptance remains pending. Old cc89 launchers are explicit-replay-only;
  reinstall with the cc94 checkpoint to use new Protpardelle requests.

The supported input slice and directory-bundle output are documented under
[`dvbfixer model`](commands/model.md#experimental-diffusion-scope).
Before inference, inspect static profile status or run the no-inference
environment handshake with `dvbfixer doctor`; see the
[`doctor` command](commands/doctor.md). A successful handshake verifies the
selected wrapper and checkpoint environment but does not replace hardware
cohort acceptance or prove broader scientific support.

## Optional Probe diagnostics

`dvbfixer diagnose` can use MolProbity Probe **2.18.211005**, discovered on `PATH`
or in the active Python environment's `bin`. It remains optional; unavailable or
unsupported versions use Python diagnostics with explicit engine/fallback
provenance. `dvbfixer doctor` reports executable availability. See
[Probe behavior and radius conventions](commands/diagnose.md#optional-probe-engine).

For a reproducible optional build, use the
[audited upstream source](https://github.com/rlabduke/probe/tree/091d117f5bceef87389002dccfe26df78516fdda),
run `make` with a C compiler available, and place the resulting `probe` executable
on `PATH`. Confirm its banner with `probe -version`.

xTB is no longer installed by `environment.yml` or exposed as a preparation
refinement backend. Existing locally installed binaries are left untouched.
The evidence and rationale are in [ADR 0012](adr/0012-remove-xtb-refinement-and-repair-probe.md).

## Optional RESP charge engines

**PySCF** is already included in `environment.yml`:

PySCF is pip-installed automatically by `environment.yml`. It enables the recommended free RESP path:

```bash
dvbfixer parametrize molecule.pdb -c resp --qm-engine pyscf -v
```

The default AM1-BCC charge method needs no QM engine. Select RESP when your
parameterization protocol requires it; see [parametrize](commands/parametrize.md).

**Optional: PSI4/psiresp alternative.** Keep this in a separate environment;
the validated recipe uses its own Python 3.9 and numerical stack, not the
project's supported Python range. PySCF remains included in Standard/Desktop.

```bash
python3 scripts/install_dvbfixer.py --runtime psi4 --manager micromamba \
  --prefix /path/to/new-psi4-runtime --non-interactive
dvbfixer parametrize molecule.sdf -c resp --qm-engine psi4 \
  --psi4-env /path/to/new-psi4-runtime/environment
# Manual alternative:
micromamba create -n psi4 -f environments/psi4.yml
```

Do not install DVBFixer into that QM environment. The main CLI invokes the worker
through the chosen conda-family manager; `--psi4-env` accepts either an existing
name or an explicit prefix. The generated recipe pins compatible Pydantic and
Libint versions: an import-only solve was insufficient, and unbounded Libint
segfaulted in the real native kernel. Linux acceptance exercised HF/6-31G* and
the public psiresp two-stage API on water and complete ethanol parameterization.

The separate recipe also passed native macOS arm64 HF, two-stage RESP water
integration and complete ethanol export through micromamba 2.9.0. PySCF remains
the included route; PSI4 is still an explicitly selected optional runtime.

## NumPy/SciPy BLAS backend

The main conda environment selects OpenBLAS for all three BLAS/LAPACK wrappers.
If an older Linux environment aborts with missing `libmkl_gf_lp64.so.3`, repair it:

```bash
micromamba install -n dvbfixer -c conda-forge \
  'libblas=*=*openblas' 'libcblas=*=*openblas' 'liblapack=*=*openblas'
micromamba run -n dvbfixer python scripts/check_scientific_linalg.py
```

Run the smoke check from the checkout. It exercises real NumPy/SciPy solves.
This policy does not change pip wheels or isolated RESP/research environments.

## Other installation problems

- Python 3.14 breaks the supported PROPKA version: retain the Python upper bound.
- NumPy remains pinned below 2.5 for the Python 3.11 development baseline.
- On macOS Docker bind mounts, create the environment on the container's native
  filesystem if ncurses symlink creation fails.

See [environment troubleshooting](known-issues.md#environment-and-developer-notes)
for symptoms and repair details.

## Next steps

- Try the [quick start](quick-start.md) or [GUI guide](gui.md)
- Browse the [command index](commands/index.md)
- Read the [pipelines](pipelines.md) page for end-to-end recipes
- Skim the [known issues](known-issues.md)
