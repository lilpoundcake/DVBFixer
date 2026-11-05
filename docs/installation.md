# Installation

[User guide](README.md) · [Quick start](quick-start.md) · [Troubleshooting](known-issues.md)

## 1. Get the source

DVBFixer's primary delivery is **source plus local installation**, using the
Standard/Desktop helper below. Manual/editable setup is a separate alternative
for development, not an additional step after the helper. Desktop builds
the GUI locally; it is not a pre-populated scientific environment or container.
Dependencies are obtained from their sources and retain their own terms; this
does not bypass restricted-tool licenses or copyleft obligations for combined
redistribution. Obtain an appropriate MODELLER license before using its workflows.
See the [licensing and distribution policy](licensing-status.md).

Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
and Git first, plus a Python **3.11+** bootstrap interpreter. Mamba and conda
are also supported managers subject to the platform limitations below.
The commands assume a terminal on macOS or Linux; on Windows use WSL.

```bash
git clone https://github.com/lilpoundcake/DVBFixer.git
cd DVBFixer
```

<a id="2-create-the-environment"></a>

## Standard/Desktop source-checkout helper

The helper and the two main profiles have Linux clean-install evidence with micromamba 2.9.0, mamba
2.10.0 and conda 26.9.4. See the [dependency matrix](installation-dependency-matrix.md)
and [acceptance record](installation-profile-acceptance.md) for implemented scope,
resolved builds and remaining acceptance gates. The manual `environment.yml`
input is the Standard alias. `environments/legacy-full.yml` retains the
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

### 2. Install one profile

Run `python3 scripts/install_dvbfixer.py` from the repository root. Select
**New Standard/Desktop installation** from the numbered menu, then select the
profile and manager by number and supply an **environment name**.
The suggested name includes the profile and release; a free numeric suffix is
offered if it is occupied. The helper uses the first writable environment
directory reported by the manager. It never asks you to type an installation path.
To automate the same name-based workflow, choose **one**:

```bash
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --name dvbfixer-standard --non-interactive --skip-modeller
# Or Desktop:
python3 scripts/install_dvbfixer.py --profile desktop --manager micromamba \
  --name dvbfixer-desktop --desktop-shortcut --non-interactive --skip-modeller
```

These unattended examples defer MODELLER licensing and can return exit **2**;
use interactive hidden entry or the protected descriptor below to license it.
For advanced/custom locations, explicit prefixes remain supported. The activation
examples later on this page use these explicit locations:

```bash
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --prefix "$HOME/micromamba/envs/dvbfixer-standard"

# Or Desktop, with an optional application shortcut:
python3 scripts/install_dvbfixer.py --profile desktop --manager micromamba \
  --prefix "$HOME/micromamba/envs/dvbfixer-desktop" --desktop-shortcut
```

The no-argument helper first offers a numbered action menu: **new installation**,
**update**, **add profile/GUI to an existing environment**, or **additional
backend**. Profile, manager, backend, optional add-on and confirmation prompts are
also numbered; typing implementation identifiers is unnecessary. Names and
operator-supplied asset/license inputs remain ordinary text or protected entry.
Backend choices are separate and do not install the main application.
The helper installs a non-editable
package and all selected dependencies; Desktop also builds its own GUI assets.
Do **not** follow this with `pip install -e .`, `npm ci` or manual environment
creation for normal use. The bootstrap Python is separate from the managed
scientific environment, whose supported Python range stays `>=3.11,<3.14`.

Choose a **new prefix without whitespace**: MAFFT's upstream shell scripts
failed the real spaced-prefix smoke. Source/input paths may contain spaces.
`--desktop-shortcut` explicitly opts into a per-user macOS application or Linux
application-menu entry. It is optional; interactive installation asks before
creating it. Each shortcut is bound to its installation and does not overwrite
another application's files. The launcher opens the browser only after server
readiness. Native shortcut behavior still requires platform-specific desktop
acceptance; the command-line launcher is tested separately.

The helper refuses existing environments; `--resume` is limited to installations
it previously recorded and continues an interrupted installation. It reads the
existing profile/backend from that directory's record rather than asking for a
new target. Use `--update` for routine Standard/Desktop updates instead.
`--dry-run` does not change the selected environment or publish ownership. Profile
extension also runs native/pip solver previews and retains preview logs; solver
caches and temporary build metadata may be written outside the target environment.
Names are not paths: slash-separated paths, whitespace and reserved base/root
names are rejected. `--name` and `--prefix` are mutually exclusive. Only the
advanced `--prefix` option accepts a filesystem path; relative explicit paths are
resolved against the current working directory and printed, never interpreted as
environment names.

Every long-running command shows a named stage, live output and an elapsed-time
notice every ten seconds while it runs. Complete command output is retained in
the printed `<prefix>.installation-logs/<run-id>/` directory beside the
installation. Protected MODELLER configuration/import output is never logged.
The final status distinguishes ready, action-required and failure; `--json` also
prints the technical plan/result when needed. Cancelling retains diagnostics.

<a id="3-configure-modeller"></a>

### 3. Configure MODELLER

Gap rebuilding and comparative modeling require an appropriate
[MODELLER license](#modeller-installation-and-license). The interactive helper
offers hidden key entry. Empty input defers licensing; it does not complete a
licensed installation. Other tools can be checked without it.

For automation, supply the license through an open protected file descriptor,
not a command-line key. To configure a deferred license without reinstalling:

```bash
# KEY_FILE is a user-owned 0600 file containing only the key, not its value in argv.
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --prefix "$HOME/micromamba/envs/dvbfixer-standard" --non-interactive --resume --configure-only \
  --modeller-key-fd 3 3< "$KEY_FILE"
```

Use the same profile, prefix and add-ons as the original installation (for
Desktop, replace `standard` with `desktop` in both profile and prefix).
For a new unattended install, omit `--resume --configure-only`. Explicitly
deferring licensing with `--non-interactive --skip-modeller` can return **2**;
do not treat that as a complete licensed setup.

<a id="4-verify-the-installation"></a>

### 4. Verify and use the installation

The helper runs real functional smokes and records their commands, results and
versions below `<prefix>/share/dvbfixer/installation-smoke/`. Finish with
`dvbfixer doctor --profile standard` (or `desktop`). Exit codes are **0** for
passed installation prerequisites, **2** for required manual action such as a
deferred MODELLER license, and **1** for failure. Scientific/hardware acceptance
remains scoped separately in the acceptance record. Activate the installed
prefix in each new terminal; choose the matching block:

```bash
micromamba activate "$HOME/micromamba/envs/dvbfixer-standard"
dvbfixer --version
# dvbfixer 0.10.3
dvbfixer doctor --profile standard
dvbfixer msa --list-engines
```

For Desktop instead:

```bash
micromamba activate "$HOME/micromamba/envs/dvbfixer-desktop"
dvbfixer doctor --profile desktop
dvbfixer gui
```

Alternatively open the helper-created DVBFixer application shortcut, without
shell activation. Run `dvbfixer gui --stop` to stop the server and its managed
jobs; Ctrl+C also stops a terminal-launched session. Repeating `dvbfixer gui`
opens the same running installation/data directory rather than a second server.
See [GUI launch options and troubleshooting](commands/gui.md).

Without activation, use `micromamba run -p` with the same prefix, for example:

```bash
micromamba run -p "$HOME/micromamba/envs/dvbfixer-standard" dvbfixer doctor --profile standard
```

Continue to the [quick start](quick-start.md) or [GUI guide](gui.md).

Doctor also checks PDBFixer's `soft.xml` and canonical template resources in an
isolated subprocess. An importable package can still be unusable if those files
are missing/corrupt; `doctor --profile` reports `pdbfixer-resources` as an error.
Preserve the affected prefix and create a fresh complete installation rather
than copying individual assets into it. This resource probe is not a full conda
ownership/integrity audit or a substitute for scientific acceptance.

### AmberTools auxiliary metadata compatibility

The normal Standard/Desktop scientific paths work with the recorded NumPy 2.x
environments, but AmberTools auxiliary distributions (`proprep`, `ndfes`,
`fetkutils`, `edgembar`, `packmol-memgen`) declare narrower Python dependencies.
For a **new** environment that also satisfies those metadata requirements, use:

```bash
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --addon ambertools-compat --prefix /path/to/new-compatible-standard \
  --non-interactive --modeller-key-fd 3 3< "$KEY_FILE"
```

The explicit add-on intersects the normal Python bounds with `<3.13`, selects
NumPy `>=1.26,<2` and Biopython `>=1.83,<1.86`, and installs `pdb2pqr`. These
constraints come from the auxiliary distributions, not a disabled dependency
check. It keeps **every Standard scientific capability, including PySCF/RESP**;
Desktop may select the same add-on. Do not downgrade an existing environment in
place. The default composition and separate diffusion/PSI4 runtimes are unchanged.

The native Apple candidate and a real helper-created Standard passed `pip check`
and 28 extended licensed scientific checks. A helper-created Desktop also passed
`pip check`, 14 licensed science checks, GUI407/typecheck/build and Chrome.
This is not yet a new Linux manager-matrix acceptance or a fix for
native conda's unrelated gettext/Tk path collisions. The catalog-generated
manual input is `environments/standard-ambertools-compat.yml`; normal `environment.yml`
remains the Standard alias. Verify the selected environment with
`dvbfixer doctor --profile standard --addon ambertools-compat` and `python -m pip check`.

### Native Apple Probe alternative

The managed `--addon probe` recipe currently has no bioconda osx-arm64 package.
The exact supported Probe can instead be built natively from its pinned upstream
source using Apple's Command Line Tools. Use a **new** source directory:

```bash
git clone https://github.com/rlabduke/probe.git /path/to/new-probe-source
git -C /path/to/new-probe-source checkout --detach 091d117f5bceef87389002dccfe26df78516fdda
make -C /path/to/new-probe-source probe CFLAGS=-O2
/path/to/new-probe-source/probe -version
# Expected: probe.2.18.211005
PATH="/path/to/new-probe-source:$PATH" micromamba run -p /path/to/new-standard \
  dvbfixer doctor --profile standard --addon probe
```

Keep the source directory/binary and put it on PATH for scientific invocations
that should use Probe. This manual route passed native arm64 build, exact version,
profile doctor and real electron/nuclear-radius diagnostics. Without the binary
on PATH, diagnose retains its explicit Python fallback. Neither a different Probe
version nor Rosetta is substituted. This does **not** mean that the managed
`--addon probe` installation can solve on Apple; omit that helper add-on when
using this separately built executable. See the
[acceptance record](installation-profile-acceptance.md) for commands/results.

## Update an existing installation

Close the selected installation's GUI and running jobs. From the source checkout:

```bash
git pull
python3 scripts/install_dvbfixer.py --update
```

The interactive helper queries the selected manager and displays a **numbered
list of existing installations**, labelled by profile or manual/package-only
status, installed version and environment name. Choose a number, not a path.
Managed targets are verified against their ownership record; older/manual targets
are detected from package metadata or an existing CLI entry point. The helper
shows installed and checkout versions and asks for confirmation
before changing packages. Multiple versions may coexist; select only the one you
intend to update. An update uses the current checkout, not an automatic download
of a release. To inspect without changing packages, add `--dry-run`.

For an older/manual installation the helper first asks **what to update**:
**1 — only DVBFixer**, **2 — the full installation**, or **0 — cancel**. Choosing
the second item leads to Standard/Desktop profile selection and its dependency
preview; it is not a disguised cancellation. Package-only application has its own
separate final confirmation. A full update covers the selected application
profile (and GUI for Desktop), not arbitrary unrelated packages or automatic
replacement of frozen model checkpoints.

Discovery reads manager `info`/`env list` JSON plus ownership records directly in
the reported environment directories, current directory and source checkout.
For composite backends it also checks parents of registered sampling/refinement
environments. It does not recursively scan the filesystem or execute discovered
environment Pythons just to build the list. Wheel metadata, legacy `egg-info` and
editable `egg-link` layouts are supported. A retained `bin/dvbfixer` entry point
also identifies metadata-incomplete installations, labelled with an unknown
version rather than silently disappearing. Unusual locations not found this way remain
accessible through advanced `--prefix`.

For automation, select an unambiguous environment name and manager:

```bash
python3 scripts/install_dvbfixer.py --update --manager micromamba \
  --name dvbfixer-desktop --non-interactive
```

If several installations share that name, use the interactive numbered list or
an explicit advanced prefix:

```bash
python3 scripts/install_dvbfixer.py --update --manager micromamba \
  --prefix "$HOME/micromamba/envs/dvbfixer-desktop" --non-interactive
```

For **helper-managed installations**, the existing Standard/Desktop profile and add-ons are preserved. Dependencies and
the DVBFixer package are reconciled, Desktop assets are rebuilt, and the same
functional/Doctor checks run. Existing MODELLER configuration is read privately
before package updates and restored afterward; it is not printed or logged.
Workspace data and unrelated environments are not update targets. Existing
shortcuts remain bound to the same prefix; no new shortcut is required.

Updates are **in-place, not transactional**: a failed package operation may leave
the environment partially updated, and no automatic rollback is promised. Keep
your license and important project backups. Retain the printed logs and use
`--resume` with the same directory to recover an interrupted operation. Exit 2
still means action required, not a fully ready environment.

### Older and manually maintained installations

These are selectable in the same numbered update list, explicitly labelled
**existing/manual; package-only update**. After confirmation the helper reinstalls
only DVBFixer from this checkout using `pip install --upgrade --no-deps`, then runs
Doctor. It does not solve native dependencies, modify MODELLER configuration,
create an ownership record, or silently convert the environment to Standard/Desktop.
Known editable installs stay editable, bound to the checkout shown before
confirmation. Neither discovery nor a dry run changes the environment.

For unattended updates, explicit package-only consent is required:

```bash
python3 scripts/install_dvbfixer.py --update --package-only --manager micromamba \
  --name dvbfixer --non-interactive
```

The current checkout still requires a compatible Python/scientific stack; missing
or incompatible dependencies need separate attention. A package update is not a
claim of profile readiness. Manual GUI assets are not rebuilt by this operation:
source Desktop users can run `npm ci` and `npm run build` in `gui/` separately.
If interrupted, rerun the same package-only update, not `--resume`, which remains
reserved for recorded installations. Frozen runtime interpreters are excluded
from manual application updates even if a DVBFixer package is present in them.

Frozen diffusion/QM runtimes are separate: use their
[backend resume instructions](installation-runtime-provisioning.md#evidence-and-resume),
not `--update`, which is reserved for the main application.

## Add a profile or GUI to an existing environment

```bash
python3 scripts/install_dvbfixer.py --extend
```

Select the manager and an existing environment from numbered lists, then choose:

1. **Standard**: complete the scientific CLI profile, including PySCF/RESP.
2. **Desktop**: Standard plus Node and the prefix-local built GUI.
3. **GUI only**: add/update Node, the DVBFixer package and built GUI while retaining
   the existing scientific stack. Requires an existing DVBFixer installation.
4. **Diffusion backend**: choose a platform-appropriate maintained model runtime,
   then create a new isolated runtime by name or continue an existing recorded
   runtime selected by number. Its manager can differ from the scientific
   environment's manager. The selected application environment is not modified.

The diffusion transition reuses the existing
[runtime provisioning workflow](installation-runtime-provisioning.md), frozen
source/asset checks and preflight-before-publication rules. It never installs
PyTorch into Standard/Desktop. Launchers use the existing **per-user** discovery,
so they may be visible to other compatible installations; this is not an exclusive
binding to the selected application prefix. Weights are still operator-supplied,
and an empty checkpoint leaves action required rather than a ready backend.
Provisioning does not upgrade the selected application or establish scientific
qualification. PSI4 remains available through the main additional-backend menu.

Existing DVBFixer installations are shown first. Environments without DVBFixer
are available through **Other existing environments…**, rather than filling the
default list. Manager bases (including other discoverable conda installations)
and macOS application bundles such as PyMOL are excluded.

Full profiles can be attached to an existing named conda-family environment even
if it predates the helper or does not yet contain DVBFixer. The manager base,
frozen backend roots/interpreters and invalid ownership records are not eligible.
Existing Desktop ownership and recorded add-ons are preserved when attaching a
full profile. GUI-only does not convert a manual/Standard environment to a
qualified Desktop scientific profile.

Unattended equivalents use names instead of installation paths:

```bash
python3 scripts/install_dvbfixer.py --extend --profile desktop \
  --manager micromamba --name dvbfixer --non-interactive --skip-modeller

# Or add just the GUI without changing the scientific stack:
python3 scripts/install_dvbfixer.py --extend --gui-only \
  --manager micromamba --name dvbfixer --non-interactive
```

`--gui-only` and `--profile` are separate scopes; do not combine them. Explicit
`--prefix` remains available for advanced custom locations. `--desktop-shortcut`
works for Desktop or GUI-only extension. Interactive installation offers a
numbered Yes/No choice; no manual npm setup is required.

### Preview, apply and failure behavior

The helper requests a channel-isolated **native dry-run transaction** showing
packages to add, replace/downgrade or remove, then a **pip dry-run report**. It
displays both before interactive confirmation. Dependency conflicts stop before
package installation; no `--clobber` or relaxed safety policy is used. The
existing Python minor version is held fixed so pip can be previewed against the
same ABI. If a full profile needs a different Python minor, create a new profile
environment. Existing conda Open Babel is incompatible with the full profile's
wheel-only route; the helper refuses to overlay the two. For an interactive full
profile request in such an existing DVBFixer environment, it offers a numbered
choice to **add GUI only** or **stop without changes** instead of exiting with a
dead-end error. GUI-only requires an explicit scope change, gets a fresh preview
and confirmation, and skips newly requested scientific add-ons/license changes.
The full profile still requires a compatible separate environment; automated
requests never silently fall back to a different scope. To request only GUI
directly, use `--extend --gui-only` instead of `--profile desktop`.
PyTorch/diffusion stacks
remain separate from full Standard/Desktop profiles. GUI-only can retain that
scientific stack, but refuses any native plan that changes existing packages
other than Node itself.

After confirmation the helper checks that the environment and native plan have
not changed, applies native dependencies, verifies their identities and rechecks
pip. Newly proposed pip packages/versions outside the shown plan are refused.
Native and pip are **separate transactions**: a later pip/build/check failure can
leave native changes applied. No automatic rollback is promised. Review the
retained logs and repeat `--extend` to preview a recovery plan; do not use
`--resume` to publish an uncompleted profile attachment.

The existing MODELLER license value is preserved privately for full-profile changes;
GUI-only neither rewrites nor requires a new license. Workspace files and `.env`
settings survive GUI rebuilding. Known editable DVBFixer installations retain
editable mode for GUI-only updates, bound to the checkout shown in the preview.
Full-profile ownership is published only after
functional, Doctor and requested GUI checks, not at the start. Exit 2 remains
action-required (for example a deferred license), not ready. GUI-only runs its
launcher/diagnostic checks without claiming Standard/Desktop science readiness.
Failed or cancelled previews retain logs under the checkout's
`tmp/installer-profile-preview-*/`; applied operations retain their printed
per-installation logs. These are workflow safeguards, not new platform/scientific
qualification evidence.

## Manual installation and source development

This is an **alternative** to the helper, for contributors or operators managing
their own environment. It does not perform the helper's licensing, smoke evidence
or application-shortcut setup. Run from the repository root:

```bash
micromamba create -f environments/standard.yml
micromamba activate dvbfixer-standard
python -m pip install -e .
```

`environment.yml` is the generated Standard alias (with environment name
`dvbfixer` instead of `dvbfixer-standard`). For GUI development, use
`environments/desktop.yml` and activate `dvbfixer-desktop` instead. Keep the
checkout in place for editable installs and configure the
[MODELLER license](#modeller-installation-and-license) manually. Then check
`dvbfixer doctor --profile standard` (or `desktop`). Profile specifications can
be checked with `PYTHONPATH=src python scripts/gen_installation_profiles.py --check`.

Manual Desktop builds run in the checkout's `gui/`: follow
[source development](../gui/README.md#source-development) for `npm ci` and
`npm run dev`, or `npm run build` followed by `npm start` for a built server.
Doctor recognizes an editable source build, but **`dvbfixer gui` requires the
helper-installed assets inside its prefix**; it does not launch that checkout.

Both profiles use **openbabel-wheel**. Do not overlay conda `openbabel`; ACPYPE's
Python metadata requires the wheel. The `obminimize` entrypoint launches the
wheel's native tool with its plugin paths initialized. A wheel alone does not
provision the full scientific environment or Desktop's Node/GUI assets.

For an already helper-built Desktop, the advanced diagnostic fallback is
`npm start` in `<prefix>/share/dvbfixer/gui/` with that prefix active; normal
users should use `dvbfixer gui` instead.

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
Babel. Both are installed by Standard and Desktop, and pip installations declare
`gemmi` plus the prebuilt `openbabel-wheel` explicitly. Verify the active
environment with `python -c "import gemmi; from openbabel import openbabel"`.
If only PDB input is used, the conversion layer is not invoked.

## AmberTools and MolProbity Reduce

The `tleap-reduce` preparation backend requires both `tleap` (AmberTools) and
the `reduce` executable (MolProbity). Standard and Desktop declare `ambertools`
and `reduce` separately; an existing environment containing `tleap` is not proof
that Reduce is installed. Repair an incomplete environment with:

```bash
# With the intended Standard or Desktop environment active:
micromamba install -p "$CONDA_PREFIX" --override-channels -c conda-forge -c bioconda ambertools reduce
dvbfixer doctor --profile standard  # Use desktop for Desktop.
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
# In the installed Standard or Desktop environment:
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
