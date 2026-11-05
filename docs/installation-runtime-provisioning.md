# Separate diffusion and QM runtime provisioning

The source-checkout helper provisions the maintained **cc94 CUDA**, **cc94 MPS**
and **Protenix v1 CUDA** routes independently of Standard/Desktop. Composition
and deployment-input references live in `src/dvbfixer/data/installation-profiles.toml`.
It uses the same micromamba/mamba/conda adapters as the main installer.
The same helper also provisions the isolated optional PSI4/psiresp environment.

Diffusion is also reachable after selecting an existing application: choose
**Add a profile or GUI**, select the installation, then choose **Diffusion
backend**. The submenu shows maintained routes for the host platform and offers
new named runtime or existing recorded runtime continuation. It reuses the already
selected manager; **Change runtime manager** is a separate explicit option.
Runtime selection has its own heading and shows each full path and retained
provisioning state, even when the application and runtime share a basename.
This is a transition
to separate provisioning, not a co-installation into the scientific environment.
Per-user launcher registration can make a validated backend visible to multiple
compatible installations; no exclusive per-prefix binding or application upgrade
is implied. The required checkpoints/assets and scientific acceptance limits below
remain unchanged.

## Prepare a new runtime

From a source checkout, with a Python 3.11+ bootstrap and the selected manager
already available:

```bash
python3 scripts/install_dvbfixer.py --runtime protpardelle-cuda \
  --manager micromamba --prefix /path/to/new-cc94-runtime --non-interactive
```

The prefix is a **new runtime root**, containing `engine/`, `sampling/`,
`refinement/` where needed, and `evidence/`. Use a prefix without whitespace.
The helper clones the pinned engine into this root and applies the maintained
patch. `--engine-source /path/to/existing/git/checkout` can use an existing local
clone as a read-only Git source; the operator's source is never reset, cleaned
or patched. Git LFS smudging is disabled, so source cloning does not fetch weights.

Without weights, provisioning still runs native accelerator and refinement
kernel checks. It returns **action-required** (exit 2), lists the missing assets
and publishes no launcher. Engine imports that require those assets are deferred.
No model download happens without explicit consent. Protenix auxiliary data and
Kalign remain operator-provided; they are not included in the weights download.

Use `--dry-run` to inspect the route and required manual inputs without creating
files. In the interactive helper select **additional model/QM backend** from the
numbered action menu, then select the runtime by number; Standard/Desktop are not
mixed into this list. Deferred asset
entry is allowed. Unsupported platforms are rejected before creating a runtime root.

## Managed weights storage

The interactive installer shows required filenames, official source, transfer
size and model terms before offering three choices: **download into shared
models/**, **use a local weights file**, or **configure later**. A verified cache
is reused without a network request. There is no implicit download during normal
model execution, and dry runs do not download or create the cache.

The shared directory belongs to the source checkout and is excluded from Git:

```text
models/
  protpardelle-cc94/model_params/
    weights/cc94_epoch3100.pth
    configs/cc94.yaml
  protenix-v1/model_params/weights/protenix_base_default_v1.0.0.pt
```

CUDA/MPS cc94 runtimes reuse the same weights/config. The official
[Zenodo archive](https://zenodo.org/records/16817230) is approximately 4.3 GB
(CC BY 4.0); only the required two regular files are retained, not other models
in the archive. Its published archive checksum and both frozen file SHA256
digests are verified. Protenix v1 weights are approximately 1.48 GB from the
official upstream checkpoint host, under its published Apache 2.0 terms, and
must match the existing frozen SHA256. This does not support Protenix v2 or cc89.

Non-interactive opt-in requires both flags after reviewing the displayed terms:

```bash
python3 scripts/install_dvbfixer.py --runtime protpardelle-mps \
  --manager micromamba --name cc94-runtime \
  --download-weights --accept-model-terms --non-interactive
```

Use the corresponding CUDA runtime on Linux. A preview with `--dry-run` prints
the plan without a transfer. Local `--checkpoint` input remains supported and
cannot be combined with `--download-weights`. Interrupted downloads are removed
from private staging; a retry restarts that transfer. Previously verified files
are reused. Corrupt/incomplete existing cache entries fail closed and are not
silently overwritten: preserve and inspect them, or supply a separate verified
local file. Symlink cache roots and required archive links are refused.

## Supply assets and finish registration

cc94 requires `cc94_epoch3100.pth` and its digest-matching
`../configs/cc94.yaml` in the layout described by the
[cc94 deployment guide](../deploy/protpardelle-1c/cc94.md):

```bash
python3 scripts/install_dvbfixer.py --runtime protpardelle-cuda \
  --manager micromamba --prefix /path/to/new-cc94-runtime --non-interactive --resume \
  --checkpoint /path/to/model_params/weights/cc94_epoch3100.pth
```

On native Apple Silicon select `--runtime protpardelle-mps`. This uses only the
Apple portability patch, native PyTorch 2.6.0 with MPS fallback disabled, and a
separate OpenMM 8.5.2 CPU refiner. Native micromamba 2.9.0 provisioning, real
MPS/CPU kernel probes and a public 1UBQ model execution passed on macOS arm64.
The digest-verified bundle had one independently accepted candidate after 500
MPS updates with fallback disabled; this case needed no refinement. This is
installation/execution evidence, not broad cc94 MPS scientific qualification.
Representative multi-seed/context research acceptance remains pending; CUDA or
historical cc89 evidence does not transfer to it. See the
[acceptance record](installation-profile-acceptance.md).

Protenix requires all assets from its existing
[frozen profile lock](../deploy/protenix-v1/profile-lock.json):

```bash
python3 scripts/install_dvbfixer.py --runtime protenix-cuda \
  --manager conda --prefix /path/to/new-protenix-runtime --non-interactive \
  --checkpoint /path/to/protenix_base_default_v1.0.0.pt \
  --kalign /path/to/pinned/kalign \
  --auxiliary-root /path/to/data-root
```

The data root must contain the listed `common/` files. The helper validates
their hashes and copies them into its own engine checkout. A wrong supplied
asset is an error, not a deferred success. The existing Protenix runtime uses
Python 3.13.15, PyTorch 2.13.0/CUDA 12.9 and the frozen A100/cgroup contract.
Its maintained runner samples and refines in **one backend process**, isolated
from Standard; provisioning preserves that process layout. Protpardelle uses
separate sampling/refinement interpreters.

The pinned sources are loaded through the existing runner's source-path route.
The helper does not install an upstream editable engine with conflicting
requirements or install Standard/DVBFixer's chemistry dependency set into an
ML runtime. In particular, Protenix's upstream requirements conflict with the
already accepted `profile-lock.json`; the provisioning input preserves the
locked identities instead. See the [acceptance record](installation-profile-acceptance.md).

Registration calls the existing `scripts/install_diffusion_backend.py` into a
private staging directory. The existing bounded runner preflight must pass
before the launcher is published. A successful install does not imply general
scientific acceptance: it is reported separately as `not-tested` by provisioning.

By default publication uses DVBFixer's usual per-user runner directory. Override
it with `--runner-directory /path/to/private/runners` for a contained installation
or test. Put that directory on PATH for discovery by Standard/Desktop. Existing
different launchers are preserved unless replacement is explicit. After preflight,
interactive setup shows the existing launcher path and asks whether to replace it;
declining keeps the existing launcher and returns action-required. Automation uses
`--replace-runner`, or a separate `--runner-directory`. Replacing a per-user launcher
affects discovery by all compatible installations of that user, not just the
selected application. A launcher changed during confirmation is not overwritten;
identical private launchers are idempotent on resume.

## Evidence and resume

`evidence/installation.json` records per-step commands, return codes, spec/patch
digests, logs, source identity and canonical preflight facts. Package versions,
conda builds and interpreter prefixes are retained in the `versions-*` logs.
Failed or interrupted installations retain diagnostics and source files.
Named stages, live command output, ten-second elapsed-time notices and per-step
log paths are shown while provisioning runs. Only the innermost active stage
emits elapsed-time notices. While a menu or other input is waiting for an answer,
all enclosing progress timers and live-log display are paused; the prompt states
that an answer is required. Waiting time is excluded from reported active time,
and retained log output resumes after the answer without being discarded.
The default final output is a short
readable summary, not the full JSON report; `--json` additionally prints that
report. `action-required` explicitly means dependencies may be installed but the
backend is **not ready**, with the required operator inputs listed. Scientific
acceptance remains separate from provisioning success.

The final result distinguishes success, preview, cancellation, required action
and failure. Missing weights list exact filenames, their official source and
the concrete resume command. A completed `versions-*` snapshot is not called an
uncompleted step or mined for unrelated network-repair suggestions. Publication
failures are labelled **Publish validated backend launcher**, with explicit
replacement/private-directory recovery instead of a raw file-exists traceback.

`--resume` accepts only roots recorded by the helper. It verifies source integrity
and the original recipe/manager identity, then reconciles native and pip packages.
It does not reset source edits or clean untracked modules. A changed recipe needs
a new prefix, preserving the old runtime and its evidence.
Interactive `--resume` first separates **application** from **backend** recovery
with a numbered choice,
then lists matching helper-managed installations for number-based selection.
It infers the backend from the selected record, without asking for an installation
path. New backend installations ask for a name and use the selected manager's
environment directory; `--name` also supports automation. Explicit `--prefix`
remains an advanced option for custom locations. Routine Standard/Desktop updates
use the separate [`--update` workflow](installation.md#update-an-existing-installation).

Current Linux evidence covers both CUDA runtimes through all three managers,
plus checkpoint-backed public `model` execution through the newly provisioned
micromamba runtimes. The first Protenix real sample exposed missing Click/Requests
and a relative refinement-output path bug; these were fixed and the failed
attempts were retained. These are installation/execution checks, not a new
scientific cohort or a change to the accepted profile identities.

## Optional PSI4/psiresp runtime

```bash
python3 scripts/install_dvbfixer.py --runtime psi4 --manager micromamba \
  --prefix /path/to/new-psi4-runtime --non-interactive
dvbfixer parametrize molecule.sdf -c resp --qm-engine psi4 \
  --psi4-env /path/to/new-psi4-runtime/environment
```

This route has no engine clone, checkpoint or launcher. It creates only its own
QM environment, records versions/builds and checks an actual native HF kernel.
The main CLI keeps its normal Python/environment and sends a worker to the
selected QM prefix (existing named environments remain supported).

The verified recipe pins PSI4 1.9.1, psiresp 0.4.2, Pydantic <2 and Libint 2.9.0.
The worker uses public QCSchema/AtomicResult and psiresp RESP APIs, carrying the
input bond graph and computing on its Job-owned orientation. Linux and native
Apple micromamba acceptance included actual water charges and full ethanol
parameterization. PySCF/RESP is still included in both main profiles
and is the preferred ordinary free route.
