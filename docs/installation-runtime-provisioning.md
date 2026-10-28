# Separate diffusion-runtime provisioning

The source-checkout helper provisions the maintained **cc94 CUDA**, **cc94 MPS**
and **Protenix v1 CUDA** routes independently of Standard/Desktop. Composition
and deployment-input references live in `src/dvbfixer/data/installation-profiles.toml`.
It uses the same micromamba/mamba/conda adapters as the main installer.

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
No checkpoint or Protenix auxiliary dataset is downloaded by the helper.

Use `--dry-run` to inspect the route and required manual inputs without creating
files. The interactive helper also offers runtime choices and deferred asset
entry. Unsupported platforms are rejected before creating a runtime root.

## Supply assets and finish registration

cc94 requires the operator-provided `cc94_epoch3100.pth` and its digest-matching
`../configs/cc94.yaml` in the layout described by the
[cc94 deployment guide](../deploy/protpardelle-1c/cc94.md):

```bash
python3 scripts/install_dvbfixer.py --runtime protpardelle-cuda \
  --manager micromamba --prefix /path/to/new-cc94-runtime --non-interactive --resume \
  --checkpoint /path/to/model_params/weights/cc94_epoch3100.pth
```

On native Apple Silicon select `--runtime protpardelle-mps`. This uses only the
Apple portability patch, native PyTorch 2.6.0 with MPS fallback disabled, and a
separate OpenMM 8.5.2 CPU refiner. **Native MPS provisioning/execution has not been
tested by this Linux implementation increment.** cc94 MPS scientific acceptance
remains pending; CUDA or historical cc89 evidence does not transfer to it.

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
different launchers are preserved unless `--replace-runner` is explicit;
identical private launchers are idempotent on resume.

## Evidence and resume

`evidence/installation.json` records per-step commands, return codes, spec/patch
digests, logs, source identity and canonical preflight facts. Package versions,
conda builds and interpreter prefixes are retained in the `versions-*` logs.
Failed or interrupted installations retain diagnostics and source files.

`--resume` accepts only roots recorded by the helper. It verifies source integrity
and the original recipe/manager identity, then reconciles native and pip packages.
It does not reset source edits or clean untracked modules. A changed recipe needs
a new prefix, preserving the old runtime and its evidence.

Current Linux evidence covers both CUDA runtimes through all three managers,
plus checkpoint-backed public `model` execution through the newly provisioned
micromamba runtimes. The first Protenix real sample exposed missing Click/Requests
and a relative refinement-output path bug; these were fixed and the failed
attempts were retained. These are installation/execution checks, not a new
scientific cohort or a change to the accepted profile identities.
