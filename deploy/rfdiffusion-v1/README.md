# RFdiffusion v1 benchmark environment

This environment is isolated from the DVBFixer core environment. It is pinned
to RFdiffusion v1 commit
`bf42b54c20a99dd7350456c85985ed4d83b95d48` and the independently measured
`Base_ckpt.pt` SHA-256
`0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca`.

Keep reusable external artifacts under the project-local ignored
`.artifacts/` directory by default. Set `DVBFIXER_ARTIFACT_ROOT` when a larger
or shared local filesystem is preferable; the directory is operational cache,
not reviewed source and must never be committed:

```bash
ARTIFACT_ROOT="${DVBFIXER_ARTIFACT_ROOT:-$PWD/.artifacts}"
RF_ROOT="$ARTIFACT_ROOT/rfdiffusion-v1"
RF_ENV="$RF_ROOT/env"
RF_SOURCE="$RF_ROOT/source"
RF_CHECKPOINT="$RF_ROOT/Base_ckpt.pt"
mkdir -p "$RF_ROOT"
```

Create the environment only when `$RF_ENV` is absent, clone the pinned source
only when `$RF_SOURCE` is absent, and install its vendored SE(3)-Transformer
package:

```bash
if [ ! -x "$RF_ENV/bin/python" ]; then
  micromamba create -f deploy/rfdiffusion-v1/environment.yml -p "$RF_ENV"
fi
if [ ! -d "$RF_SOURCE/.git" ]; then
  git clone https://github.com/RosettaCommons/RFdiffusion.git "$RF_SOURCE"
fi
git -C "$RF_SOURCE" checkout bf42b54c20a99dd7350456c85985ed4d83b95d48
if ! "$RF_ENV/bin/python" -c 'import se3_transformer' >/dev/null 2>&1; then
  micromamba run -p "$RF_ENV" python -m pip install --no-build-isolation \
    "$RF_SOURCE/env/SE3Transformer"
fi
```

Download the checkpoint only from the upstream URL and verify it before use:

```bash
if [ ! -f "$RF_CHECKPOINT" ]; then
  curl -fL \
    https://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt \
    -o "$RF_CHECKPOINT.partial"
  mv "$RF_CHECKPOINT.partial" "$RF_CHECKPOINT"
fi
printf '%s  %s\n' \
  0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca \
  "$RF_CHECKPOINT" | sha256sum --check
```

The checksum is independently measured because upstream does not publish a
signed checksum manifest. The upstream license was amended after v1.0.0 to say
that it covers the model weights linked from the README; redistribution still
requires release review before a checkpoint is embedded in a container.

The internal adapter is `python -m dvbfixer.model.diffusion.rfdiffusion_v1`.
It is not a public modeling backend. A future service image must preserve the
same request/result protocol, use an immutable base-image digest, and fetch or
mount the verified checkpoint according to the final redistribution decision.

## Minimal container runner

`Dockerfile` reproduces the tested two-environment boundary: Python 3.11.16
with OpenMM runs the DVBFixer adapter, while pinned Python 3.9.19/CUDA 11.1
runs RFdiffusion. The adapter environment intentionally contains only the
diffusion dependencies and does not install or invoke PROPKA; the main DVBFixer
environment retains its documented `>=3.11,<3.14` constraint. The image's
linux/amd64 micromamba base is pinned by digest, RFdiffusion is checked out at
the audited commit during the build, and the checkpoint is not copied into the
image.

Build from the repository root on a host with Docker:

```bash
docker build \
  --file deploy/rfdiffusion-v1/Dockerfile \
  --tag dvbfixer/rfdiffusion-v1:bf42b54 .
docker image inspect dvbfixer/rfdiffusion-v1:bf42b54 \
  --format '{{.Id}}'
```

Run one already-built benchmark workspace with no network access and a
read-only checkpoint mount:

```bash
WORKSPACE="$PWD/.artifacts/diffusion-benchmarks/7k8s-insertion-seed7-run1"
docker run --rm --gpus device=0 --network none \
  --user "$(id -u):$(id -g)" \
  --volume "$WORKSPACE:/workspace" \
  --volume "$RF_CHECKPOINT:/models/Base_ckpt.pt:ro" \
  dvbfixer/rfdiffusion-v1:bf42b54
```

The entrypoint refuses a missing or digest-mismatched checkpoint before model
loading. Record the resulting image ID/digest and repeatability output before
using the image as deployment evidence.

The following are deliberately deferred until an actual deployment requires
them: registry publication/signing, multi-architecture builds, embedded
checkpoint redistribution, Kubernetes manifests, and remote scheduling. They
do not solve a current benchmark-runner problem and are not part of Phase 2's
scientific acceptance claim.

The adapter enables DVBFixer's seeded post-sampling OpenMM boundary refinement
by default. Set the internal backend option `boundary_refinement=false` only
for the declared raw-backbone ablation. This pass freezes all source atoms,
constrains generated heavy-atom bonds, and independently rechecks the complete
candidate; it is not per-denoising-step reinjection. The acceptance path uses
OpenMM's deterministic `Reference` platform and seeded hydrogen placement.
Faster CPU/CUDA refinement is not accepted as a repeatability substitute.

Frozen withheld-coordinate workspaces and their evidence are reproducible with:

```bash
python scripts/build_diffusion_benchmark_request.py CASE_ID WORKSPACE
# Run the internal adapter twice in separate workspaces with the pinned paths above.
python scripts/analyze_diffusion_benchmark_pair.py FIRST_WORKSPACE SECOND_WORKSPACE
python scripts/analyze_modeller_benchmark.py MODELLER_WORKSPACE MODEL_1 MODEL_2
```

The MODELLER analyzer validates the complete target sequence before mapping
its output residues by sequence ordinal back to request identities. This is
needed when MODELLER sequentially renumbers a withheld insertion-code region;
the normalized PDB exists only in a temporary contained directory, and source
MODELLER coordinates and artifacts are never modified.

Reusable generated workspaces may live under
`$ARTIFACT_ROOT/diffusion-benchmarks/`. They are ignored because they contain
large backend outputs and machine-local paths; accepted measurements belong in
the reviewed research inventory, not as committed runtime directories.

The builders use the reviewed coordinate-fragment sequence as the sampler
target so unrelated natural terminal/internal gaps in a deposited FASTA do not
silently turn a one-gap benchmark into a different request. The inventory still
records and validates the corresponding full-FASTA interval.

Cases whose source fixture contains unrelated chains or heterogens may declare
`request_structure_scope = "target-protein-chain-only"`. The builder then emits
both masked input and native reference from only that canonical ATOM chain. Such
cases measure isolated loop reconstruction only and must not be cited as
interface- or ligand-context evidence.

Interface-adjacent cases instead declare
`request_structure_scope = "target-and-context-protein-chains"` plus explicit,
case-sensitive `context_chains`. The adapter maps the target to synthetic chain
`A`, maps partner chains deterministically to `B..Z`, and emits RFdiffusion's
documented receptor-contig form (`target/0 B1-N`). The pinned base checkpoint is
intentional: upstream documents separate-chain motif inpainting with this form;
its complex checkpoint is selected for hotspot-driven binder design, which this
adapter does not perform. Partner backbone atoms participate in the shared
post-sampling Kabsch fit and must remain within the frozen `0.5 Angstrom` RMSD
and `2.0 Angstrom` maximum-displacement limits before source records are
reinserted exactly. The `7x35-chain-a-interface-5` case has completed pinned A100
same-seed repeatability, independent validation, and MODELLER comparison; exact
measurements and the recorded partner-threshold revision are in the research
inventory.
