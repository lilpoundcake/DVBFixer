# RFdiffusion v1 benchmark environment

This environment is isolated from the DVBFixer core environment. It is pinned
to RFdiffusion v1 commit
`bf42b54c20a99dd7350456c85985ed4d83b95d48` and the independently measured
`Base_ckpt.pt` SHA-256
`0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca`.

Create the environment, clone the pinned source, and install its vendored
SE(3)-Transformer package:

```bash
micromamba create -f deploy/rfdiffusion-v1/environment.yml -p "$RF_ENV"
git clone https://github.com/RosettaCommons/RFdiffusion.git "$RF_SOURCE"
git -C "$RF_SOURCE" checkout bf42b54c20a99dd7350456c85985ed4d83b95d48
micromamba run -p "$RF_ENV" python -m pip install --no-build-isolation \
  "$RF_SOURCE/env/SE3Transformer"
```

Download the checkpoint only from the upstream URL and verify it before use:

```bash
curl -fL \
  https://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt \
  -o Base_ckpt.pt
printf '%s  %s\n' \
  0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca \
  Base_ckpt.pt | sha256sum --check
```

The checksum is independently measured because upstream does not publish a
signed checksum manifest. The upstream license was amended after v1.0.0 to say
that it covers the model weights linked from the README; redistribution still
requires release review before a checkpoint is embedded in a container.

The internal adapter is `python -m dvbfixer.model.diffusion.rfdiffusion_v1`.
It is not a public modeling backend. A future service image must preserve the
same request/result protocol, use an immutable base-image digest, and fetch or
mount the verified checkpoint according to the final redistribution decision.

The adapter enables DVBFixer's seeded post-sampling OpenMM boundary refinement
by default. Set the internal backend option `boundary_refinement=false` only
for the declared raw-backbone ablation. This pass freezes all source atoms,
constrains generated heavy-atom bonds, and independently rechecks the complete
candidate; it is not per-denoising-step reinjection.
