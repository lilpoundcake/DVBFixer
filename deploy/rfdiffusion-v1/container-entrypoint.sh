#!/bin/sh
set -eu

EXPECTED_CHECKPOINT_SHA256="0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca"
CHECKPOINT="${DVBFIXER_RF_CHECKPOINT:-/models/Base_ckpt.pt}"

if [ ! -f "$CHECKPOINT" ]; then
    printf '%s\n' "missing mounted RFdiffusion checkpoint: $CHECKPOINT" >&2
    exit 2
fi

ACTUAL_CHECKPOINT_SHA256="$(sha256sum "$CHECKPOINT" | cut -d ' ' -f 1)"
if [ "$ACTUAL_CHECKPOINT_SHA256" != "$EXPECTED_CHECKPOINT_SHA256" ]; then
    printf '%s\n' \
        "RFdiffusion checkpoint digest mismatch: expected $EXPECTED_CHECKPOINT_SHA256, got $ACTUAL_CHECKPOINT_SHA256" >&2
    exit 2
fi

exec /opt/conda/envs/adapter/bin/python \
    -m dvbfixer.model.diffusion.rfdiffusion_v1 \
    --repository /opt/conda/rfdiffusion-source \
    --python /opt/conda/envs/rfdiffusion/bin/python \
    --checkpoint "$CHECKPOINT" \
    "$@"
