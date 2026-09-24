FROM mambaorg/micromamba@sha256:008e06cd8432eb558faa4738a092f30b38dd8db3137a5dd3fca57374a790825b

LABEL org.opencontainers.image.title="DVBFixer RFdiffusion v1 benchmark runner" \
      org.opencontainers.image.source="https://github.com/lilpoundcake/DVBFixer" \
      org.opencontainers.image.revision="bf42b54c20a99dd7350456c85985ed4d83b95d48"

COPY --chown=$MAMBA_USER:$MAMBA_USER \
    deploy/rfdiffusion-v1/environment.yml /tmp/rfdiffusion-environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    deploy/rfdiffusion-v1/adapter-environment.yml /tmp/adapter-environment.yml

RUN micromamba create --yes --prefix /opt/conda/envs/rfdiffusion \
        --file /tmp/rfdiffusion-environment.yml && \
    micromamba create --yes --prefix /opt/conda/envs/adapter \
        --file /tmp/adapter-environment.yml && \
    micromamba clean --all --yes

RUN /opt/conda/envs/adapter/bin/git clone \
        https://github.com/RosettaCommons/RFdiffusion.git \
        /opt/conda/rfdiffusion-source && \
    /opt/conda/envs/adapter/bin/git -C /opt/conda/rfdiffusion-source \
        checkout --detach bf42b54c20a99dd7350456c85985ed4d83b95d48 && \
    test "$(/opt/conda/envs/adapter/bin/git -C /opt/conda/rfdiffusion-source rev-parse HEAD)" = \
        "bf42b54c20a99dd7350456c85985ed4d83b95d48" && \
    /opt/conda/envs/rfdiffusion/bin/python -m pip install --no-build-isolation \
        /opt/conda/rfdiffusion-source/env/SE3Transformer && \
    rm -rf /opt/conda/rfdiffusion-source/.git

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml /opt/dvbfixer/pyproject.toml
COPY --chown=$MAMBA_USER:$MAMBA_USER src /opt/dvbfixer/src
RUN /opt/conda/envs/adapter/bin/python -m pip install --no-deps /opt/dvbfixer

COPY --chown=$MAMBA_USER:$MAMBA_USER \
    deploy/rfdiffusion-v1/container-entrypoint.sh /usr/local/bin/dvbfixer-rfdiffusion-v1

ENV DVBFIXER_RF_CHECKPOINT=/models/Base_ckpt.pt
WORKDIR /workspace
ENTRYPOINT ["/bin/sh", "/usr/local/bin/dvbfixer-rfdiffusion-v1"]
