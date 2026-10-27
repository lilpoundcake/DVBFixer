"""Versioned production identity for the multichain Protpardelle checkpoint."""

CUDA_PROFILE = "protpardelle-1c-cc94-cuda-v1"
MPS_PROFILE = "protpardelle-1c-cc94-mps-v1"
PROFILES = frozenset({CUDA_PROFILE, MPS_PROFILE})
CHECKPOINT_NAME = "cc94_epoch3100.pth"
CHECKPOINT_SHA256 = "a58c009ead22b6f79ef61d44e258b432b0fada1dfd8d44cb7a137d7459038a15"
CONFIG_SHA256 = "680096576f6bcec1da7e9c30a89c7278c04cf1ab16ed9a74c2cb99e37fa4db74"
ENGINE_REVISION = "ee378400f25b801fa481028000f9060183d7fb4c"
CUDA_ENVIRONMENT = "linux-amd64;cc94-ode-v1;python=3.12;torch=2.6.0+cu124;refinement=openmm-8.6.1-reference"
MPS_ENVIRONMENT = "macos-arm64;cc94-ode-v1;python=3.12;torch=2.6.0;refinement=openmm-8.5.2-cpu"
