# ADR 0011: Experimental diffusion CLI

Date: 2026-10-02

Status: Accepted.

## Context

ADR 0010 accepted an isolated, versioned diffusion boundary but deliberately did
not accept a user-facing backend. Since then, the narrow canonical-protein scope,
independent validator, atomic directory publication, Apple MPS runner, and
portable Linux/CUDA runner have been implemented. The frozen 231-case Apple M3
Pro portability evaluation completed with fallback disabled, and the public
protocol wrapper passed a native smoke. The Linux wrapper still requires
checkpoint-backed acceptance on the pinned NVIDIA host.

The two runners do not have equivalent sampling semantics or evidence. Protenix
uses fixed-coordinate reinjection after every denoising update and was selected
through the frozen confirmatory benchmark. Protpardelle uses template
conditioning followed by final fixed-coordinate restoration; its training-data
membership remains unresolved and its evidence is descriptive. Treating these
profiles as interchangeable would overstate the Apple result.

## Decision

- `dvbfixer model --backend diffusion` is an opt-in experimental public CLI.
  MODELLER remains the default and production-supported backend.
- The public profiles are explicit and never selected automatically:
  - `protenix-v1-cuda` targets the pinned Linux/NVIDIA runner and records
    per-step fixed-coordinate reinjection;
  - `protpardelle-1c-mps` targets the pinned native Apple Silicon runner, records
    final restoration, and must not claim per-step reinjection.
- There is no automatic fallback between diffusion, profiles, devices, or
  MODELLER. Unsupported scope, failed preflight, runner failure, and validation
  failure terminate without publishing output.
- The initial public scope is one canonical protein chain with one unambiguous
  two-anchor internal gap of 3-12 residues. Heterogens, noncanonical chemistry,
  terminal or multiple gaps, unsupported links, batch, GUI, `zbs`, and
  `homology` remain outside this promotion.
- Diffusion output is a new directory bundle committed with one same-parent,
  no-replace rename after independent validation. Candidate `.dat` files remain
  engine-neutral; engine, checkpoint, device, and validation evidence belongs
  in separate provenance artifacts.
- Engine packages, ML frameworks, platform runtimes, and checkpoints remain
  outside the core environment. DVBFixer does not download or redistribute
  checkpoints.
- The Apple profile remains experimental/descriptive with unresolved training
  membership. The Linux profile remains experimental until the production
  wrapper reproduces the accepted frozen scope on the pinned NVIDIA host.
- Every successful candidate references a contained, digest-verified sampler
  trace artifact. The trace records profile semantics, atom/fixed-mask digests,
  callback evidence, final restoration, refinement, device/fallback state, and
  bounded resources; atomic publication copies it into the candidate bundle.
  `doctor` reports immutable profile status and can perform a bounded,
  no-inference handshake with a selected production wrapper.
- An operator-supplied runner is trusted executable code. Core path checks
  contain and validate returned artifacts but do not form an OS/filesystem
  sandbox and cannot prevent writes elsewhere. Maintained wrappers are
  contractually required to write only inside their dedicated workspace;
  deployments requiring enforcement must add an external sandbox or container.

## Consequences

Users can explicitly invoke the narrow diffusion path while existing MODELLER
commands retain their file-output contract and behavior. Core CI remains free of
Torch, CUDA, MPS, engine packages, and checkpoints; dependency-light protocol
tests use fakes and mocks, while hardware acceptance stays in separate lanes.

This decision does not promote diffusion to the default, deprecate MODELLER,
authorize arbitrary chemistry, prove Apple training independence, accept the
Linux hardware profile, or add diffusion to batch, GUI, `zbs`, or `homology`.
Future scope expansion and production support require separate evidence-backed
decisions.
