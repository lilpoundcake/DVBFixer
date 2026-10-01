# ADR 0010: Experimental diffusion gap-reconstruction boundary

Date: 2026-09-22

Status: Accepted.

## Context

DVBFixer is evaluating sequence-guided diffusion for internal protein gaps while
retaining deposited coordinates outside explicitly generated or movable regions.
The candidate engines require large, engine-specific Python/CUDA environments,
model checkpoints, and hardware that do not belong in the core DVBFixer runtime.
Their scores and failure modes are not interchangeable with MODELLER, and a
scientifically silent fallback would make provenance and benchmark conclusions
ambiguous.

PATCHR describes a relevant fixed-coordinate sampling strategy, but DVBFixer must
own its implementation and validation boundary. The repository must not imply
that proposed research has replaced the supported Salilab MODELLER, PDBFixer, or
`tleap-reduce` paths.

## Decision

- Diffusion engines run through a versioned subprocess protocol in an isolated
  workspace. Engine-specific Python packages, CUDA libraries, model checkpoints,
  and containers remain outside the core DVBFixer dependency environment.
- DVBFixer validates protocol manifests, paths, digests, identities, geometry,
  connectivity, clashes, and chirality independently of the external runner.
  Candidate artifacts are published only after every hard gate passes.
- Diffusion environment and model evidence is stored in a separate versioned
  provenance manifest. `DatRecord` remains the downstream preparation sidecar and
  does not absorb engine, checkpoint, device, or CUDA metadata.
- There is no automatic scientific fallback between diffusion and MODELLER.
  Unsupported input, unavailable engines, runner failure, and failed validation
  remain explicit outcomes rather than changing algorithms silently.
- Salilab MODELLER remains the default and supported `model` and `homology`
  backend. No public diffusion backend is added until the implementation plan's
  prerequisite phases and acceptance gates pass.
- DVBFixer implements the PATCHR-like constrained-sampling method independently.
  It does not import, vendor, copy, or wrap PATCHR as its implementation. PATCHR
  may be used only as an external scientific comparator unless a later ADR
  changes this policy.
- CIF normalization remains at the CLI boundary. The scientific diffusion stage
  consumes normalized PDB and preserves case-sensitive chain IDs, insertion
  codes, atom identities, and retained explicit links.

## Consequences

The CPU-testable contract, masking, validation, provenance, and publication
layers can be developed without an ML dependency or public CLI change. Each
external engine requires its own pinned revision, artifact hashes, license
inventory, and runner environment. GPU inference remains a separate benchmark
lane and cannot make the core test suite depend on CUDA availability.

Diffusion results cannot be compared to MODELLER by treating backend-native
scores as a shared scale. Promotion from proposed research to an experimental
public backend, and later to production support, require separate evidence-backed
decisions. This ADR accepts the integration boundary only; it does not claim that
any diffusion engine or user-facing diffusion capability is implemented.
