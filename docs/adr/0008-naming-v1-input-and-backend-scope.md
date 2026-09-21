# ADR 0008: Naming V1 input and prep-backend scope

Date: 2026-09-21

Status: Accepted.

## Context

The workspace API already provides bounded upload, download, artifact ownership,
and path containment. Naming V1 accepts exact residue overrides inline, while
large reusable files would require a second artifact schema and lifecycle.
The opt-in `tleap-reduce` backend emits AMBER-family names after all tleap work,
but its early-return adapters previously skipped the public `--atom-naming`
output policy.

The pure naming matrix uses minimal synthetic PDB records because byte-level
fixed-column behavior is easier to review there. Existing tracked protein,
nucleic-acid, cap, insertion-code, and regression fixtures cover pipeline input
classes; no missing behavior currently justifies another generated golden file.

## Decision

- Naming V1 remains workspace-artifact-only. Clients use the existing workspace
  upload and download routes; the naming route does not duplicate direct file
  transfer.
- Variant overrides remain an inline bounded array. V1 does not accept an
  override artifact. A future reusable-override format requires its own schema,
  media type, provenance, and ADR.
- `tleap-reduce` honors `--atom-naming` for its final public output in both
  `prepare` and `protonate`. Conversion runs after tleap, Reduce, variant-H
  patching, and chirality inspection. The backend remains AMBER-only, so the
  target naming family is AMBER; `gromacs` applies shifts and `standard` keeps
  native AMBER atom names.
- Add new golden fixtures only when a real structure exposes behavior that the
  independent mapping matrix and current tracked fixtures cannot express.

## Consequences

The naming endpoint stays small and has one artifact transaction model. Inline
overrides remain easy to validate and include in provenance. The two prep
backends now agree on the documented naming option without feeding GROMACS names
back into tleap. `.dat` atom identities are computed after final naming and
therefore match the emitted PDB.
