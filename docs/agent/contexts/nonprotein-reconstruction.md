# Nonprotein Reconstruction

Status: partial

Verified on: 2026-09-30

Verified at commit: `11611b8f688325afcf831c6269f3aaafe47eba5e`

## Purpose

This context owns the opt-in geometry-only reconstruction boundary for one exact
nonprotein component and its handoff to a separate parameterization decision.

## Scope

V1 supports known isolated Class A chemistry through typed Python APIs. It does
not alter prepare/minimize defaults or expose a CLI/GUI command. Class B linked
chemistry, Class C complex chemistry, and online authority are explicit refusals.

## Capabilities

- Exact model, chain, residue, insertion-code, altloc, occurrence, and atom identity.
- Exclusive pinned-local or per-instance user-mapped authority.
- Deterministic explicit-pH microstate records.
- Fixed-observed geometry generation, validation, provenance, and atomic bundles.
- A separate candidate route decision that never claims MD readiness.

## Entry Points

- `src/dvbfixer/nonprotein_reconstruction.py::reconstruct_nonprotein`
- `src/dvbfixer/nonprotein_reconstruction.py::decide_parameterization`
- `src/dvbfixer/nonprotein_reconstruction.py::reconstruct_and_publish`
- `src/dvbfixer/domain/nonprotein_reconstruction.py::ReconstructionRequest`

## Contracts

- One request carries one authority mode; there is no fallback across sources.
- Local input is a bounded `dvbfixer-ccd-snapshot-v1` JSON extraction with an
  exact SHA-256 pin, not an unrestricted raw CCD parser.
- User-mapped chemistry is a complete, uniquely named graph for one exact
  component instance. Existing ligand chemistry remains the owner of any
  SMILES/RDKit valence interpretation before constructing that graph.
- Non-success has no coordinates or output digest and cannot be published.

## Invariants

- Observed atoms do not move in V1.
- Geometry success leaves `geometry_approved=False` and `md_ready=False`.
- Classes B/C and online authority refuse before coordinate generation.
- Parameterization consumes the successful geometry digest and explicit approval.

## Callers

There are no production CLI, GUI, prepare, minimize, or topology callers. The
typed API is the smallest public integration and is covered directly by tests.

## Adapters

The local snapshot loader and atomic directory publisher are the only adapters.
There is no online client, raw monomer-CIF adapter, reconstruction subprocess,
or force-field generator in this boundary.

## Side Effects

Reconstruction and route decisions are pure except for reading a selected local
snapshot. Publication stages two JSON files under the caller-owned root, fsyncs
them, and renames the complete directory. Failed operations clean staging and do
not publish a visible result.

## Known Divergences

- No whole-structure connected-component inventory exists yet.
- Raw CCD CIF, mapped SMILES/SDF parsing, link dictionaries, planarity/angle
  restraint sets, environment-aware state scoring, and backend negotiation are absent.
- Geometry uses authority ideal coordinates and a rigid anchor fit; no local
  optimization engine is integrated.
- Route selection does not run or validate actual MD parameter generation.

## Proposed Work

Gate any CLI/GUI exposure on benchmark acceptance. Add B/C or online adapters
only with their complete chemistry, threat model, provenance, and refusal tests;
never widen Class A fallback to cover them.

## Focused Verification

```bash
PYTHONPATH=src pytest -q tests/test_nonprotein_reconstruction.py tests/test_scientific_domain.py
PYTHONPATH=src mypy src/dvbfixer/nonprotein_reconstruction.py \
  src/dvbfixer/domain/nonprotein_reconstruction.py \
  src/dvbfixer/domain/structure_identity.py
python scripts/check_agent_docs.py
```
