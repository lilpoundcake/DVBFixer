# Nonprotein Reconstruction

Status: partial

Verified on: 2026-10-06

Verified at commit: `de2069369494a6e71cd17d02888ec6d6beb9f57a`

## Purpose

This context owns the opt-in geometry-only reconstruction boundary for one exact
nonprotein component and its handoff to a separate parameterization decision.

## Scope

V1 supports known isolated Class A chemistry through typed Python APIs and an
opt-in single-component PDB CLI. It does not alter prepare/minimize defaults or
provide a dedicated GUI workflow. Class B linked chemistry, Class C complex chemistry, and
online reconstruction authority are explicit refusals.

## Capabilities

- Exact model, chain, residue, insertion-code, altloc, occurrence, and atom identity.
- Exclusive pinned-local or per-instance user-mapped authority.
- Deterministic explicit-pH microstate records.
- Fixed-observed geometry generation, validation, provenance, and atomic bundles.
- A separate candidate route decision that never claims MD readiness.
- Bounded local/explicit-online CCD information lookup with content digest.
- Local raw CCD heavy-atom conversion and atomic whole-PDB/provenance bundles.

## Entry Points

- `src/dvbfixer/nonprotein_reconstruction.py::reconstruct_nonprotein`
- `src/dvbfixer/nonprotein_reconstruction.py::decide_parameterization`
- `src/dvbfixer/nonprotein_reconstruction.py::reconstruct_and_publish`
- `src/dvbfixer/ccd.py::fetch_ccd_component`
- `src/dvbfixer/pdb_component_reconstruction.py::reconstruct_pdb_component`
- `src/dvbfixer/reconstruct_component.py::main`
- `src/dvbfixer/domain/nonprotein_reconstruction.py::ReconstructionRequest`

## Contracts

- One request carries one authority mode; there is no fallback across sources.
- Local input is a bounded `dvbfixer-ccd-snapshot-v1` JSON extraction with an
  exact SHA-256 pin, not an unrestricted raw CCD parser.
- User-mapped chemistry is a complete, uniquely named graph for one exact
  component instance. Existing ligand chemistry remains the owner of any
  SMILES/RDKit valence interpretation before constructing that graph.
- Non-success has no coordinates or output digest and cannot be published.
- Online CCD lookup is read-only discovery; reconstruction requires an explicit
  local CIF and does not fall back across authority sources.
- Whole-PDB materialization adds only absent heavy atoms and preserves every
  source ATOM/HETATM line byte-for-byte.

## Invariants

- Observed atoms do not move in V1.
- Geometry success leaves `geometry_approved=False` and `md_ready=False`.
- Classes B/C and online authority refuse before coordinate generation.
- Parameterization consumes the successful geometry digest and explicit approval.

## Callers

The opt-in `component-info` and `reconstruct-component` CLIs are production
callers and are discoverable through the generated generic command schema.
There are no dedicated GUI, prepare, minimize, ZBS, or topology callers.

## Adapters

Adapters include the local snapshot loader, bounded raw component-CIF parser,
explicit read-only HTTPS lookup, strict PDB instance mapper/materializer, and
atomic directory publishers. There is no online reconstruction authority,
reconstruction subprocess, or force-field generator in this boundary.

## Side Effects

Reconstruction and route decisions are pure except for reading a selected local
snapshot. Publication stages two JSON files under the caller-owned root, fsyncs
them, and renames the complete directory. Failed operations clean staging and do
not publish a visible result.

## Known Divergences

- No whole-structure connected-component inventory exists yet.
- Mapped SMILES/SDF parsing, link dictionaries, planarity/angle restraint sets,
  environment-aware state scoring, and backend negotiation are absent.
- Geometry uses authority ideal coordinates and a rigid anchor fit; no local
  optimization engine is integrated.
- Route selection does not run or validate actual MD parameter generation.

## Proposed Work

Gate any CLI/GUI exposure on benchmark acceptance. Add B/C or online adapters
only with their complete chemistry, threat model, provenance, and refusal tests;
never widen Class A fallback to cover them.

## Focused Verification

```bash
PYTHONPATH=src pytest -q tests/test_nonprotein_reconstruction.py \
  tests/test_ccd_component_reconstruction.py tests/test_scientific_domain.py
PYTHONPATH=src mypy src/dvbfixer/nonprotein_reconstruction.py \
  src/dvbfixer/domain/nonprotein_reconstruction.py \
  src/dvbfixer/domain/structure_identity.py
python scripts/check_agent_docs.py
```
