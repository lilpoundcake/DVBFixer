# WP0: shared private chemical contract v1

Status: private WP0 implemented and checked on `feature/chemical-reconstruction`;
scientific reconstruction and production qualification remain incomplete.
Scope authorized by the active handoff: contract/admission/adapters/tests, not
coordinate generation, new native sampling, production routing or dependency changes.

## Decision before implementation

Use a private `dvbfixer.chemical_reconstruction` namespace with independent,
versioned request/admission types. Reuse the existing domain `ChemicalGraph`,
`ChemicalAtom`, `ChemicalBond`, `AuthorityProvenance`, `ObservedAtom`,
`ExternalLink` and exact structural identity. Do not change old dataclasses,
evidence readers, policies or serializers in place.

The old graph is a bounded chemical record, not the whole future solver:
ideal positions remain authority data with declared units; external covalent/
coordination connections and explicit geometry/stereo constraints are separate.
No rotor library or molecular toolkit is imported by this boundary. All new
contract/admission/adapter modules must be stdlib-only, including their imports.

Documentation order: this decision, private source/tests, active handoff/task-map
update and WP0 validation record (four documentation/map files total), plus a
narrow ignore rule for retained test artifacts. Pending links: none. User's implementation
request covers these bounded records; no generic capability is advertised as shipped.

## Immutable request and authority

- A component has exact identity, a declared chemistry category and zero or more
  explicit graph states. Zero states means unknown/not assessed, not complete;
  multiple states require an explicit state-id choice, never a priority/pH guess.
- Every state binds graph, authority provenance, ideal-coordinate units and a
  separate explicit/unspecified H-inventory flag. Presence of some H atoms does
  NOT prove the inventory is complete. A user-mapped graph digest must match
  the existing Class A canonical graph digest. Other provenance remains caller-
  supplied resolved metadata, not independently verified file/source authorship.
- Observations reuse domain atom identity/element/position with explicit units
  and occupancy. Identity spelling, component occurrence and source ordering
  are preserved. Coordinates are not normalized or converted inside the DTO.
- Generated, temporary and mobile masks are separate. Generated/temporary atoms
  must be absent observations; mobile atoms must already be observed and may be
  admitted only for explicit repacking. Completion never silently mobilizes an
  observation. Temporary atoms cannot be called scientific generated output.
- The request carries original source SHA and a deterministic revision-bound
  digest over ALL declared context, states, units, constraints, links and masks.
  This binds supplied content, not proof that a caller parsed the entire source
  file correctly. Admission never claims full-source inventory verification.

The representation must cover heavy completion, H placement and explicit
repacking as distinct operations. H placement requires a complete declared H
inventory and present/materialized heavy anchors; no protonation/state selection
is performed here. State pH ranges and graph priorities are retained but not used
to select alternatives or override an authoritative microspecies.

## Constraints and admission

Explicit distance, ring-closure, exclusion, angle, planarity, torsion and ordered
stereo constraints carry their units. Stereo is a specified orientation of
ordered neighbours, not an automatic translation of every chemical centre to
L or S. Labels on graph atoms must have an explicit corresponding stereo mapping;
absence is not a stereochemical pass. Future authority adapters must establish
the chemical correctness of that mapping (including CIP), not infer it from a
label alone at admission. Invalid references/units/masks refuse before dispatch.

Evaluate constraints whose referenced atoms are all fixed observations, using
independent Cartesian measurements only. Missing/mobile references defer the
check explicitly; there is no candidate generation or hidden observation repair.
This does not claim global geometry quality: it checks only declared constraints.

Capability descriptors are passed explicitly, NOT globally registered backends.
They declare operations, component categories, elements, constraint types and
link/context/mobility/temporary support. With no matching descriptor, report
unsupported even for a structurally valid request. An admitted request means
ONLY that this declared strategy contract matches — no solver has run.

Malformed immutable DTOs raise `ContractError` at construction. Well-shaped but
semantically incompatible requests return typed admission findings; adapters
raise `AdapterError` with a stable refusal code for transport they cannot perform.
These are not instructions to retry another physical engine. Findings and request
hashes are deterministic across hash seeds; bounded per-component inventories and
anchor counts avoid re-scanning the whole retained complex per capability.

Admission outcomes distinguish invalid input, unresolved/ambiguous chemistry,
unsupported capability and admitted planning. Coordinates-produced,
geometry-approved, MD-ready and source-coverage-verified remain false for every
WP0 outcome. Invalid or unresolved chemistry must not trigger a different engine.

Opaque unmodified environment components may remain present with an explicit
allow-opaque-context descriptor; list them as unassessed. Never drop their
observations/connections or claim their chemistry was reconstructed. Unresolved
target chemistry/link order is not an opaque-context exemption.

## Adapters and compatibility

- Protein adapter: caller supplies the intended graph state/authority; keep all
  old observations, exact missing identities, source SHA and declared links.
  Links lacking an authoritative order stay unresolved, not guessed peptide bonds.
  Context components without provided chemistry become explicit opaque records.
- Class A adapter: initially only the existing in-memory `UserMappedAuthority`.
  Verify its canonical graph digest, preserve locked state and source observations,
  and require a separate original-source SHA. Local/online unresolved authorities
  return a typed unsupported-adapter error; do not open CCD files/network. Reject
  linked target use rather than widen isolated Class A implicitly. Explicit full
  context must include every original observation unchanged.
  The initial transport envelope is deliberately narrower than legacy Class A:
  connected CHNO organic graphs, no redox/coordination metadata or known complex
  cofactor. Other organic elements are an explicit adapter refusal, not a claim
  that they are invalid Class A chemistry. This avoids copying the legacy solver's
  full classification knowledge into a second private implementation. Extracting
  a shared classification policy is a later compatibility extension; the generic
  contract itself has no CHNO-only or protein-only restriction.
- Old geometry-policy/source-request information is retained as origin metadata
  and request digest. It is not silently translated into a calibrated new solver
  policy. Passing new admission is not legacy scientific acceptance.
  The Class A origin stores the requested class and each old geometry-policy value
  explicitly; changing those values changes the origin and request binding even
  though WP0 does not execute or adopt the old solver's geometry policy.

Tests use synthetic graph authorities, mixed-complex context, H and heavy masks,
state ambiguity, fixed stereo conflicts, unit/case/icode/occurrence binding,
forged graph digests, unsupported links/classes and pure-import/no-I/O controls.
Old protein/Class A evidence and production routes remain unchanged. No numerical
experiment or empirical data is needed to implement WP0.

## Implemented boundaries and examples

| Owner | Role |
|---|---|
| [`contracts.py`](../../src/dvbfixer/chemical_reconstruction/contracts.py) | Immutable DTOs, authority-shape checks, canonical graph and source-bound request hashes |
| [`admission.py`](../../src/dvbfixer/chemical_reconstruction/admission.py) | Semantic masks/state/link checks and explicitly supplied capability matching |
| [`measurement.py`](../../src/dvbfixer/chemical_reconstruction/measurement.py) | Independent fixed Cartesian measurements; no generator or state chooser |
| [`adapters.py`](../../src/dvbfixer/chemical_reconstruction/adapters.py) | Existing protein and in-memory Class A transport; unresolved authority refuses |
| [`synthetic tests`](../../tests/test_chemical_reconstruction.py) | Contract controls, not scientific reconstruction claims |

Private planning usage, assuming a caller already has a typed legacy request
and an explicitly resolved `ChemicalState`:

```python
from dvbfixer.chemical_reconstruction.adapters import from_protein
from dvbfixer.chemical_reconstruction.admission import admit

planning_request = from_protein(legacy_request, state=resolved_state)
result = admit(planning_request)  # no declared strategy => unsupported
assert not result.coordinates_produced
assert not result.geometry_approved
```

**Do not** call admission with a synthetic capability and describe the result as
successful reconstruction. **Do** retain the original request and use the returned
findings, checked/deferred constraint indices and unassessed context to distinguish
planning compatibility from unexecuted scientific work. An unresolved protein link
becomes an explicit unknown connection; only caller-provided authoritative link
semantics can resolve it. No peptide order is inferred from spatial proximity.

Measurements normalize lengths to angstrom without mutating observations.
`stereo-center-frame` references `(centre, neighbour1, neighbour2, neighbour3)`
and measures the normalized scalar triple product. `stereo-four-ligands` references
`(centre, ligand1, ligand2, ligand3, ligand4)`: normalize centre-to-ligand vectors
first, then measure the normalized determinant of the differences from ligand4.
The sign and nondegeneracy threshold are caller-declared; this is not a CIP
implementation. Planarity uses the first three atoms as ordered plane anchors;
collinear anchors refuse. Torsion uses a signed degree interval in `[-180,180]`;
a wrapping/disjunctive interval is not represented by reversed bounds. Unsupported
scientific validation is deferred to later versioned validators, not silently
claimed by a planning success.

No source-file reader, byte publication, `.dat` writer, GUI, CLI, backend registry,
protonation selection, minimizer or reconstruction generator was added. The
[active handoff](chemical-reconstruction-handoff.md) remains the forward roadmap.

The [WP0 validation record](chemical-reconstruction-wp0-validation.json) binds
the final source/test hashes and retained checks. Its status, not the presence
of this design document, establishes which verification has actually completed.
