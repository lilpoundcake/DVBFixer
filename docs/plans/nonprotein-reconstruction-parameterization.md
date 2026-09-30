# Nonprotein reconstruction and parameterization architecture plan

- **Status:** partially implemented. The offline Class A production slice is
  implemented as an opt-in typed Python application service; the CLI/GUI,
  Classes B/C, online enrichment, and MD backend gates remain unimplemented.
- **Scope:** reconstruct and prepare retained *nonprotein* components for a later,
  explicitly requested MD parameterization step. Components include isolated small
  molecules, glycans, PTMs, covalently attached ligands, cofactors, and metal
  systems.
- **Non-claim:** the implemented slice does not add a CLI option, remote service,
  force-field template, actual MD parameter generation, or MD-ready output. It does not replace PDBFixer,
  Salilab MODELLER, current ligand handling, or existing force-field adapters.
- **Core rule:** coordinate reconstruction/geometry regularization and MD
  parameterization are separate stages. Passing the former is never evidence that
  a molecule has valid MD charges, bonded terms, nonbonded parameters, redox
  treatment, or metal-coordination treatment.

## Goals and boundaries

### Implemented production slice (2026-09-30)

- `dvbfixer.domain.nonprotein_reconstruction` defines immutable backend-neutral
  requests, results, exact component/atom identities, authority modes, graph,
  microstate, geometry provenance, findings, and parameterization decisions.
- `dvbfixer.nonprotein_reconstruction` implements one opt-in Python boundary.
  It supports per-instance user-mapped complete graphs and a bounded JSON
  extraction from a checksum-pinned local CCD snapshot. The extraction format is
  `dvbfixer-ccd-snapshot-v1`; raw CCD monomer CIF ingestion is not claimed.
- V1 reconstructs only Class A known isolated components. It keeps every
  observed coordinate fixed, rigidly fits authority ideal geometry from at least
  three non-collinear observed anchors, generates absent atoms, and validates
  graph consistency, charge, bond geometry, stereochemistry, clashes, movement,
  and exact identity. Geometry success sets neither approval nor MD readiness.
- Local automatic state selection requires explicit pH and non-overlapping
  authority-declared pH ranges/priorities. User-mapped graphs are locked
  microspecies. Conflicts are `ambiguous`; missing coverage is `unsupported`.
- Class B links, Class C metal/redox/coordination/multi-component evidence, and
  online authority return explicit `unsupported` results before coordinate
  generation. They do not fall back to residue names, proximity, another source,
  or generic isolated-ligand chemistry.
- `decide_parameterization` consumes only a successful result plus independent
  approval and explicit route evidence. It reuses
  `domain.parameterization.classify_parameterization`, returns a route candidate,
  and always reports `md_ready=False`; it runs no force-field adapter.
- `reconstruct_and_publish` emits nothing for non-success. Successful geometry
  and provenance are staged, fsynced, and published together by one directory
  rename beneath an existing caller-owned root. Inputs, authority records, names,
  counts, nesting, pH, coordinates, and publication names are bounded.

- Build a fail-closed path from an incomplete observed nonprotein component to a
  validated geometry candidate, then—only when separately approved—to an MD
  parameterization decision.
- Make chemical authority, microstate choice, atom mapping, generated coordinates,
  validation, parameterization route, and publication fully reproducible.
- Preserve deposited information unless an explicit movable-atom policy permits a
  change, and report all permitted displacement.
- Preserve full structural identity throughout. The minimum atom key is
  `(model, chain, resid, icode, altloc, atom)`; residue/component-instance identity
  is `(model, chain, resid, icode, altloc)` plus the component occurrence. Chain
  IDs and insertion codes remain case-sensitive and must not be normalized.
- Preserve explicit covalent links and distinguish observed, reconstructed, and
  inferred information. Coordinate proximity alone is not chemical authority.
- Do not infer chemistry from residue names or incomplete coordinates alone.
  Bond order, aromaticity, formal charge, stereochemistry, protonation, redox
  state, attachment chemistry, and metal coordination must be established by
  allowed evidence or cause a safe refusal.
- Keep CIF handling at the existing CLI boundary. Scientific stages consume the
  normalized identity-preserving representation and must neither add ad-hoc CIF
  readers nor silently truncate identifiers to PDB limits.

## Explicit chemical-authority selection

- The user chooses exactly one authority mode for each reconstruction request;
  there is no silent source fallback:
  - **Local pinned CCD:** resolve only against a named, checksum-verified local
    Chemical Component Dictionary snapshot and its declared monomer/link data.
  - **Online API/enrichment:** resolve through a configured remote service or
    enrichment provider, subject to the security and reproducibility controls
    below. The exact response, endpoint identity, version/ETag where available,
    timestamp, and content digest become input provenance.
  - **User-mapped chemistry:** accept mapped SMILES, SDF, or monomer CIF supplied
    by the user. The mapping identifies the intended component instance(s), rather
    than applying a residue-name-wide guess.
- Source precedence is explicit, not an opportunistic retry order:
  - An accepted user-mapped graph and explicitly selected microstate are
    authoritative for their supported mapped component instance.
  - A local pinned CCD is authoritative only for the selected snapshot and exact
    component/link records.
  - Online enrichment is authoritative only for the captured response set of the
    explicit request; a later remote response cannot revise an already published
    result.
  - Explicit per-instance user overrides outrank automated decisions within the
    selected authority mode. Conflicts among equally authoritative records, an
    unmapped record, or a non-unique atom map are `ambiguous`, not a tie broken by
    heuristics.
- A source record must expose enough graph and state evidence for the requested
  class: atom and bond graph, bond orders, formal charges, stereochemistry where
  applicable, ideal/restraint geometry where needed, and explicit link definitions
  for attached components.
- User-mapped SMILES/SDF/monomer CIF remains authoritative only for an isolated,
  uniquely mapped, single-residue molecule unless the selected staged class and
  explicit link model support broader chemistry. Incompatible mappings,
  unexplained external bonds, and covalent attachment outside that scope fail
  rather than guess.

## Deterministic microstate policy

- Every automatic microstate decision requires an explicit pH and a declared
  state-selection policy; no hidden default pH is permitted.
- The selector receives the resolved chemical graph, allowed protonation/tautomer
  and stereochemical state space, explicit links, relevant local environment
  inputs, a deterministic algorithm version, and a declared seed when sampling is
  unavoidable.
- It produces a stable ranked candidate list and one selected state using a
  deterministic tie-break rule. The result records pH, units/conventions,
  algorithm and database versions, candidate identifiers and scores, tie-break
  reason, seed, and all input graph digests.
- An explicit user microstate override supersedes automatic selection and is
  recorded as an override, not silently relabelled as an automatic result. A
  user-supplied graph that already specifies a microspecies may lock that state;
  requesting automatic replacement of it requires an explicit opt-in.
- If state enumeration, scoring, compatibility with an attachment, or mapping to
  observed atoms is incomplete or ambiguous, return `unsupported` or `ambiguous`.
  Do not neutralize, protonate, deprotonate, alter stereochemistry, or fabricate a
  metal state by default.

## Staged component classes and fail-closed routing

- Classification happens before coordinate generation and is based on connected
  component inventory, complete full-identity links, selected chemical authority,
  graph completeness, and state evidence—not residue name alone.
- Proposed classes, in increasing evidence burden:
  - **Class A — known isolated organic component:** one residue, unique atom map,
    complete graph and microstate, no metal, and no external covalent bond outside
    the component. Eligible for constrained coordinate generation and later
    isolated-ligand parameterization evaluation.
  - **Class B — linked carbohydrate/PTM/covalent ligand:** complete component and
    link dictionaries, unique attachment mapping, linkage stereo, and a
    linkage-aware geometry objective. Not eligible for isolated-ligand fallback.
  - **Class C — cofactor, metal, redox, or multi-component assembly:** requires a
    validated complete graph, charge/redox state, coordination and link model, and
    an explicitly approved geometry and MD route. Generic GAFF-style isolated
    parameterization is prohibited.
  - **Class D — unknown, incomplete, ambiguous, or unsupported chemistry:** no
    reconstruction or parameterization. Return a structured refusal with the
    missing evidence and preserve input unchanged.
- Existing complex-cofactor protection remains a minimum guard: known complex
  cofactors must not be routed to generic isolated-ligand GAFF unless an explicit,
  validated user template covers the actual graph and links.
- An explicit strip/exclude request is a distinct route, not evidence that a
  retained component is parameterizable or that a protein-only result relaxed the
  interface.

## Proposed request, result, and provenance boundary

- Introduce backend-neutral typed request/result contracts before integrating any
  scientific adapter.
- A reconstruction request contains:
  - normalized structure coordinates and full atom/residue/component identity;
  - observed-atom masks, alternate-location and MODEL policy, explicit links, and
    fixed/movable coordinate policy;
  - selected authority mode, source locator/reference, expected checksums or
    captured online response, and mapped user inputs where applicable;
  - component class request, explicit pH/state policy and overrides;
  - candidate count, deterministic seed/policy, geometry tolerances, and output
    location/publication intent.
- A reconstruction result contains:
  - `succeeded`, `unsupported`, `ambiguous`, `invalid-input`, or `failed` status;
  - selected class, resolved graph/state/link model, unique observed-to-authority
    atom mapping, and full generated-atom provenance;
  - candidate coordinates, movement report for observed atoms, validation findings,
    warnings, backend/tool/version/dictionary provenance, and immutable input and
    output digests;
  - no public coordinate artifact for any non-success status.
- A separate parameterization request consumes only a successful, approved
  geometry result plus a concrete force-field family, compatible template/user
  parameter evidence, intended simulation conditions, and a parameterization
  policy. It returns a route decision and validation report; it does not mutate or
  reinterpret reconstruction provenance.
- Artifact publication is atomic: stage private files, validate all declared
  outputs and manifests, then publish coordinates, reports, and provenance
  together. Failure registers no partial visible result.

## Geometry stage: complete before MD evaluation

- **1. Inventory and authority resolution**
  - Enumerate connected nonprotein components and explicit links without losing
    full identity.
  - Resolve the user-selected authority source; verify pin/checksum or capture an
    online snapshot; validate source syntax and allowed component scope.
  - Classify the component and reject missing required evidence before invoking a
    coordinate engine.
- **2. State and graph resolution**
  - Apply deterministic pH-dependent microstate selection or a recorded user
    override.
  - Establish the exact atom/bond graph, bond orders, formal charges, stereo, and
    attachment definitions; require a unique mapping to observed atoms.
- **3. Coordinate materialization**
  - Generate only absent coordinates under constraints from retained observed
    coordinates, authority geometry, and explicit links.
  - Attribute every atom and coordinate to observed, source-template,
    algorithmically generated, or user-overridden provenance.
- **4. Local geometry regularization**
  - Regularize using graph- and link-aware geometry restraints while honoring the
    declared movable selection. Do not use a generic MD force field as a substitute
    for missing chemical authority or geometry reconstruction.
  - Validate graph, valence, bond/angle/planarity/stereo, attachment geometry,
    clashes, identity, and observed-atom displacement before publication.
- **5. Geometry approval gate**
  - A component passes to MD evaluation only if all class-specific validations and
    acceptance thresholds pass, the provenance is complete, and an authorized
    reviewer/policy approves the result. A geometry pass does not auto-launch MD
    parameterization.

## MD parameterization stage: separate decision and approval

- Re-inventory the approved geometry and verify that its component graph, links,
  identity keys, and provenance digest match the geometry-stage result.
- Route using an integrated decision service rather than residue-name heuristics:
  - explicit exclusion/strip;
  - validated user template covering the actual graph and external links;
  - exact native template after graph compatibility validation;
  - approved class-specific cofactor/metal system;
  - isolated-ligand candidate only for Class A; or
  - fail closed.
- Force-field matching, charge derivation, bonded/nonbonded term generation,
  cross-residue links, and system construction are MD-specific validations. They
  must report their own source/tool versions, inputs, warnings, and failures.
- No automatic downgrade from a failed explicit parameterization request to a
  stripped or generic route. A non-explicit optional route may decline only by
  returning a visible structured result that says the interface remains
  unrelaxed—not by claiming MD readiness.
- **MD approval gate:** publish a parameterized system only after template/graph
  compatibility, link coverage, charge/state consistency, force-field construction,
  and class-specific validation pass. Geometry artifacts remain independently
  available and correctly labelled if this gate fails.

## Validation, tests, and benchmark design

- Establish versioned fixtures with source provenance, input masks, expected
  class/status, authority source snapshot, pH/state decision, links, and withheld
  reference coordinates where available. Do not add generated production outputs
  unless a test requires a golden artifact.
- Unit tests cover:
  - case-distinct chain IDs, insertion codes, models, alternate locations, and
    identity-preserving atom maps;
  - source-mode exclusivity, source precedence, checksum/version pinning, and no
    silent offline/online/user-input fallback;
  - deterministic microstate ranking at explicit pH, tie handling, and override
    precedence/provenance;
  - graph/map ambiguity, unsupported chemistry, covalent attachment, metals,
    redox, incomplete records, and no-output-on-failure behavior;
  - class routing and prohibition of isolated-ligand treatment for linked and
    complex systems;
  - separation of a geometry-valid result from an MD-parameterized result.
- Integration tests cover:
  - isolated charged, aromatic, and stereogenic ligands with local CCD and
    user-mapped inputs;
  - glycans/PTMs/covalent ligands with explicit linkage records;
  - supported and intentionally refused cofactors/metals;
  - repeat-run byte/provenance stability in deterministic mode;
  - atomic publication, stale-manifest conflict handling where applicable, and
    absence of visible partial artifacts after adapter, validation, or network
    failure.
- Benchmark metrics are stratified by class and authority mode:
  - atom recovery and withheld-coordinate RMSD on observed anchors;
  - graph/link and stereochemical correctness;
  - bond, angle, planarity, attachment, clash, and observed-coordinate movement
    metrics;
  - correct `unsupported`/`ambiguous` classification;
  - MD template/system construction success only after the independent MD gate;
  - runtime, memory, cache behavior, repeatability, and external-service failure
    behavior.
- Fix numerical thresholds and expected classifications before comparing engines.
  Evaluate against declared PDBFixer/legacy and existing parameterization
  baselines; completing more inputs alone is not evidence of safer chemistry.

## Security, privacy, and licensing controls

- Treat user-supplied chemistry files and remote responses as untrusted input:
  enforce size, atom/bond/count, nesting, decompression, schema, parser-time, and
  subprocess-resource limits before chemistry processing.
- Online enrichment is opt-in only. Require HTTPS, a configured allowlist, bounded
  redirects/timeouts/retries, response-size limits, and private cache storage.
  Never put credentials in URLs, logs, provenance visible to other users, or
  localStorage. Redact authorization headers and secrets from diagnostics.
- Cache local and online authority records by content digest, source/version, and
  access scope; do not let an unpinned mutable response masquerade as a reproducible
  local dictionary.
- Validate all artifact paths and workspace ownership/authorization at the adapter
  boundary. Resolve inputs through approved artifact IDs or contained paths; never
  trust caller-supplied filesystem paths for publication.
- Keep parsers and external chemistry tools sandboxed/contained where deployment
  supports it; capture bounded diagnostics and clean private staging files on
  failure.
- Review source-code, data, service-term, binary, and model/weight licenses
  independently before adopting an adapter, CCD redistribution path, online API,
  force-field distribution, or benchmark fixture. Record the accepted license and
  redistribution obligations in release evidence.

## Delivery sequence and approval gates

- **Gate 0 — design approval:** approve the typed identity/provenance schema,
  authority-mode semantics, source precedence, class taxonomy, explicit pH/state
  policy, error taxonomy, security model, and licensing inventory. No production
  routing changes before this decision.
- **Gate 1 — offline evidence prototype:** implement no public behavior until the
  local pinned-CCD and user-mapped paths demonstrate complete provenance,
  fail-closed mapping, deterministic state selection, and geometry-only fixtures.
  Keep online enrichment out of the default path.
- **Gate 2 — geometry acceptance:** approve class-specific reconstruction results
  only after predeclared benchmark thresholds, identity/link preservation, no
  partial-publication failures, and independent scientific review pass.
- **Gate 3 — MD route acceptance:** separately demonstrate force-field graph/link
  coverage, charge/state compatibility, system construction, and explicit-route
  failures. Do not infer this approval from Gate 2.
- **Gate 4 — online operation approval:** approve remote enrichment only after
  reproducibility capture, offline behavior, cache isolation, credential handling,
  resilience tests, provider terms, and threat-model review pass.
- **Gate 5 — experimental exposure:** expose an opt-in experimental selector with
  complete provenance and no hidden fallback. Preserve existing supported
  PDBFixer, MODELLER, and parameterization behavior until a separate,
  evidence-backed product decision changes it.

Implementation status: Gate 0 and the bounded offline Class A portion of Gate 1
are represented by the typed Python service and tests. Gates 2-5 are not passed.
There is no public command selector, backend capability negotiation, accepted
benchmark threshold set, online adapter, B/C adapter, or MD output publication.

## Explicit non-goals

- Universal reconstruction of arbitrary unknown HETATM chemistry.
- Guessing covalent bonds, protonation, tautomers, redox state, stereochemistry,
  or metal coordination from coordinates or residue names.
- Treating CCD ideal geometry, RDKit/Open Babel output, UFF-style refinement, or a
  successful geometry optimization as MD-ready parameterization.
- Replacing existing production baselines, changing their defaults, or silently
  widening accepted chemistry while this plan is under evaluation.
