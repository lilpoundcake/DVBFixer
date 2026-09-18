# LLM-oriented domain documentation roadmap

Completion status is tracked in the canonical
[implementation checklist](implementation-checklist.md).

Status: implemented; future additions are demand-driven maintenance.

This roadmap applies selected Domain-Driven Design practices to create a compact
project index for coding agents. The documentation is optimized for navigation,
change safety, delegation, and verification rather than for teaching DDD or
describing an aspirational architecture.

## Objective

An agent starting a change should be able to answer quickly:

1. Which symbol owns the behavior?
2. Which callers and artifacts depend on it?
3. Which invariants and known gaps constrain the change?
4. Which files require one coordinating owner?
5. Which focused checks prove the change?

The map must distinguish shipped behavior from vocabulary, proposals, research,
and known gaps. It must not present a conceptual bounded context as a completed
code migration.

## Implemented foundation

The canonical agent entry point is [`../agent/README.md`](../agent/README.md).
The first slice contains:

- [`../agent/tasks.toml`](../agent/tasks.toml): task-to-symbol navigation and
  coupled change groups;
- [`../agent/contracts.toml`](../agent/contracts.toml): PDB naming and `.dat`
  variant contracts, producers, consumers, side effects, and divergences;
- [`../agent/invariants.toml`](../agent/invariants.toml): naming and identity
  rules with status, enforcement, evidence, and exceptions;
- [`../agent/contexts/force-field-naming.md`](../agent/contexts/force-field-naming.md):
  the first bounded-context explanation;
- `scripts/check_agent_docs.py`: dependency-free validation;
- `tests/test_agent_docs.py`: repository-level validation coverage.

The maps use TOML because Python 3.11 includes `tomllib`. Adding a YAML parser
solely for documentation would complicate the lightweight CI lane.

## Source precedence

The complete rule is in the agent README. In summary:

1. Executable code.
2. Tests and tracked fixtures.
3. Hard rules in `AGENTS.md`.
4. Accepted ADRs.
5. Agent maps.
6. Current-state documentation.
7. Proposed roadmaps and research.
8. Historical notes.

Contradictions are recorded as gaps; agents must not silently choose a roadmap
over production behavior.

## Deliberate reductions

The earlier plan proposed a large glossary, a separate context map, a separate
invariant catalog, six long context pages, and several initial ADRs. That shape
would duplicate `AGENTS.md`, `ARCHITECTURE.md`, `docs/domain-model.md`, command
guides, and API plans.

The revised design deliberately avoids:

- a repository-wide DDD migration;
- one Python package per conceptual context;
- manually maintained sequence diagrams;
- line-number references that drift after refactors;
- ADRs for decisions that have not been accepted;
- a representation matrix for every glossary term;
- generated human documentation before the structured maps prove useful;
- cross-language AST symbol validation in the first checker version.

Context pages are added only when a task or contract needs explanation beyond
the structured records.

## Record model

### Tasks

Each recurring task records:

- context page;
- starting symbols;
- active callers;
- contracts and invariants;
- coupled change group;
- focused test files and commands.

Tasks are added on demand. The goal is not to predict every possible change.

### Contracts

Each cross-stage artifact records:

- authoritative owner;
- producers and consumers;
- identity keys;
- source mutation policy;
- commit and failure behavior;
- tests and known divergences.

High-value contracts include PDB identity, `.dat`, CONECT/LINK/SSBOND, CIF
normalization, Homology template plans, topology outputs, and workspace
manifests.

### Invariants

Each invariant has a stable ID, per-claim status, owner, enforcement points,
failure mode, tests, and known exceptions. A page-level `partial` label never
substitutes for these individual statuses.

### Change groups

Cross-cutting files use either:

- `exclusive-owner`: one agent coordinates all edits; research may still be
  delegated;
- `parallel-safe`: independent agents may edit the listed leaves.

Identity, variant naming, `.dat`, CONECT, chirality, generated CLI surfaces,
ZBS propagation, and workspace revisions default to `exclusive-owner`.

## Context page template

Every context page uses these sections, validated by the checker:

1. Purpose
2. Scope
3. Capabilities
4. Entry Points
5. Contracts
6. Invariants
7. Callers
8. Adapters
9. Side Effects
10. Known Divergences
11. Proposed Work
12. Focused Verification

Pages should remain around 100-200 lines and link to authoritative records
instead of copying them.

## Delivery phases

### Phase 0: format and governance

Status: complete.

- Define source precedence and status vocabulary.
- Choose dependency-free TOML maps.
- Define reference syntax and delegation rules.
- Add validation to the lightweight CI lane.

Exit criterion: the maps parse, references resolve, cross-record IDs exist, and
implemented records have owner/test evidence.

### Phase 1: Force-field Naming vertical slice

Status: complete.

- Map the current in-place naming helper and its callers.
- Separate naming conversion from canonical `rename` and parameterization.
- Record PDB and `.dat` contracts.
- Record insertion-code, explicit-variant, idempotence, collision, and
  structure-preservation invariants.
- Record known backend and parser gaps without presenting proposed API behavior
  as shipped.

Exit criterion: an agent can scope a naming change without first reading the
entire architecture document.

### Phase 2: existing domain modules

Status: complete.

Two agents may independently research Structure Identity and Parameterization.
One integrator updates shared TOML maps.

Required outputs:

- `contexts/structure-identity.md`;
- `contexts/parameterization.md`;
- task, contract, invariant, and change-group records supported by active code;
- explicit lists of exported vocabulary with no production consumers.

Exit criterion: every export from `dvbfixer.domain` has meaning, active callers
or an explicit unused/vocabulary-only status, and tests.

### Phase 3: artifact contracts

Status: complete.

Evidence gathering may run in parallel for:

- complete `.dat` lifecycle;
- PDB/mmCIF identity normalization;
- CONECT/LINK/SSBOND;
- workspace manifest and artifact identity;
- Homology template-plan/PIR/mosaic;
- topology output directories.

One integrator writes `contracts.toml`. Each record must include mutation,
atomicity, partial-output, and cleanup behavior in addition to schema shape.

Exit criterion: agents can follow an artifact from producer through consumer
and understand failure side effects.

### Phase 4: pipeline contexts

Status: complete.

Independent context-page packets:

- Structure Preparation;
- Topology Generation;
- Diagnostics;
- Workflow Execution.

Agents own one page each and do not edit shared TOML files. The integrator adds
cross-context records after reconciling contradictions.

Exit criterion: pages identify boundaries and contracts without reproducing
pipeline algorithms.

### Phase 5: task expansion

Status: complete for the current high-risk task set; further entries are
demand-driven.

Likely records:

- change `.dat` schema;
- change missing-atom rebuild;
- change CONECT inference;
- change topology matching;
- change workspace persistence;
- add API endpoint;
- add fixture.

Exit criterion: only repeated, high-risk work patterns are indexed.

### Phase 6: stronger validation

Status: baseline complete; advanced symbol/freshness automation remains
optional.

The checker now validates TOML syntax, metadata, statuses, path and Markdown
references, task-to-contract/invariant/change-group links, required evidence,
required context sections, and unindexed context pages.

Add the following only when maintenance evidence justifies it:

- symbol existence checks;
- stale `verified_at_commit` reporting;
- generated human views;
- pull-request checks for changed contracts or invariants.

Do not implement a general Python/TypeScript AST index as part of the initial
documentation rollout.

## Subagent execution model

Use this pattern for every new context:

1. Parallel read-only evidence agents inspect policy, callers, tests, and
   adapters.
2. One writer owns the context page.
3. One integrator owns shared TOML maps and contradiction resolution.
4. One reviewer verifies every `implemented` claim against code and tests.

Subagents return unresolved contradictions rather than choosing silently. They
must not edit generated `docs/reference/*.md` or shared map files unless assigned
as the sole owner.

## Validation

Run after every map change:

```bash
python scripts/check_agent_docs.py
pytest -q tests/test_agent_docs.py
git diff --check
```

Context-specific focused tests remain listed in `tasks.toml` and the context
page. Documentation validation checks structure and references; it does not
claim to validate scientific truth automatically.

## Definition of done

- A new agent finds the owner, callers, contracts, invariants, known gaps, and
  focused tests for the indexed task without broad repository discovery.
- Every `implemented` record has code ownership and test evidence.
- Proposed behavior is visibly separate from current behavior.
- Cross-cutting changes declare a parallelization policy.
- Artifact records include failure and filesystem side effects.
- Shared maps have one coordinating editor per change.
- Links and referenced paths pass `scripts/check_agent_docs.py`.
- `docs/domain-model.md`, `ARCHITECTURE.md`, and the agent map do not contradict
  the current limited DDD integration.

## Related documents

- [Agent knowledge map](../agent/README.md)
- [DVBFixer API roadmap](api-roadmap.md)
- [API and DDD continuation context](api-ddd-continuation-context.md)
- [Current scientific domain model](../domain-model.md)
- [Repository architecture](../../ARCHITECTURE.md)
