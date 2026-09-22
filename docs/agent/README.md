# Agent knowledge map

This directory is the compact, current-state index used by coding agents. It
answers four questions before an agent edits code:

1. Which symbol owns the behavior?
2. Which contracts and invariants constrain it?
3. Which files must change under one coordinating owner?
4. Which focused checks prove the change?

It is not a replacement for user documentation, generated CLI reference, or
accepted architecture decisions. Context pages explain only the details that
cannot be represented clearly in the structured maps.

## Source precedence

When sources disagree, use this order and record the contradiction:

1. Executable code.
2. Tests and tracked fixtures.
3. Hard rules in [`AGENTS.md`](../../AGENTS.md).
4. Accepted ADRs.
5. The structured maps in this directory.
6. Current-state architecture and command documentation.
7. Proposed roadmaps and research notes.
8. Historical design notes and changelog entries.

A roadmap is never evidence that a capability has shipped.

## Files

- [`tasks.toml`](tasks.toml) maps recurring tasks to entry symbols, callers,
  contracts, invariants, change groups, and focused tests.
- [`contracts.toml`](contracts.toml) records cross-stage artifacts and their
  producers, consumers, identity keys, side effects, and known divergences.
- [`invariants.toml`](invariants.toml) records enforceable rules, failure modes,
  evidence, and known gaps.
- [`contexts/`](contexts/) contains short explanations for Structure Identity,
  Structure Preparation, Force-field Naming, Parameterization, Topology
  Generation, Diagnostics, and Workflow Execution. Start with the page named by
  the task record rather than reading every context.
- [`docs/adr/`](../adr/) contains accepted architecture decisions. Unaccepted
  alternatives stay in [`docs/plans/`](../plans/) or
  [`docs/research/`](../research/).

The maps use TOML rather than YAML so validation works with Python 3.11's
standard `tomllib` and adds no dependency to the lightweight CI lane.

## Context index

| Context | Use it for |
|---|---|
| [Structure Identity](contexts/structure-identity.md) | Chain/residue/atom identity, PDB/mmCIF mapping, allocator policy |
| [Structure Preparation](contexts/structure-preparation.md) | Model, prepare, protonation, missing atoms, chirality, `.dat`, ZBS |
| [Force-field Naming](contexts/force-field-naming.md) | AMBER/CHARMM/GROMACS residue and atom names |
| [Parameterization](contexts/parameterization.md) | Native/user/GAFF routing and complex-cofactor guard |
| [Topology Generation](contexts/topology-generation.md) | RTP, ACPYPE, topology matching, output bundles |
| [Structure Diagnostics](contexts/diagnostics.md) | Report-only findings, reports, exit status, runtime summary boundary |
| [Workflow Execution](contexts/workflow-execution.md) | CLI, batch, CIF boundary, jobs, workspaces, API adapters |

Use `tasks.toml` as the primary router. This table is for browsing, not for
deciding which files may change independently.

## Status vocabulary

| Status | Meaning |
|---|---|
| `implemented` | Active production behavior with an owner and test evidence |
| `partial` | Some listed behavior is active, but the complete policy is not |
| `missing` | Required behavior has no implementation |
| `known-gap` | Current implementation is known to violate or omit the rule |
| `proposed` | Planned behavior; not available to callers |
| `research` | Investigated alternative with no delivery commitment |
| `deprecated` | Existing behavior retained temporarily and not for new callers |

Status belongs to each capability or invariant. A page-level status is only a
summary and must not hide a mixture of implemented and proposed behavior.

## Reference syntax

Code references use `path::symbol`, for example:

```text
src/dvbfixer/ffutils/ff_names.py::apply_variants_to_pdb_text
```

The path is validated automatically. Symbols are intentionally not resolved by
the first validator version because a cross-language AST index would add more
complexity than value. Prefer stable symbols over line numbers.

## Agent workflow

Before editing:

1. Find the task in `tasks.toml`.
2. Read its context, contracts, invariants, and change group.
3. Verify the map against current code and tests.
4. Report contradictions instead of silently choosing one source.

While editing:

1. Give one coordinating agent ownership of every `exclusive-owner` change
   group.
2. Parallelize evidence gathering and independent leaf files, not a shared
   cross-cutting contract.
3. Keep scientific policy in Python; transport and workspace adapters must not
   duplicate it.
4. Update a map only when ownership, contract, invariant, known gap, or focused
   verification changes.

After editing:

```bash
python scripts/check_agent_docs.py
git diff --check
```

Run the task's focused checks before broader suites.

## Delegation policy

Usually safe to parallelize:

- evidence inventories for unrelated contexts;
- tests for already-defined independent behavior;
- separate context pages;
- external-adapter research that does not alter shared contracts.

Require one coordinating owner:

- residue or atom identity;
- `.dat` schema or merge semantics;
- force-field variants and naming;
- CONECT inference;
- missing-atom rebuild and chirality;
- argparse, command registry, generated CLI docs, and GUI command schema;
- ZBS option propagation and intermediate placement;
- workspace manifest revisions and artifact registration;
- shared files in this directory.

## Maintenance rules

- Keep context pages under roughly 200 lines unless a real contract requires
  more detail.
- Do not copy generated CLI help.
- Do not copy full implementation algorithms into these files.
- `implemented` requires at least one owner and one test path.
- Proposed names and types must be clearly marked `proposed`.
- Update `verified_on` and `verified_at_commit` only after checking the cited
  code and tests.
- Run `python scripts/check_agent_docs.py` in CI and locally.
