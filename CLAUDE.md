# CLAUDE.md

Compatibility entry point for coding harnesses that discover this filename.
[`AGENTS.md`](AGENTS.md) is the canonical repository instruction file; read it
before making changes. If this file and `AGENTS.md` disagree, follow
`AGENTS.md` and report the drift.

## Project

**dvbfixer** is a Python package and CLI for preparing PDB/PDBx structural
biology files, with a React/Node workspace in `gui/`. The current preparation
backends are:

- `legacy` (default): broad-input PDBFixer/OpenMM-Modeller path supporting the
  input classes documented by `prepare`.
- `tleap-reduce` (opt-in): deterministic pure-protein AmberTools/Reduce path
  that rejects unsupported chemistry.

PDBFixer and Salilab MODELLER are supported production dependencies. Research
into atom-reconstruction, loop-modeling, homology-modeling, diffusion, or
geometry-regularization alternatives does not deprecate them and must not be
described as shipped behavior.

## Required workflow

Use the DDD agent map rather than reading the whole codebase first:

1. Find the task in [`docs/agent/tasks.toml`](docs/agent/tasks.toml).
2. Read the referenced context, contracts, invariants, change group, and focused
   checks.
3. Verify those records against the relevant current code and tests.
4. Report contradictions instead of silently choosing one source.
5. Give one coordinating owner to each `exclusive-owner` change group.
6. After editing, run the task's focused checks, then:

```bash
python scripts/check_agent_docs.py
git diff --check
```

Source precedence and status vocabulary are defined in
[`docs/agent/README.md`](docs/agent/README.md). A roadmap or research note is
never evidence that a capability has shipped.

## Navigation

| Need | Read |
|---|---|
| Canonical agent rules | [`AGENTS.md`](AGENTS.md) |
| Task/contract/invariant routing | [`docs/agent/README.md`](docs/agent/README.md) |
| Commands and quick start | [`README.md`](README.md), [`docs/commands/`](docs/commands/) |
| Current module/data-flow architecture | [`ARCHITECTURE.md`](ARCHITECTURE.md) |
| Implemented scientific policies and limits | [`docs/domain-model.md`](docs/domain-model.md) |
| Installation and external tools | [`docs/installation.md`](docs/installation.md) |
| Known issues and workarounds | [`docs/known-issues.md`](docs/known-issues.md) |
| Accepted architecture decisions | [`docs/adr/`](docs/adr/) |
| Proposed plans and research | [`docs/plans/`](docs/plans/), [`docs/research/`](docs/research/) |
| Generated exact CLI help | [`docs/reference/`](docs/reference/) — do not edit by hand |
| Fixture provenance/checksums | [`tests/fixtures/README.md`](tests/fixtures/README.md) |

[`docs/DESIGN_NOTES.md`](docs/DESIGN_NOTES.md) is a historical migration stub,
not a current source of implementation guidance.

## Dangerous-to-miss rules

These are only the compact cross-cutting subset. The task-specific rules in
`AGENTS.md` and `docs/agent/` still apply.

- Preserve case-sensitive chain IDs and insertion codes. Atom identity is
  `(chain, resid, icode, atom)`; residue/variant identity is
  `(chain, resid, icode)`. Never normalize chain IDs by case.
- Preserve `keepIds=True` on every `PDBFile.writeFile` call.
- Use `dvbfixer.ffutils.dat.DatRecord` for `.dat` files; do not hand-roll JSON
  readers or writers.
- CIF normalization belongs at the CLI boundary. Do not add scientific-stage
  CIF readers or silently truncate identifiers/data to PDB limits.
- Never call `fixer.addMissingAtoms()` directly. Use
  `ffutils.geometry.rebuild_missing_atoms_with_retry`, and preserve the legacy
  preparation order documented in the Structure Preparation context.
- Never call `PDBFixer.addMissingHydrogens`. The legacy path uses variant-aware
  OpenMM `Modeller.addHydrogens` followed by `repair_misplaced_hydrogens`; the
  `tleap-reduce` path uses Reduce.
- A stage claiming hard L-chirality output must run `assert_all_l` after its
  final heavy-atom coordinate change. The DDD map records where that guarantee
  is still partial; do not overstate standalone backend coverage.
- Supplied ligand SMILES is authoritative only for its supported mapped,
  isolated, single-residue molecule. Fail rather than guess on incompatible,
  ambiguous, or covalently attached chemistry.
- Do not reimplement force-field variant handling, ligand-valence chemistry,
  GLYCAM system construction, CONECT policy, or shared sequence placement in a
  local pipeline.
- Generated CLI reference and GUI command schema must follow argparse changes:
  run `python scripts/gen_cli_reference.py` and
  `python scripts/gen_gui_spec.py`; never hand-edit `docs/reference/*.md`.
- Keep credentials out of URLs, logs, manifests, and localStorage. Preserve
  workspace containment, authorization, revision, and atomic-publication
  contracts.
- Keep reviewed structural inputs under `tests/fixtures/`, document provenance,
  and regenerate `tests/fixtures/MANIFEST.sha256`. `/test/` is local scratch.
- For development `minimize`/`zbs` iterations, always use `--no-solvent`.
  Use `--strip-heterogens` when the iteration is not explicitly evaluating
  retained heterogens, ligands, glycans, cofactors, or their interfaces.

## Environment and verification

Python remains `>=3.11,<3.14`; NumPy remains `<2.5` while mypy targets Python
3.11. On macOS Docker/VirtioFS, place the micromamba root on native container
overlay storage (for example `/opt/mamba`), not the host bind mount. Full setup
and rationale live in [`docs/installation.md`](docs/installation.md) and
[`docs/known-issues.md`](docs/known-issues.md).

Common checks:

```bash
pytest -m 'not slow' -q
ruff check src/dvbfixer
mypy src/dvbfixer/cli.py src/dvbfixer/ffutils src/dvbfixer/pdbutils src/dvbfixer/align.py
python scripts/check_agent_docs.py
python scripts/gen_cli_reference.py --check
python scripts/gen_gui_spec.py --check
```

For GUI changes, from `gui/` run `npm run typecheck` and `npm test -- --run`.
Prefer the task record's narrow checks before broad suites.
