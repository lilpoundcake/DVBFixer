# Design notes (migrated)

> This file was ~85 KB of per-subcommand implementation notes moved
> out of `CLAUDE.md` in Phase 4b of the revision plan. All of its
> content has since migrated to more targeted destinations:

- **Per-subcommand "How it works" prose** →
  [`commands/`](commands/), one page per subcommand under a
  `## How it works` heading.
- **CLI dispatcher, GROMACS topology, GLYCAM integration, PDB format
  notes** → [`../ARCHITECTURE.md`](../ARCHITECTURE.md).
- **`Package Structure` and `Environment & Installation`** → replaced
  by [`../CLAUDE.md`](../CLAUDE.md) (scannable index) and
  [`installation.md`](installation.md).

Nothing lives here anymore. Kept as a stub so git-history links from
old commits still resolve to a real file that explains where the
content went.

For readers looking at old commits or discussions that reference this
file: run `git log --follow docs/DESIGN_NOTES.md` — the migration
tranches are commits `38cddb7` (initial move from CLAUDE.md),
`b8c9657` (Notes → ARCHITECTURE), and the follow-up that split the
per-subcommand prose across `docs/commands/*.md`.
