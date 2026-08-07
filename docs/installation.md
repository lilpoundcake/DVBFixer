# Installation

[← README](../README.md)

Create the environment and install:

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
```

`environment.yml` pins `python >=3.11,<3.14`. The upper bound is required:
propka 3.5.1 (used by `protonate` / `prepare`) reads the dataclass attribute
`self.__annotations__` at the instance level, which Python 3.14's PEP 649/749
change makes raise `AttributeError`, crashing the PROPKA step. Do not loosen
the cap.

## Multiple-sequence alignment executables

`dvbfixer msa`, `dvbfixer homology --template-plan`, and the GUI Homology
workspace run an external alignment
program. The full `environment.yml` installs all three supported engines:

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
dvbfixer msa --list-engines
```

For an existing conda/micromamba environment, install one or more engines
from conda-forge/bioconda. MAFFT is the default and is sufficient by itself:

```bash
micromamba install -n dvbfixer -c conda-forge -c bioconda mafft
# Optional alternatives:
micromamba install -n dvbfixer -c conda-forge -c bioconda "muscle>=5" clustalo
```

The commands must be available on `PATH` under these exact names:

| Engine | Required executable | Check |
|---|---|---|
| MAFFT | `mafft` | `mafft --version` |
| MUSCLE 5 | `muscle` | `muscle -version` |
| Clustal Omega | `clustalo` | `clustalo --version` |

On macOS, Homebrew users can install the native packages with
`brew install mafft muscle clustal-omega`. On Debian/Ubuntu, distribution
packages are commonly available as `mafft`, `muscle`, and `clustalo`, but
verify that the `muscle` package is version 5: MUSCLE 3 uses a different CLI
and is not supported. Native Windows users should use WSL with the Linux
instructions, or download upstream binaries and put their directory on
`PATH`. When installing a standalone MUSCLE binary, rename it to `muscle`
(`muscle.exe` on Windows); rename a standalone Clustal Omega binary to
`clustalo` (`clustalo.exe` on Windows).

After any standalone installation, restart the GUI/dev server so it inherits
the updated `PATH`, then confirm detection with:

```bash
dvbfixer msa --list-engines
```

`--engine auto` selects the first installed engine in this order: MAFFT,
MUSCLE 5, Clustal Omega. See the official download pages for
[MAFFT](https://mafft.cbrc.jp/alignment/software/),
[MUSCLE 5](https://drive5.com/muscle5/manual/install.html), and
[Clustal Omega](http://www.clustal.org/omega/).

**macOS Docker (VirtioFS) users:** if your `MAMBA_ROOT_PREFIX` is under a host
bind mount (e.g. `/home/agent`, a `fakeowner` mount of `/Users`), `micromamba
create` will abort at `Linking 'ncurses'` with
`filesystem error: cannot copy symlink: Invalid argument` — the bind mount is
case-insensitive and rejects ncurses's case-variant terminfo symlinks, and no
`always_copy`/`--copy` flag fixes it. Create the env on the container's native
overlay filesystem instead:

```bash
sudo mkdir -p /opt/mamba && sudo chown -R agent:agent /opt/mamba
export MAMBA_ROOT_PREFIX=/opt/mamba
micromamba create -f environment.yml -n dvbfixer -y
micromamba run -n dvbfixer pip install -e ".[dev]"
```

Persist `MAMBA_ROOT_PREFIX=/opt/mamba` and keep `/opt/mamba/envs/dvbfixer/bin`
on `PATH` in your shell rc. See [known issues](known-issues.md) for the full
diagnosis.

**Modeller license:** The `model` command requires Modeller, which needs a free academic license key. Register at https://salilab.org/modeller/registration.html, then set the key in `<env>/lib/modeller-10.8/modlib/modeller/config.py`.

**Free RESP charges (PySCF)** — already included in `environment.yml`:

PySCF is pip-installed automatically by `environment.yml`. It enables the recommended free RESP path:

```bash
dvbfixer parametrize molecule.pdb -c resp --qm-engine pyscf -v
```

The default AM1-BCC charge method needs no QM engine and is the right choice for ~95% of small-molecule use cases.

**Optional: PSI4 backend for RESP (alternative to PySCF).** The `--qm-engine psi4` path is supported but requires a **separate conda env** (PSI4 ships its own BLAS/MKL stack that conflicts with OpenMM's, and PSI4 conda-forge builds pull Python 3.9 which dvbfixer doesn't accept):

```bash
# One-time setup (~5 minutes):
micromamba create -n psi4 -c conda-forge psi4 psiresp
```

**Do NOT pip-install dvbfixer in the psi4 env** — dvbfixer stays in its main env and shells out to the psi4 env via `micromamba run -n psi4 python …` when `--qm-engine psi4` is invoked. Override the env name via `--psi4-env <name>`.

Fragile on macOS arm64 (libint2 SONAME mismatches). On macOS, prefer `--qm-engine pyscf`.

After installation, `dvbfixer` is available as a CLI command:

```bash
dvbfixer <command> [options]
```

Or without activating the environment:

```bash
micromamba run -n dvbfixer dvbfixer <command> [options]
```

## Next steps

- Browse the [command index](commands/index.md)
- Read the [pipelines](pipelines.md) page for end-to-end recipes
- Skim the [known issues](known-issues.md)
