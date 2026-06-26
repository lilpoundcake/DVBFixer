# Installation

[← README](../README.md)

Create the environment and install:

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
```

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
