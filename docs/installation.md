# Installation

[← README](../README.md)

Create the environment and install:

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
```

**Modeller license:** The `model` command requires Modeller, which needs a free academic license key. Register at https://salilab.org/modeller/registration.html, then set the key in `<env>/lib/modeller-10.8/modlib/modeller/config.py`.

**Optional: PSI4 for free RESP charges.** The `parametrize -c resp --qm-engine psi4` path needs `psi4` + `psiresp` in a **separate conda env** (PSI4 ships its own BLAS/MKL stack that conflicts with OpenMM's; trying to install psi4 alongside dvbfixer's main env fails with libmamba "Could not solve" or pulls Python 3.9 which dvbfixer doesn't accept).

```bash
# One-time setup (~5 minutes):
micromamba create -n psi4 -c conda-forge psi4 psiresp
```

That's it. **Do NOT pip-install dvbfixer in the psi4 env** — dvbfixer stays in its own env and shells out to the psi4 env via `micromamba run -n psi4 python …` when `--qm-engine psi4` is invoked. Override the env name via `--psi4-env <name>` if you called it something else.

Only needed for `--qm-engine psi4`. The default AM1-BCC path and the Gaussian RESP path don't require it.

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
