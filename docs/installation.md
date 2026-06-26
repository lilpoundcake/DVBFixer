# Installation

[← README](../README.md)

Create the environment and install:

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
```

**Modeller license:** The `model` command requires Modeller, which needs a free academic license key. Register at https://salilab.org/modeller/registration.html, then set the key in `<env>/lib/modeller-10.8/modlib/modeller/config.py`.

**Optional: PSI4 for free RESP charges.** The `parametrize -c resp --qm-engine psi4` path needs `psi4` + `psiresp`. These ship their own BLAS/MKL stack that conflicts with OpenMM's, so they are NOT in the main `environment.yml`. Install them either in a SEPARATE env (recommended) or carefully alongside (may need libmamba solver, ~10 min):

```bash
# Recommended — separate env:
micromamba create -n dvbfixer-psi4 -c conda-forge psi4 psiresp
micromamba run -n dvbfixer-psi4 pip install -e .

# Or single-env (fragile):
micromamba install -n dvbfixer -c conda-forge psi4 psiresp
```

Only needed if you use `--qm-engine psi4`. The default AM1-BCC path and the Gaussian RESP path don't require it.

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
