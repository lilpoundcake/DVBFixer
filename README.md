# dvbfixer

Prepare structural biology files from PDB and CIF inputs: repair missing atoms,
rebuild protein gaps, assign protonation states, minimize structures, and export
GROMACS topologies. Use the command line or the browser-based GUI with Mol*.

Current release: **0.10.3**.

## Install

Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html),
Git and a Python **3.11+** bootstrap interpreter. The installer offers two profiles:

- **Standard**: scientific CLI, including PySCF/RESP, MDAnalysis/Plotly and MAFFT.
- **Desktop**: everything in Standard plus Node and the built browser GUI.

```bash
git clone https://github.com/lilpoundcake/DVBFixer.git
cd DVBFixer
python3 scripts/install_dvbfixer.py
```

Follow the prompts to choose your profile, environment manager and environment
name. The manager determines its location; no installation path entry is needed.
For Desktop, you can also opt into an application shortcut.
Actions, profiles and managers are selected from numbered lists.

The helper installs the package and dependencies, builds Desktop when selected,
offers hidden MODELLER license entry, and runs functional checks. No separate
`pip install` or npm setup is needed. MODELLER requires an appropriate license
for gap rebuilding and comparative modeling; deferring it leaves action required
(exit **2**), not a fully licensed installation.

See [installation](docs/installation.md#standarddesktop-source-checkout-helper)
for environment activation, explicit installation paths, supported manager/platform
versions, automation, licensing and optional add-ons. Native Apple Silicon users
should use micromamba **2.9.0** or mamba
**2.10.0**; native conda has unresolved upstream package collisions. Manual
environment files are an [alternative for development](docs/installation.md#manual-installation-and-source-development),
not additional installation steps. Platform evidence and remaining gates are
in the [acceptance record](docs/installation-profile-acceptance.md).

### Update an existing installation

From your source checkout:

```bash
git pull
python3 scripts/install_dvbfixer.py --update
```

Choose your installation from the numbered list. The helper shows the
installed and checkout versions before confirmation. Managed installations retain
their profile/add-ons and rebuild Desktop's GUI when applicable. Older/manual
installations offer an explicitly confirmed package-only update, leaving their
scientific dependencies unchanged. Close that installation's
GUI and jobs first. See [updating](docs/installation.md#update-an-existing-installation)
for non-interactive commands and manual environments.

### Add a profile or GUI to an existing environment

```bash
python3 scripts/install_dvbfixer.py --extend
```

Select the environment by number, then choose **Standard**, **Desktop**, or
**GUI only**. The same menu also offers **Diffusion backend**, provisioned in a
separate runtime rather than adding PyTorch to the scientific environment.
The helper previews profile/GUI dependency changes before confirmation. GUI-only
keeps the existing scientific stack; incompatible full-profile plans are refused
rather than forced. See [adding components](docs/installation.md#add-a-profile-or-gui-to-an-existing-environment).

## First structure

In your active installed environment, repair atoms and hydrogens for a protein
with no missing sequence regions:

```bash
dvbfixer prepare input.pdb -o prepared.pdb
dvbfixer diagnose prepared.pdb
```

This writes `prepared.pdb` and a `prepared.dat` sidecar for later minimization.
Inspect the diagnostic findings before using the structure downstream.
Preparation is not a guarantee of MD readiness for arbitrary ligands or cofactors.

For gap rebuilding and minimization, follow the [quick start](docs/quick-start.md)
or choose an [end-to-end recipe](docs/pipelines.md).

## Open the GUI

Open the installer-created Desktop application shortcut—no terminal activation
needed. Alternatively, from your active Desktop environment, run:

```bash
dvbfixer gui
```

The command opens the browser after the local server is ready. Stop with Ctrl+C
or `dvbfixer gui --stop`. See [Desktop launch](docs/commands/gui.md) for details.

Follow the [GUI guide](docs/gui.md) to import files, run tools, and inspect
results. Only contributors need the separate [source-development setup](gui/README.md#source-development).

## User manual

| I want to… | Read |
|---|---|
| Find my starting point | [User guide](docs/README.md) |
| Process my first structure | [Quick start](docs/quick-start.md) |
| Prepare proteins, glycans, or topology inputs | [Recipes](docs/pipelines.md) |
| Find a particular tool | [Command index](docs/commands/index.md) |
| Choose alternate conformations | [Altloc selection](docs/commands/altlocs.md) |
| Look up exact CLI options | [CLI reference](docs/reference/index.md) |
| Process a folder or save logs | [Batch mode and logging](docs/batch-mode.md) |
| Choose a force field | [Force fields](docs/force-fields.md) |
| Understand CIF conversion | [CIF input](docs/cif-input.md) |
| Resolve an error or limitation | [Troubleshooting](docs/known-issues.md) |

## Developers and administrators

- [Contributing](CONTRIBUTING.md)
- [Software citation](CITATION.cff)
- [GUI setup and server configuration](gui/README.md)
- [API](docs/api.md) and [public deployment](docs/deployment.md)
- [Architecture](ARCHITECTURE.md) and [domain policies](docs/domain-model.md)
- [Agent instructions](AGENTS.md) and [task map](docs/agent/README.md)

Scientific inputs are tracked in [tests/fixtures](tests/fixtures/README.md).
After installing development extras with `pip install -e ".[dev]"`, run
`pytest -m 'not slow' -q`. The full suite needs the external scientific tools
declared in the selected profile. Optional engine tests require their separate
runtimes; keep PyTorch outside Standard. The real PSI4 test uses an explicit
`DVBFIXER_TEST_PSI4_PREFIX`.

## License

For distribution scope and delivery rules, see
[licensing and distribution policy](docs/licensing-status.md).

DVBFixer's original code is available under [BSD-3-Clause](LICENSE).
Copyright (c) 2026, Stanislav Skripkin. Commercial use, modification and
redistribution are permitted under that license; author notification is not
required.

Third-party libraries, bundled force fields, datasets and model weights retain
their own terms and are not relicensed by this grant. See
[third-party notices](THIRD_PARTY_NOTICES/README.md),
[model terms](docs/model-licenses.md) and the
[recorded dependency review](docs/research/third-party-license-audit.md).
The project license does not establish permission for every combined
distribution or restricted scientific backend.

The primary delivery route is source plus the existing Standard/Desktop installer,
which obtains dependencies into the user's environment under their own terms.
This licensing work does not publish a prebuilt scientific environment/container.
Bundled FF sets remain included with upstream attribution; wheel support remains
available but is not a complete scientific-stack installer.
