# dvbfixer doctor — Capability Report

[← command index](index.md) · [← README](../../README.md)

Reports installed Python backends, external chemistry executables, and
available OpenMM platforms without starting a structure-processing job.
Modeller is probed in an isolated subprocess, so an installed package with an
invalid or missing license is reported as unavailable without crashing Doctor.

```bash
dvbfixer doctor
dvbfixer doctor --format json
```

Missing optional tools are reported rather than treated as an error. Use the
JSON form in CI or before submitting a large folder-input run.

## Batch mode

`doctor` does not support directory batch input. It checks the DVBfixer
installation once and has no per-structure input. See the
[batch support matrix](../batch-mode.md#support-by-tool).
