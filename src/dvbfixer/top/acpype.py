"""``dvbfixer top --acpype`` — OpenMM + ParmEd + ACPYPE pipeline.

Alternative to the RTP-based topology build. Runs when the user
passes ``--acpype``; ignores ``--ff`` (always AMBER14+GLYCAM),
``--water`` (always TIP3P), ``--ignh``, and ``--merge``. Respects
``--ss``, ``--his``, ``--protonate`` (per-residue subset), and
``--keep-all-hydrogens``.

Delegates the heavy lifting to :func:`dvbfixer.acpype_export.export_gromacs`,
which handles glycoprotein preprocessing (CYX renaming, GLYCAM H
stripping, chain reordering), the OpenMM → ParmEd handoff, and the
ACPYPE call that produces the GROMACS ``topol.top`` + ``.gro`` with
``[ pairs_nb ]`` for mixed AMBER/GLYCAM 1-4 scaling.

Split out of ``top/pipeline.py`` in the Phase 2.4 follow-up work.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# AMBER variant names that OpenMM's addHydrogens understands. The
# --protonate / --his flags accept CHARMM-style names too; the
# _NAME_TO_AMBER lookup normalises them before the AMBER-variants
# validation gate.
_AMBER_VARIANTS = {"ASH", "GLH", "HIE", "HID", "HIP", "CYX", "LYN"}

_NAME_TO_AMBER = {
    "ASPP": "ASH", "ASPH": "ASH",
    "GLUP": "GLH", "GLUH": "GLH",
    "HSP": "HIP", "HSE": "HIE", "HSD": "HID",
}


def _parse_ss_pairs(ss_specs: list[str]) -> set[tuple[str, int]]:
    """Parse ``--ss CHAIN1:NUM1:CHAIN2:NUM2`` specs into a flat set
    of ``(chain, resseq)`` residues to force CYX renaming on.
    """
    extra_ss: set[tuple[str, int]] = set()
    for ss_spec in ss_specs:
        parts = ss_spec.split(":")
        if len(parts) == 4:
            extra_ss.add((parts[0], int(parts[1])))
            extra_ss.add((parts[2], int(parts[3])))
    return extra_ss


def _parse_protonation_overrides(
    protonate_spec: str | None,
    his_specs: list[str],
) -> dict[tuple[str, int], str]:
    """Parse ``--protonate`` (comma-separated CHAIN:NUM[:STATE]) and
    ``--his CHAIN:NUM:STATE`` into a merged ``{(chain, resseq):
    amber_variant}`` map.

    ``--protonate all`` isn't a residue selector — the caller handles
    that mode by feeding every ASP/GLU/HIS into the pipeline; this
    function returns an empty dict for it. Unknown / malformed specs
    are silently dropped (matches the historical inline behaviour).
    """
    prot_overrides: dict[tuple[str, int], str] = {}

    if protonate_spec and protonate_spec != "all":
        if ":" not in protonate_spec:
            print(
                f"ERROR: Invalid --protonate value '{protonate_spec}'.",
                file=sys.stderr,
            )
            sys.exit(1)
        for spec in protonate_spec.split(","):
            parts = spec.split(":")
            if len(parts) == 3:
                state = parts[2].upper()
                state = _NAME_TO_AMBER.get(state, state)
                if state in _AMBER_VARIANTS:
                    prot_overrides[(parts[0], int(parts[1]))] = state

    for his_spec in his_specs:
        parts = his_spec.split(":")
        if len(parts) == 3:
            state = parts[2].upper()
            state = _NAME_TO_AMBER.get(state, state)
            if state in _AMBER_VARIANTS:
                prot_overrides[(parts[0], int(parts[1]))] = state

    return prot_overrides


def run_acpype_mode(input_path: Path, args: argparse.Namespace) -> None:
    """Execute the ``--acpype`` branch of ``dvbfixer top``.

    Called from ``top.pipeline.main`` when ``args.acpype`` is set;
    it takes over the whole invocation and returns after writing
    the ACPYPE output tree. RTP-based topology building is skipped
    in this mode.
    """
    from dvbfixer.acpype_export import export_gromacs

    if args.ff == "charmm":
        print(
            "WARNING: --acpype always uses AMBER14+GLYCAM, ignoring --ff charmm",
            file=sys.stderr,
        )

    extra_ss = _parse_ss_pairs(args.ss)
    prot_overrides = _parse_protonation_overrides(args.protonate, args.his)

    if args.output:
        out_dir = Path(args.output).parent or Path(".")
        basename = Path(args.output).stem
    else:
        out_dir = input_path.parent or Path(".")
        basename = input_path.stem

    export_gromacs(
        input_path,
        out_dir,
        basename=basename,
        extra_ss=extra_ss or None,
        prot_overrides=prot_overrides or None,
        verbose=args.verbose,
        keep_all_hydrogens=args.keep_all_hydrogens,
    )
