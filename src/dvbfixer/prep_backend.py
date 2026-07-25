"""Deterministic protein prep backend built on AmberTools + MolProbity.

Replaces the old ``PDBFixer.addMissingAtoms`` + ``Modeller.addHydrogens``
pair, both of which had known upstream bugs:

* ``PDBFixer`` places CB on the D face of the CA-N-C plane for
  backbone-only residues (openmm/pdbfixer#145 — a template-based
  internal-coordinate rebuild with no chirality constraint).
* ``Modeller.addHydrogens`` produces coincident H atoms on the same
  methyl/methylene/NH3+ parent (0.4-1.5 Å apart) — internal H-only
  minimize fails to separate them. LJ 1/r^12 → NaN in the next
  minimize.

The new pipeline:

    strip H → tleap (heavy atoms + all H per L-only templates) →
    strip H from tleap output → reduce -build -nuclear (deterministic
    H + HIS tautomers + ASN/GLN flips) → PROPKA + ``assign_amber_variants``
    (rename residues for ASH/GLH/LYN/CYM/CYX) → output

``tleap`` is deterministic and L-only by construction. ``reduce`` is
deterministic; it picks HIS tautomers based on H-bond environment and
rewrites the residue name column (HIE/HID/HIP). Neither tool produces
coincident atoms.

Both tools drop chain IDs and renumber residues on output, so this
module preserves the input's ``(chain, resseq, icode)`` metadata by
matching residues by ordinal position and re-emitting the PDB with
the original identifiers.

Ships with ``ambertools>=23`` (which provides both ``tleap`` and
``reduce``).
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


class PrepBackendError(RuntimeError):
    """Base for tleap / reduce subprocess failures with useful context."""


class TleapError(PrepBackendError):
    """tleap subprocess failed. ``log`` holds the full stdout+stderr so
    the caller can surface the offending residue name to the user."""

    def __init__(self, message: str, log: str) -> None:
        super().__init__(f"{message}\n---- tleap log ----\n{log}")
        self.log = log


class ReduceError(PrepBackendError):
    """reduce subprocess failed."""


# ---------------------------------------------------------------------------
# Binary discovery
# ---------------------------------------------------------------------------


def _require_binary(name: str) -> str:
    """Return absolute path to ``name`` or raise with an install hint."""
    path = shutil.which(name)
    if not path:
        raise PrepBackendError(
            f"'{name}' not found on PATH. Install AmberTools >= 23 "
            f"(conda install -c conda-forge ambertools) which ships "
            f"both tleap and reduce."
        )
    return path


# ---------------------------------------------------------------------------
# Metadata preservation — tleap and reduce both drop chain IDs / renumber
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ResidueMeta:
    """Original ``(chain_id, resseq, icode, resname)`` for one residue,
    keyed by ordinal position in the input file. Used to restore metadata
    after tleap/reduce strip it."""

    chain: str
    resseq: str
    icode: str
    resname: str


def _capture_residue_meta(pdb_path: Path) -> list[ResidueMeta]:
    """Walk the ATOM/HETATM records of ``pdb_path`` and record one
    ``ResidueMeta`` per residue in file order.

    Uniqueness is determined by ``(chain, resseq, icode)`` tuple —
    each unique tuple appears once regardless of how many atoms it has.
    TER records reset the "last seen" so N-terminal residues of a new
    chain don't merge with the C-terminal of the previous one when the
    resseq happens to match.
    """
    meta: list[ResidueMeta] = []
    last_key: tuple[str, str, str] | None = None
    for raw in pdb_path.read_text().splitlines():
        if raw.startswith("TER"):
            last_key = None
            continue
        if not raw.startswith(("ATOM  ", "HETATM")):
            continue
        if len(raw) < 27:
            continue
        chain = raw[21]
        resseq = raw[22:26].strip()
        icode = raw[26].strip()
        resname = raw[17:20].strip()
        key = (chain, resseq, icode)
        if key != last_key:
            meta.append(ResidueMeta(chain=chain, resseq=resseq,
                                    icode=icode, resname=resname))
            last_key = key
    return meta


def _restore_metadata(
    tleap_pdb: Path,
    out_pdb: Path,
    meta: list[ResidueMeta],
) -> None:
    """Rewrite chain ID (col 21), resseq (col 22-25), icode (col 26) in
    each ATOM/HETATM line of ``tleap_pdb`` using ``meta[residue_ordinal]``.

    tleap emits residues in the same order it read them; a new residue
    starts whenever the resseq column changes or a TER record is seen.
    We track residue ordinal by that same rule and index into ``meta``.
    """
    lines = tleap_pdb.read_text().splitlines(keepends=True)
    out_lines: list[str] = []
    residue_idx = -1
    last_tleap_key: tuple[str, str] | None = None

    for line in lines:
        if line.startswith("TER"):
            # TER emitted after last residue; carry residue_idx unchanged
            # (tleap emits TER between chains and at end).
            last_tleap_key = None
            out_lines.append(line)
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            out_lines.append(line)
            continue
        if len(line) < 27:
            out_lines.append(line)
            continue

        tleap_resseq = line[22:26].strip()
        tleap_resname = line[17:20].strip()
        key = (tleap_resseq, tleap_resname)
        if key != last_tleap_key:
            residue_idx += 1
            last_tleap_key = key

        if residue_idx >= len(meta):
            # More residues than we captured — leave the line as-is.
            out_lines.append(line)
            continue

        m = meta[residue_idx]
        # Rewrite cols 21 (chain), 22-25 (resseq), 26 (icode).
        # Pad resseq into a 4-char right-justified field; icode is a
        # single char (blank if empty).
        resseq_field = f"{m.resseq:>4s}"[:4]
        icode_field = m.icode[:1] if m.icode else " "
        new_line = line[:21] + m.chain + resseq_field + icode_field + line[27:]
        out_lines.append(new_line)

    out_pdb.write_text("".join(out_lines))


# ---------------------------------------------------------------------------
# H-strip (used before tleap and before reduce; tleap's H isn't kept)
# ---------------------------------------------------------------------------


def _strip_hydrogens(pdb_path: Path, out_path: Path) -> None:
    """Copy ``pdb_path`` to ``out_path`` dropping every ATOM/HETATM line
    whose element (cols 77-78) is ``H`` or whose atom name starts with
    ``H``. Non-ATOM lines pass through unchanged.

    Uses element column first (authoritative); falls back to atom-name
    prefix for files that don't fill the element column.
    """
    out: list[str] = []
    for raw in pdb_path.read_text().splitlines(keepends=True):
        if not raw.startswith(("ATOM  ", "HETATM")):
            out.append(raw)
            continue
        elem = raw[76:78].strip() if len(raw) >= 78 else ""
        name = raw[12:16].strip()
        if elem == "H" or (not elem and name and name[0] == "H"):
            continue
        # Also catch names like "HA", "HB2", "HG1" etc. when element
        # column present but wrong; but keep atoms named e.g. "HG" that
        # are element Hg (mercury) if elem says so.
        if elem and elem != "H":
            out.append(raw)
            continue
        if name and name[0] == "H":
            continue
        out.append(raw)
    out_path.write_text("".join(out))


# ---------------------------------------------------------------------------
# tleap subprocess wrapper
# ---------------------------------------------------------------------------


def run_tleap(
    input_pdb: Path,
    output_pdb: Path,
    ff: str = "leaprc.protein.ff14SB",
    extra_leaprc: list[str] | None = None,
    verbose: bool = False,
) -> None:
    """Run ``tleap`` on ``input_pdb`` and write the L-only heavy-atom-complete
    result to ``output_pdb``.

    Preserves chain IDs, resseq, and insertion codes from the input by
    capturing metadata before tleap runs and re-injecting it into tleap's
    output (tleap otherwise drops chain IDs and renumbers residues).

    Raises :class:`TleapError` on subprocess failure; the exception's
    ``.log`` attribute has the full tleap stdout+stderr so callers can
    surface the offending residue name to the user.
    """
    tleap_bin = _require_binary("tleap")

    meta = _capture_residue_meta(input_pdb)

    with tempfile.TemporaryDirectory(prefix="dvbfixer_tleap_") as tmp:
        tmpdir = Path(tmp)
        script = tmpdir / "tleap.in"
        raw_out = tmpdir / "raw.pdb"

        leaprc_lines = [f"source {ff}"]
        if extra_leaprc:
            leaprc_lines.extend(f"source {x}" for x in extra_leaprc)

        script.write_text("\n".join([
            *leaprc_lines,
            f"mol = loadpdb {input_pdb}",
            f"savepdb mol {raw_out}",
            "quit",
            "",
        ]))

        proc = subprocess.run(
            [tleap_bin, "-s", "-f", str(script)],
            capture_output=True, text=True,
        )
        log = (proc.stdout or "") + (proc.stderr or "")
        if proc.returncode != 0 or not raw_out.exists():
            raise TleapError(
                f"tleap exited with code {proc.returncode}; no output PDB.",
                log,
            )
        if "FATAL" in log or "Failed to generate parameters" in log:
            raise TleapError("tleap logged FATAL / parameter failure.", log)

        if verbose:
            print(f"  [prep] tleap OK: {input_pdb.name} -> "
                  f"{output_pdb.name} ({len(meta)} residues)")

        _restore_metadata(raw_out, output_pdb, meta)


# ---------------------------------------------------------------------------
# reduce subprocess wrapper
# ---------------------------------------------------------------------------


def run_reduce(
    input_pdb: Path,
    output_pdb: Path,
    build: bool = True,
    nuclear: bool = True,
    verbose: bool = False,
) -> None:
    """Run MolProbity ``reduce`` on ``input_pdb`` and write H-complete
    PDB to ``output_pdb``.

    * ``build=True``  → ``-build`` (optimize HIS tautomers + ASN/GLN flips).
    * ``nuclear=True`` → ``-nuclear`` (neutron-derived X-H bond lengths).

    Reduce is deterministic and does NOT move heavy atoms; it only adds
    or repositions H atoms. It rewrites the residue name column to
    HID / HIE / HIP based on the picked tautomer, so the output is a
    directly AMBER-compatible PDB.

    Reduce preserves chain IDs when they're present in the input. If
    the input's chain column is blank (e.g. straight tleap output),
    the output will also be blank — the caller is responsible for
    re-injecting chain IDs upstream (usually by calling ``run_tleap``
    with metadata restoration first).
    """
    reduce_bin = _require_binary("reduce")

    args = [reduce_bin]
    if build:
        args.append("-build")
    if nuclear:
        args.append("-nuclear")
    args.extend(["-quiet", str(input_pdb)])

    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode not in (0, 1):
        # reduce exits 1 on success too (writes to stdout, "success" via
        # a diagnostic on stderr). Only >1 is a real failure.
        raise ReduceError(
            f"reduce exited with code {proc.returncode}: {proc.stderr}"
        )
    if not proc.stdout:
        raise ReduceError(
            f"reduce produced no output PDB. stderr:\n{proc.stderr}"
        )
    output_pdb.write_text(proc.stdout)

    if verbose:
        print(f"  [prep] reduce OK: {input_pdb.name} -> "
              f"{output_pdb.name}")


# ---------------------------------------------------------------------------
# AMBER variant assignment (HIS handled by reduce; the rest by PROPKA)
# ---------------------------------------------------------------------------


# pKa thresholds — protonate ↔ deprotonate transition point per residue.
# Standard values from PROPKA papers; used only to decide the variant
# label. Actual H placement is by reduce.
_STD_PKA = {
    "ASP": 3.86,   # ASH if PROPKA > pH (protonated)
    "GLU": 4.34,   # GLH if PROPKA > pH
    "LYS": 10.34,  # LYN if PROPKA < pH (deprotonated)
    "CYS": 8.55,   # CYM if PROPKA < pH (deprotonated)
}


def _classify_variant(
    resname: str, pka: float | None, ph: float, is_ss: bool,
) -> str | None:
    """Return the AMBER variant name for this residue, or None if the
    residue should stay as its standard name.

    HIS is not handled here (reduce writes HID/HIE/HIP directly into
    the residue name column). SS-bonded CYS → CYX (highest priority,
    overrides pKa-based CYM classification).
    """
    if resname == "CYS":
        if is_ss:
            return "CYX"
        if pka is not None and pka < ph:
            return "CYM"
        return None
    if resname == "ASP":
        if pka is not None and pka > ph:
            return "ASH"
        return None
    if resname == "GLU":
        if pka is not None and pka > ph:
            return "GLH"
        return None
    if resname == "LYS":
        if pka is not None and pka < ph:
            return "LYN"
        return None
    return None


def _rename_residues_in_pdb(
    pdb_path: Path,
    positions: dict[tuple[str, int], str],
    target: str,
) -> int:
    """Rewrite the resname column (17-20) of every ATOM/HETATM whose
    ``(chain, resseq)`` is in ``positions`` to ``target`` (3-char).

    Used to temporarily rename CYX/CYM → CYS before the second tleap
    pass; tleap tolerates CYS but rejects CYX without an explicit
    ``bond`` declaration.
    """
    if not positions:
        return 0
    n = 0
    out: list[str] = []
    for raw in pdb_path.read_text().splitlines(keepends=True):
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 27:
            out.append(raw)
            continue
        chain = raw[21]
        try:
            resseq = int(raw[22:26].strip())
        except ValueError:
            out.append(raw)
            continue
        if (chain, resseq) not in positions:
            out.append(raw)
            continue
        out.append(raw[:17] + f"{target:>3s}" + raw[20:])
        n += 1
    pdb_path.write_text("".join(out))
    return n


def _rename_residues_and_strip_atoms(
    pdb_path: Path,
    positions: dict[tuple[str, int], str],
    strip_atoms: set[str],
    verbose: bool = False,
) -> int:
    """For each ``(chain, resseq)`` in ``positions``: rewrite the
    resname column to the mapped value AND drop any ATOM record whose
    atom name is in ``strip_atoms``.

    Used post-second-tleap to restore CYX/CYM names and remove HG
    (CYX/CYM have no HG per AMBER template).
    """
    if not positions:
        return 0
    renamed = 0
    stripped = 0
    out: list[str] = []
    for raw in pdb_path.read_text().splitlines(keepends=True):
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 27:
            out.append(raw)
            continue
        chain = raw[21]
        try:
            resseq = int(raw[22:26].strip())
        except ValueError:
            out.append(raw)
            continue
        key = (chain, resseq)
        if key not in positions:
            out.append(raw)
            continue
        atom_name = raw[12:16].strip()
        if atom_name in strip_atoms:
            stripped += 1
            continue
        new_resname = positions[key]
        out.append(raw[:17] + f"{new_resname:>3s}" + raw[20:])
        renamed += 1
    pdb_path.write_text("".join(out))
    if verbose and (renamed or stripped):
        print(f"  [prep] restored {renamed} CYX/CYM residue name(s), "
              f"stripped {stripped} HG atom(s)")
    return renamed


def _find_terminal_residues(
    pdb_path: Path,
) -> set[tuple[str, int]]:
    """Return the first + last (chain, resseq) tuple of each chain.

    AMBER ff14SB has NRES/CRES templates for standard residues plus
    HID/HIE/HIP, but NOT for LYN / ASH / GLH / CYX / CYM (no NLYN /
    CLYN / NASH / CASH / etc.). Assigning those variants to terminals
    would produce a residue OpenMM's template matcher can't resolve
    ("matches CLYS but missing 1 H"). Callers of
    :func:`assign_amber_variants` use this set to skip terminal
    residues for the four unsupported-terminal variants.
    """
    per_chain: dict[str, list[tuple[int, str]]] = {}
    for raw in pdb_path.read_text().splitlines():
        if raw.startswith("TER") and len(raw) >= 27:
            # TER carries the chain of the last residue; nothing to add
            # (we've already recorded it via ATOM lines above).
            continue
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 27:
            continue
        chain = raw[21]
        try:
            resseq = int(raw[22:26].strip())
        except ValueError:
            continue
        per_chain.setdefault(chain, []).append((resseq, raw[17:20].strip()))

    terminals: set[tuple[str, int]] = set()
    for chain, entries in per_chain.items():
        if not entries:
            continue
        # De-dup while keeping order (first appearance per resseq).
        seen_resseqs: list[int] = []
        seen_set: set[int] = set()
        for resseq, _resname in entries:
            if resseq not in seen_set:
                seen_set.add(resseq)
                seen_resseqs.append(resseq)
        if seen_resseqs:
            terminals.add((chain, seen_resseqs[0]))
            terminals.add((chain, seen_resseqs[-1]))
    return terminals


def assign_amber_variants(
    pdb_path: Path,
    propka_pkas: dict[tuple[str, int], tuple[str, float]] | None,
    ph: float,
    ss_pairs: set[tuple[str, int]] | None = None,
    verbose: bool = False,
) -> dict[tuple[str, int], str]:
    """Rewrite ATOM/HETATM residue names in ``pdb_path`` (in place) to
    AMBER protonation variants based on PROPKA pKa predictions.

    * ``propka_pkas``: dict ``(chain, resseq) → (resname, pka)`` from
      :mod:`dvbfixer.protonate.run_propka` (or equivalent).
    * ``ph``: reference pH (protonated when pKa > pH for ASP/GLU,
      deprotonated when pKa < pH for LYS/CYS).
    * ``ss_pairs``: set of ``(chain, resseq)`` for CYS residues that
      are in a disulfide bond; those become CYX regardless of pKa.

    HIS variants (HID/HIE/HIP) are already set by reduce and pass
    through unchanged.

    Returns a dict of ``(chain, resseq) → new_resname`` for every
    residue that was renamed; useful for the .dat file.
    """
    ss = ss_pairs or set()
    pkas = propka_pkas or {}
    # Terminals have no NLYN/CLYN/NASH/CASH/NGLH/CGLH/etc. templates in
    # ff14SB — skip variant assignment there and let the residue stay
    # standard (LYS/ASP/GLU/CYS).
    terminals = _find_terminal_residues(pdb_path)
    _NO_TERMINAL_VARIANT = {"LYN", "ASH", "GLH", "CYX", "CYM"}

    renames: dict[tuple[str, int], str] = {}
    out: list[str] = []
    n_terminal_skipped = 0
    for raw in pdb_path.read_text().splitlines(keepends=True):
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 27:
            out.append(raw)
            continue
        chain = raw[21]
        try:
            resseq = int(raw[22:26].strip())
        except ValueError:
            out.append(raw)
            continue
        resname = raw[17:20].strip()
        is_ss = (chain, resseq) in ss
        pk_entry = pkas.get((chain, resseq))
        pka = pk_entry[1] if pk_entry else None
        variant = _classify_variant(resname, pka, ph, is_ss)
        if variant is None or variant == resname:
            out.append(raw)
            continue
        if variant in _NO_TERMINAL_VARIANT and (chain, resseq) in terminals:
            # Skip: no ff14SB terminal template for this variant.
            n_terminal_skipped += 1
            out.append(raw)
            continue
        new = raw[:17] + f"{variant:>3s}" + raw[20:]
        out.append(new)
        renames[(chain, resseq)] = variant
    pdb_path.write_text("".join(out))
    if verbose and n_terminal_skipped:
        print(f"  [prep] skipped {n_terminal_skipped} terminal variant "
              f"rename(s) (no NRES/CRES template in ff14SB for "
              f"LYN/ASH/GLH/CYX/CYM)")
    if verbose and renames:
        by_kind: dict[str, int] = {}
        for name in renames.values():
            by_kind[name] = by_kind.get(name, 0) + 1
        parts = ", ".join(f"{k}={v}" for k, v in sorted(by_kind.items()))
        print(f"  [prep] applied {len(renames)} AMBER variant renames "
              f"({parts})")
    return renames


# ---------------------------------------------------------------------------
# End-to-end orchestrator
# ---------------------------------------------------------------------------


def run_prep(
    input_pdb: str | Path,
    output_pdb: str | Path,
    *,
    ph: float = 7.0,
    ff: str = "leaprc.protein.ff14SB",
    extra_leaprc: list[str] | None = None,
    assign_variants: bool = True,
    ss_pairs: set[tuple[str, int]] | None = None,
    verbose: bool = False,
) -> dict[str, object]:
    """Full deterministic prep pipeline. Returns a metadata dict:

    * ``renames``: dict ``(chain, resseq) → new_resname`` from
      :func:`assign_amber_variants`.
    * ``n_residues``: number of residues in the output.

    Raises :class:`TleapError` / :class:`ReduceError` on tool failure.
    """
    input_pdb = Path(input_pdb)
    output_pdb = Path(output_pdb)

    with tempfile.TemporaryDirectory(prefix="dvbfixer_prep_") as tmp:
        tmpdir = Path(tmp)
        step1 = tmpdir / "1_no_h.pdb"
        step2 = tmpdir / "2_tleap.pdb"
        step3 = tmpdir / "3_no_h.pdb"
        step4 = tmpdir / "4_reduce.pdb"

        _strip_hydrogens(input_pdb, step1)
        run_tleap(step1, step2, ff=ff, extra_leaprc=extra_leaprc,
                  verbose=verbose)
        _strip_hydrogens(step2, step3)
        run_reduce(step3, step4, build=True, nuclear=True, verbose=verbose)

        # Reduce may strip chain IDs; restore from input metadata.
        meta = _capture_residue_meta(input_pdb)
        _restore_metadata(step4, output_pdb, meta)

    renames: dict[tuple[str, int], str] = {}
    if assign_variants:
        propka_dict: dict[tuple[str, int], tuple[str, float]] = {}
        try:
            from dvbfixer.protonate import get_pka_results, run_propka
            mc = run_propka(str(output_pdb))
            for r in get_pka_results(mc):
                try:
                    resnum = int(r["resnum"])
                except (TypeError, ValueError):
                    continue
                propka_dict[(r["chain"], resnum)] = (r["restype"], r["pka"])
        except Exception as e:
            if verbose:
                print(f"  [prep] PROPKA skipped ({e}); "
                      f"variant assignment limited to CYX from SS bonds.")
        renames = assign_amber_variants(
            output_pdb, propka_dict, ph, ss_pairs=ss_pairs,
            verbose=verbose,
        )

        # Second tleap pass — with the variant residue names now in place,
        # tleap places the CORRECT H per each variant template: ASH gets
        # HD2, GLH gets HE2, LYN uses HZ2+HZ3 (no HZ1), HID/HIE/HIP get
        # correct sidechain H. Also re-serializes atoms cleanly (Reduce's
        # "serial 0 + new" markers are gone).
        #
        # CYX/CYM exception: tleap requires explicit `bond mol.X.SG mol.Y.SG`
        # to accept CYX. Temporarily rename CYX/CYM back to CYS for the
        # tleap pass, then post-process to restore the name and strip HG
        # (CYX/CYM have no HG per AMBER template).
        cyx_like: dict[tuple[str, int], str] = {
            k: v for k, v in renames.items() if v in ("CYX", "CYM")
        }
        if cyx_like:
            _rename_residues_in_pdb(output_pdb, cyx_like, target="CYS")
        try:
            with tempfile.TemporaryDirectory(
                    prefix="dvbfixer_prep2_") as tmp2:
                tmpdir2 = Path(tmp2)
                step5 = tmpdir2 / "5_no_h.pdb"
                step6 = tmpdir2 / "6_tleap.pdb"
                _strip_hydrogens(output_pdb, step5)
                run_tleap(step5, step6, ff=ff, extra_leaprc=extra_leaprc,
                          verbose=verbose)
                _restore_metadata(step6, output_pdb, meta)
        except TleapError as e:
            if verbose:
                print(f"  [prep] second tleap pass skipped ({e})")
        if cyx_like:
            _rename_residues_and_strip_atoms(
                output_pdb, cyx_like, strip_atoms={"HG"}, verbose=verbose,
            )

    n_res = 0
    seen: set[tuple[str, str, str]] = set()
    for raw in output_pdb.read_text().splitlines():
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 27:
            continue
        key = (raw[21], raw[22:26].strip(), raw[26].strip())
        if key not in seen:
            seen.add(key)
            n_res += 1

    return {"renames": renames, "n_residues": n_res}
