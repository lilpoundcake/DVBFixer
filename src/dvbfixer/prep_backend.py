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
    H + HIS tautomers + ASN/GLN flips) → PROPKA + ``decide_protonation``
    (variant map keyed by ``(chain, resseq, icode)``, filtered by
    PROPKA group type so N+/C- pKas don't overwrite side-chain pKas)
    → overlay Reduce's HID/HIE tautomer choice for neutral HIS →
    overlay SS bonds → CYX → apply variants → output

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


def _patch_variant_hydrogens(
    pdb_path: Path,
    renames: dict[tuple[str, int, str], str],
    verbose: bool = False,
) -> tuple[int, int]:
    """Adjust H set per residue after `assign_amber_variants` renamed
    it — the H that Reduce placed matches the ORIGINAL (standard)
    residue name, not the variant. Text-level rewrite:

    * **CYX / CYM**: drop `HG` (no HG in these templates).
    * **LYN**: drop `HZ1` (ff14SB LYN uses HZ2 + HZ3 only).
    * **ASH**: insert `HD2` at OD2 + 0.97 Å pointing away from OD1
      (protonated carboxyl, sp2 O, in the CG-OD1-OD2 plane).
    * **GLH**: insert `HE2` at OE2 + 0.97 Å opposite OE1.
    * **HIP**: insert whichever of `HD1` / `HE2` Reduce did not place
      (PROPKA-promoted HIP; ``reduce -build`` picks only a single
      tautomer H per residue, so the doubly-protonated form needs the
      second imidazole H patched in). Bond length 1.01 Å, in the
      imidazole ring plane, outward from the ring center.
    * **HID/HIE**: no-op (Reduce already picked and placed the H).

    Returns ``(added, dropped)``.
    """
    import math

    if not renames:
        return 0, 0

    drops_by_resname: dict[str, set[str]] = {
        "CYX": {"HG"},
        "CYM": {"HG"},
        "LYN": {"HZ1"},
    }
    needs_ash = {k for k, v in renames.items() if v == "ASH"}
    needs_glh = {k for k, v in renames.items() if v == "GLH"}
    needs_hip = {k for k, v in renames.items() if v == "HIP"}

    # Read positions once so we can compute geometry for the inserted H.
    coords_per_res: dict[
        tuple[str, int, str], dict[str, tuple[float, float, float]]
    ] = {}
    # Also record which H atoms Reduce already placed on HIP residues so
    # we only add the missing one.
    hip_h_present: dict[tuple[str, int, str], set[str]] = {}
    interesting = needs_ash | needs_glh | needs_hip
    for raw in pdb_path.read_text().splitlines():
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 54:
            continue
        chain = raw[21]
        try:
            resseq = int(raw[22:26].strip())
        except ValueError:
            continue
        icode = raw[26].strip()
        key = (chain, resseq, icode)
        if key not in interesting:
            continue
        atom_name = raw[12:16].strip()
        try:
            x = float(raw[30:38])
            y = float(raw[38:46])
            z = float(raw[46:54])
        except ValueError:
            continue
        coords_per_res.setdefault(key, {})[atom_name] = (x, y, z)
        if key in needs_hip and atom_name in ("HD1", "HE2"):
            hip_h_present.setdefault(key, set()).add(atom_name)

    def _place_bisector_h(center, neigh1, neigh2, bond_len):
        """Place an H at ``center`` at ``bond_len`` Å, in the plane of
        the three atoms, on the side opposite the bisector of the two
        neighbours. Used for sp2 carboxyl-O H and for imidazole ring-N
        H (both are in-plane, one-neighbour-per-side geometry)."""
        v1 = (neigh1[0] - center[0], neigh1[1] - center[1], neigh1[2] - center[2])
        v2 = (neigh2[0] - center[0], neigh2[1] - center[1], neigh2[2] - center[2])
        s = (-(v1[0] + v2[0]), -(v1[1] + v2[1]), -(v1[2] + v2[2]))
        n = math.sqrt(s[0] ** 2 + s[1] ** 2 + s[2] ** 2)
        if n < 1e-6:
            return None
        return (center[0] + bond_len * s[0] / n,
                center[1] + bond_len * s[1] / n,
                center[2] + bond_len * s[2] / n)

    added = 0
    dropped = 0
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
        icode = raw[26].strip()
        variant = renames.get((chain, resseq, icode))
        if variant is None:
            out.append(raw)
            continue
        atom_name = raw[12:16].strip()
        strip = drops_by_resname.get(variant, set())
        if atom_name in strip:
            dropped += 1
            continue
        out.append(raw)

    # Assemble the H atoms we want to insert, keyed by the file-order
    # index of the anchor atom they should follow.
    to_insert: dict[int, list[str]] = {}

    def _queue_insert(key, anchor_atom, h_name, pos):
        for i, ln in enumerate(out):
            if not ln.startswith(("ATOM  ", "HETATM")) or len(ln) < 27:
                continue
            if ln[21] != key[0]:
                continue
            try:
                rs = int(ln[22:26].strip())
            except ValueError:
                continue
            if rs != key[1] or ln[26].strip() != key[2]:
                continue
            if ln[12:16].strip() != anchor_atom:
                continue
            to_insert.setdefault(i + 1, []).append(
                _emit_h_line(ln, h_name, pos)
            )
            return True
        return False

    for key in needs_ash:
        cds = coords_per_res.get(key, {})
        if not {"CG", "OD1", "OD2"}.issubset(cds):
            continue
        pos = _place_bisector_h(cds["OD2"], cds["CG"], cds["OD1"], 0.97)
        if pos is None:
            continue
        if _queue_insert(key, "OD2", "HD2", pos):
            added += 1

    for key in needs_glh:
        cds = coords_per_res.get(key, {})
        if not {"CD", "OE1", "OE2"}.issubset(cds):
            continue
        pos = _place_bisector_h(cds["OE2"], cds["CD"], cds["OE1"], 0.97)
        if pos is None:
            continue
        if _queue_insert(key, "OE2", "HE2", pos):
            added += 1

    for key in needs_hip:
        cds = coords_per_res.get(key, {})
        present = hip_h_present.get(key, set())
        # Place whichever of HD1 / HE2 is missing. Both missing → place
        # both (Reduce couldn't decide; unusual on -build).
        if "HD1" not in present and {"ND1", "CG", "CE1"}.issubset(cds):
            pos = _place_bisector_h(cds["ND1"], cds["CG"], cds["CE1"], 1.01)
            if pos is not None and _queue_insert(key, "ND1", "HD1", pos):
                added += 1
        if "HE2" not in present and {"NE2", "CD2", "CE1"}.issubset(cds):
            pos = _place_bisector_h(cds["NE2"], cds["CD2"], cds["CE1"], 1.01)
            if pos is not None and _queue_insert(key, "NE2", "HE2", pos):
                added += 1

    if to_insert:
        rebuilt: list[str] = []
        for i, ln in enumerate(out):
            rebuilt.append(ln)
            if i + 1 in to_insert:
                rebuilt.extend(to_insert[i + 1])
        out = rebuilt

    pdb_path.write_text("".join(out))
    if verbose and (added or dropped):
        print(f"  [prep] variant H patch: +{added} "
              f"(ASH HD2 / GLH HE2 / HIP HD1|HE2), "
              f"-{dropped} (CYX/CYM HG, LYN HZ1)")
    return added, dropped


def _emit_h_line(template_line: str, atom_name: str,
                 pos: tuple[float, float, float]) -> str:
    """Build a new ATOM record for a hydrogen at `pos`, using
    template_line as the source of chain/resseq/resname/etc."""
    # PDB ATOM record layout (0-indexed):
    #  0-5   record name
    #  6-10  serial
    #  12-15 atom name
    #  17-19 residue name
    #  21    chain
    #  22-25 resseq
    #  26    icode
    #  30-37 x, 38-45 y, 46-53 z
    #  54-59 occupancy
    #  60-65 tempFactor
    #  76-77 element
    chain = template_line[21] if len(template_line) > 21 else " "
    resseq = template_line[22:26] if len(template_line) >= 26 else "    "
    icode = template_line[26] if len(template_line) > 26 else " "
    resname = template_line[17:20] if len(template_line) >= 20 else "   "
    # Right-pad atom name field to 4 chars, left-justified for 3-char names
    # per PDB v3.3 convention.
    if len(atom_name) < 4:
        name_field = f" {atom_name:<3s}"
    else:
        name_field = f"{atom_name:<4s}"
    return (
        f"ATOM  {0:5d} {name_field} {resname:>3s} "
        f"{chain}{resseq}{icode}   "
        f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}          "
        f" H\n"
    )


def _filter_altlocs(
    pdb_path: Path,
    out_path: Path | None = None,
    keep: str = "A",
    verbose: bool = False,
) -> int:
    """Keep only altloc `' '` and `keep` (default 'A'); drop others.

    Also blanks column 17 on the kept lines so downstream tools don't
    trip over the altloc marker. Emits one WARN per affected residue.

    Returns the count of dropped atoms. If ``out_path is None`` the
    file is rewritten in place.
    """
    target = out_path if out_path is not None else pdb_path
    dropped = 0
    warned: set[tuple[str, str, str]] = set()
    out: list[str] = []
    for raw in pdb_path.read_text().splitlines(keepends=True):
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 27:
            out.append(raw)
            continue
        altloc = raw[16]
        if altloc not in (" ", keep):
            resname = raw[17:20].strip()
            chain = raw[21]
            resseq = raw[22:26].strip()
            k = (chain, resseq, resname)
            if k not in warned:
                warned.add(k)
                if verbose:
                    print(f"  [altloc] {chain}/{resname}{resseq}: keeping "
                          f"'{keep}', dropping '{altloc}' (and any "
                          f"further alt-locs)")
            dropped += 1
            continue
        # Blank the altloc column on kept lines.
        if altloc != " ":
            raw = raw[:16] + " " + raw[17:]
        out.append(raw)
    target.write_text("".join(out))
    if verbose and dropped and not warned:
        print(f"  [altloc] dropped {dropped} atom(s) across {len(warned)} "
              f"residue(s)")
    return dropped


def _rename_all_his_variants_to_his(pdb_path: Path) -> int:
    """Rewrite HID/HIE/HIP → HIS in every ATOM/HETATM line so Reduce's
    ``-build`` tautomer decision fires. Reduce treats HID/HIE/HIP as
    fixed and won't touch them; it only performs the network optimisation
    on residues labelled ``HIS``. Called between tleap (which defaults
    HIS→HIE) and reduce so Reduce sees a clean HIS input.
    """
    n = 0
    out: list[str] = []
    for raw in pdb_path.read_text().splitlines(keepends=True):
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 20:
            out.append(raw)
            continue
        resname = raw[17:20].strip()
        if resname in ("HID", "HIE", "HIP"):
            out.append(raw[:17] + "HIS" + raw[20:])
            n += 1
        else:
            out.append(raw)
    pdb_path.write_text("".join(out))
    return n


def _infer_his_tautomers_from_atoms(
    pdb_path: Path,
) -> dict[tuple[str, int, str], str]:
    """Return ``{(chain, resseq, icode): HID/HIE/HIP}`` for each
    HIS/HID/HIE/HIP residue based on which imidazole H atoms are present.

    ``reduce -build -nuclear`` picks the tautomer per residue but leaves
    the residue name as ``HIS``. We inspect the placed H atoms:

    * HD1 present, HE2 absent  → ``HID`` (proton on δ-nitrogen)
    * HE2 present, HD1 absent  → ``HIE`` (proton on ε-nitrogen)
    * Both present             → ``HIP`` (doubly protonated, +1 charge —
      rare on ``-build`` alone; PROPKA is the primary HIP signal)
    * Neither                  → key omitted (deprotonated HIS is
      unusual and downstream should treat the residue as HIS)
    """
    hd1_seen: set[tuple[str, int, str]] = set()
    he2_seen: set[tuple[str, int, str]] = set()
    his_residues: set[tuple[str, int, str]] = set()
    for raw in pdb_path.read_text().splitlines():
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 27:
            continue
        resname = raw[17:20].strip()
        if resname not in ("HIS", "HID", "HIE", "HIP"):
            continue
        chain = raw[21]
        try:
            resseq = int(raw[22:26].strip())
        except ValueError:
            continue
        icode = raw[26].strip()
        key = (chain, resseq, icode)
        his_residues.add(key)
        name = raw[12:16].strip()
        if name == "HD1":
            hd1_seen.add(key)
        elif name == "HE2":
            he2_seen.add(key)

    result: dict[tuple[str, int, str], str] = {}
    for key in his_residues:
        has_hd1 = key in hd1_seen
        has_he2 = key in he2_seen
        if has_hd1 and has_he2:
            result[key] = "HIP"
        elif has_hd1:
            result[key] = "HID"
        elif has_he2:
            result[key] = "HIE"
        # else: deprotonated HIS — omit from map; caller keeps "HIS".
    return result


def _detect_ss_bonds_from_distance(
    pdb_path: Path, cutoff_a: float = 2.5,
) -> set[tuple[str, int]]:
    """Fallback SS-bond detector when CONECT records are missing.

    Returns a set of ``(chain, resseq)`` for every CYS SG within
    ``cutoff_a`` Å of another CYS SG. Distance is Angstroms (PDB units).
    """
    import math

    sg_positions: list[tuple[str, int, float, float, float]] = []
    for raw in pdb_path.read_text().splitlines():
        if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 54:
            continue
        resname = raw[17:20].strip()
        if resname not in ("CYS", "CYX", "CYM"):
            continue
        if raw[12:16].strip() != "SG":
            continue
        chain = raw[21]
        try:
            resseq = int(raw[22:26].strip())
            x = float(raw[30:38])
            y = float(raw[38:46])
            z = float(raw[46:54])
        except ValueError:
            continue
        sg_positions.append((chain, resseq, x, y, z))

    result: set[tuple[str, int]] = set()
    for i, (c1, r1, x1, y1, z1) in enumerate(sg_positions):
        for c2, r2, x2, y2, z2 in sg_positions[i + 1:]:
            d = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
            if d <= cutoff_a:
                result.add((c1, r1))
                result.add((c2, r2))
    return result


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
    variant_map: dict[tuple[str, int, str], str] | None,
    verbose: bool = False,
) -> dict[tuple[str, int, str], str]:
    """Rewrite ATOM/HETATM residue names in ``pdb_path`` (in place) to
    AMBER protonation variants using a pre-decided variant map.

    * ``variant_map``: dict ``(chain, resseq, icode) → new_resname``
      produced by the caller from PROPKA (``decide_protonation``) plus
      Reduce (HIS tautomer inference) plus SS-bond overlay.

    The map is applied verbatim except for terminal residues whose
    target variant is one of ``LYN / ASH / GLH / CYM`` — ff14SB has no
    NRES/CRES template for those four (no RESP charges were ever
    computed for terminal deprotonated/protonated variants of them),
    so we drop the rename for terminal residues to avoid template-match
    failures in downstream OpenMM. All other variants (HID/HIE/HIP,
    CYX) have NRES/CRES coverage and pass through at termini.

    Returns the *applied* variant map ``(chain, resseq, icode) → variant``
    (with terminal-skipped entries removed).
    """
    if not variant_map:
        return {}
    terminals = _find_terminal_residues(pdb_path)  # 2-tuple (chain, resseq)
    _NO_TERMINAL_VARIANT = {"LYN", "ASH", "GLH", "CYM"}

    # Prune terminal-unsupported variants up front so the file rewrite
    # loop stays a plain dict lookup.
    applied: dict[tuple[str, int, str], str] = {}
    terminal_skipped: set[tuple[str, int, str, str]] = set()
    for (chain, resseq, icode), variant in variant_map.items():
        if (variant in _NO_TERMINAL_VARIANT
                and (chain, resseq) in terminals):
            terminal_skipped.add((chain, resseq, icode, variant))
            continue
        applied[(chain, resseq, icode)] = variant

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
        icode = raw[26].strip()
        target = applied.get((chain, resseq, icode))
        if target is None:
            out.append(raw)
            continue
        resname = raw[17:20].strip()
        if target == resname:
            out.append(raw)
            continue
        out.append(raw[:17] + f"{target:>3s}" + raw[20:])
    pdb_path.write_text("".join(out))
    if terminal_skipped:
        n = len(terminal_skipped)
        detail = ", ".join(f"{c}/{r}{i}→{v}"
                           for (c, r, i, v) in sorted(terminal_skipped)[:5])
        more = "" if n <= 5 else f", ... and {n - 5} more"
        print(f"  [prep] skipped {n} terminal residue variant rename(s) "
              f"({detail}{more}) — ff14SB has no NRES/CRES template for "
              f"LYN/ASH/GLH/CYM.")
    if verbose and applied:
        by_kind: dict[str, int] = {}
        for name in applied.values():
            by_kind[name] = by_kind.get(name, 0) + 1
        parts = ", ".join(f"{k}={v}" for k, v in sorted(by_kind.items()))
        print(f"  [prep] applied {len(applied)} AMBER variant renames "
              f"({parts})")
    return applied


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

    * ``renames``: dict ``(chain, resseq, icode) → new_resname`` from
      :func:`assign_amber_variants` (icode is empty string when the
      input has no insertion code — typical for non-antibody inputs).
    * ``n_residues``: number of residues in the output.

    Raises :class:`TleapError` / :class:`ReduceError` on tool failure.
    """
    input_pdb = Path(input_pdb)
    output_pdb = Path(output_pdb)

    with tempfile.TemporaryDirectory(prefix="dvbfixer_prep_") as tmp:
        tmpdir = Path(tmp)
        step0 = tmpdir / "0_altloc_clean.pdb"
        step1 = tmpdir / "1_no_h.pdb"
        step2 = tmpdir / "2_tleap.pdb"
        step3 = tmpdir / "3_no_h.pdb"
        step4 = tmpdir / "4_reduce.pdb"

        # Drop non-A altloc atoms up front; tleap otherwise treats each
        # altloc as a distinct atom and downstream OpenMM warns about
        # duplicates (or crashes).
        _filter_altlocs(input_pdb, step0, keep="A", verbose=verbose)
        _strip_hydrogens(step0, step1)
        run_tleap(step1, step2, ff=ff, extra_leaprc=extra_leaprc,
                  verbose=verbose)
        _strip_hydrogens(step2, step3)
        # tleap converts HIS → HIE (its default); Reduce treats HIE as
        # fixed and won't decide the tautomer per residue. Rename all
        # HIS-variants back to HIS so Reduce actually runs its H-bond
        # network optimisation and picks the tautomer for each residue.
        _rename_all_his_variants_to_his(step3)
        run_reduce(step3, step4, build=True, nuclear=True, verbose=verbose)

        # Reduce may strip chain IDs; restore from altloc-filtered input.
        meta = _capture_residue_meta(step0)
        _restore_metadata(step4, output_pdb, meta)

    renames: dict[tuple[str, int, str], str] = {}
    if assign_variants:
        from dvbfixer.protonate import (
            decide_protonation,
            get_pka_results,
            run_propka,
        )

        # 1. PROPKA — one call, then decide_protonation filters by group
        #    type (ASP/GLU/LYS/CYS/HIS side chains only, N+/C- skipped)
        #    and keys on (chain, resnum, icode). This is the fix for the
        #    old bug where propka_dict[(chain, resnum)] silently
        #    overwrote a side-chain pKa with the terminal N+/C- pKa.
        pka_results: list[dict[str, object]] = []
        n_pka = 0
        try:
            mc = run_propka(str(output_pdb))
            pka_results = get_pka_results(mc)
            n_pka = len(pka_results)
        except Exception as e:
            print(f"  [prep] PROPKA skipped ({e}); variant assignment "
                  f"limited to CYX from SS bonds + HIS tautomers from "
                  f"Reduce.")

        # Note: decide_protonation may over-classify HIS→HIP when a
        # neighbouring positive charge shifts the pKa; that's the intended
        # PROPKA signal. For neutral HIS (his_default=HIE) we overlay
        # Reduce's per-residue HID/HIE tautomer choice below.
        propka_renames = decide_protonation(
            pka_results, ph, his_default="HIE", cys_ss_pka=99.99,
        )

        # 2. Reduce -build picks the HIS tautomer per H-bond network but
        #    leaves the residue name as "HIS"; read the placed H atoms.
        reduce_his = _infer_his_tautomers_from_atoms(output_pdb)

        variant_map: dict[tuple[str, int, str], str] = {}
        for key, variant in propka_renames.items():
            if variant == "HIP":
                # PROPKA-driven HIP wins: reduce -build almost never
                # places both HD1 and HE2 on its own.
                variant_map[key] = "HIP"
            elif variant in ("HIE", "HID"):
                tautomer = reduce_his.get(key)
                variant_map[key] = tautomer if tautomer in ("HID", "HIE") else variant
            else:
                variant_map[key] = variant
        # HIS residues that PROPKA didn't rename (e.g. PROPKA failed or
        # didn't emit a pKa) still get Reduce's tautomer choice.
        for key, tautomer in reduce_his.items():
            if key in variant_map:
                continue
            if tautomer in ("HID", "HIE", "HIP"):
                variant_map[key] = tautomer

        # 3. SS-bond overlay — every CYS in an SS pair becomes CYX
        #    regardless of PROPKA. Matches on (chain, resseq) because
        #    detect_ss_bonds returns 2-tuples and SS-bonded CYS with
        #    insertion codes is essentially never seen.
        effective_ss = set(ss_pairs) if ss_pairs else set()
        if not effective_ss:
            distance_ss = _detect_ss_bonds_from_distance(output_pdb)
            if distance_ss:
                effective_ss = distance_ss
                if verbose:
                    print(f"  [prep] detected {len(distance_ss)} SS-bonded "
                          f"CYS via distance (input had no CONECT)")
        if effective_ss:
            seen_cys: set[tuple[str, int, str]] = set()
            for raw in output_pdb.read_text().splitlines():
                if not raw.startswith(("ATOM  ", "HETATM")) or len(raw) < 27:
                    continue
                if raw[17:20].strip() not in ("CYS", "CYX", "CYM"):
                    continue
                chain = raw[21]
                try:
                    resseq = int(raw[22:26].strip())
                except ValueError:
                    continue
                icode = raw[26].strip()
                key = (chain, resseq, icode)
                if key in seen_cys:
                    continue
                seen_cys.add(key)
                if (chain, resseq) in effective_ss:
                    variant_map[key] = "CYX"

        # 4. Apply the map (with terminal-skip for LYN/ASH/GLH/CYM).
        renames = assign_amber_variants(
            output_pdb, variant_map, verbose=verbose,
        )

        # 5. Always-on PROPKA activity summary.
        by_kind: dict[str, int] = {}
        for name in renames.values():
            by_kind[name] = by_kind.get(name, 0) + 1
        by_kind_str = (", ".join(f"{k}={v}" for k, v in sorted(by_kind.items()))
                       if by_kind else "none")
        print(f"  [prep] PROPKA + Reduce: {n_pka} titratable "
              f"pKas scanned, {len(renames)} residues renamed "
              f"({by_kind_str}); {len(effective_ss)} SS-bonded CYS.")

        # 6. Patch H set per variant template (avoids a second tleap
        #    pass that duplicated C-terminal O/OXT atoms in 0.7.1).
        _patch_variant_hydrogens(output_pdb, renames, verbose=verbose)

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
