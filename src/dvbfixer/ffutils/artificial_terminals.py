"""Normalize protonation variants at FASTA-truncated chain boundaries."""

from __future__ import annotations

from pathlib import Path

from dvbfixer.ffutils.dat import AddedAtom, DatRecord, ResidueSummary
from dvbfixer.model.pipeline import get_atom_sequence, parse_fasta
from dvbfixer.sequence_alignment import align_observed_to_reference

# These variants are meaningful internally, but the standard protein force
# fields do not define their combination with an N/C terminal patch.  At an
# artificial boundary created by ``model --no-terminal`` they are also
# chemically wrong: the coordinate endpoint is not the biological terminus.
_PARENT_AND_VARIANT_H = {
    "ASH": ("ASP", {"HD2"}),
    "GLH": ("GLU", {"HE2"}),
    "LYN": ("LYS", set()),
    "CYM": ("CYS", set()),
}


def normalize_fasta_truncated_terminal_variants(
    pdb_path: str | Path,
    fasta_path: str | Path | None = None,
    dat_path: str | Path | None = None,
    *,
    force_field: str = "amber",
    verbose: bool = False,
) -> list[tuple[str, str, str, str]]:
    """Revert unsupported variants at uncapped N/C termini.

    With a FASTA, only boundaries lying strictly inside the complete sequence
    are considered.  Without one, all physical chain ends are considered.
    ACE/NME-capped ends are always left unchanged.
    """
    if force_field.lower() != "amber":
        return []

    pdb_path = Path(pdb_path)
    lines = pdb_path.read_text().splitlines(keepends=True)
    references = parse_fasta(fasta_path) if fasta_path else None

    residues: dict[str, list[tuple[str, str, str]]] = {}
    all_residues: dict[str, list[tuple[str, str, str]]] = {}
    seen: set[tuple[str, str, str]] = set()
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            continue
        key = (line[21], line[22:26].strip(), line[26].strip())
        if key in seen:
            continue
        seen.add(key)
        entry = (key[1], key[2], line[17:20].strip())
        all_residues.setdefault(key[0], []).append(entry)
        from dvbfixer.model.cli import AA3TO1
        if entry[2] in AA3TO1:
            residues.setdefault(key[0], []).append(entry)

    replacements: dict[tuple[str, str, str], tuple[str, str, set[str]]] = {}
    chains = references if references is not None else residues
    for chain in chains:
        reference = references[chain] if references is not None else None
        chain_residues = residues.get(chain)
        if not chain_residues:
            continue
        names = [entry[2] for entry in all_residues.get(chain, [])]
        candidates = []
        if references is not None:
            if not isinstance(reference, str):
                continue
            observed = get_atom_sequence(lines, chain)
            alignment = align_observed_to_reference(observed, reference)
            if alignment is None:
                continue
            if alignment.start > 0 and "ACE" not in names:
                candidates.append(chain_residues[0])
            if alignment.end < len(reference) and "NME" not in names:
                candidates.append(chain_residues[-1])
        else:
            if "ACE" not in names:
                candidates.append(chain_residues[0])
            if "NME" not in names:
                candidates.append(chain_residues[-1])
        for resid, icode, variant in candidates:
            normalized = _PARENT_AND_VARIANT_H.get(variant)
            if normalized:
                parent, hydrogens = normalized
                replacements[(chain, resid, icode)] = (variant, parent, hydrogens)

    if not replacements:
        return []

    removed_serials: set[int] = set()
    output: list[str] = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27:
            key = (line[21], line[22:26].strip(), line[26].strip())
            replacement = replacements.get(key)
            if replacement:
                _variant, parent, variant_h = replacement
                if line[12:16].strip() in variant_h:
                    try:
                        removed_serials.add(int(line[6:11]))
                    except ValueError:
                        pass
                    continue
                line = f"{line[:17]}{parent:>3s}{line[20:]}"
        elif line.startswith("TER   ") and len(line) >= 27:
            key = (line[21], line[22:26].strip(), line[26].strip())
            replacement = replacements.get(key)
            if replacement:
                line = f"{line[:17]}{replacement[1]:>3s}{line[20:]}"
        output.append(line)

    if removed_serials:
        output = [
            line for line in output
            if not (line.startswith("CONECT") and any(
                field.strip().isdigit() and int(field) in removed_serials
                for field in (line[i:i + 5] for i in range(6, min(len(line), 31), 5))
            ))
        ]
    pdb_path.write_text("".join(output))

    sidecar = Path(dat_path) if dat_path else pdb_path.with_suffix(".dat")
    if sidecar.exists():
        record = DatRecord.load(sidecar)
        overrides = dict(record.variant_overrides or {})
        for (chain, resid, icode), (_old, _parent, _hydrogens) in replacements.items():
            overrides.pop(f"{chain}:{resid}:{icode}", None)
        record.variant_overrides = overrides or None
        cleaned: list[AddedAtom] = []
        for atom in record.added_atoms:
            key = (
                atom["chain"], atom["resid"],
                str(atom.get("icode", "")).strip(),
            )
            replacement = replacements.get(key)
            if replacement:
                _old, parent, hydrogens = replacement
                if atom["atom"] in hydrogens:
                    continue
                updated_atom = atom.copy()
                updated_atom["resname"] = parent
                atom = updated_atom
            cleaned.append(atom)
        record.added_atoms = cleaned
        summary: dict[str, ResidueSummary] = {}
        for atom in cleaned:
            summary_key = f'{atom["chain"]}/{atom["resname"]}{atom["resid"]}'
            bucket = summary.setdefault(summary_key, {"heavy": 0, "hydrogen": 0})
            bucket["hydrogen" if atom.get("element") == "H" else "heavy"] += 1
        record.residue_summary = summary
        record.save(sidecar, verbose=False)

    result = [(chain, resid, old, parent) for (chain, resid, _), (old, parent, _) in replacements.items()]
    if verbose:
        for chain, resid, old, parent in result:
            print(f"  Artificial FASTA-truncated terminus {chain}:{resid}: {old} -> {parent}")
    return result
