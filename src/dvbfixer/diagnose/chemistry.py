"""Chemistry / bond-geometry checks for ``dvbfixer diagnose``.

Uses hardcoded canonical values (AMBER14 protein defaults + MolProbity
guidance). All checks operate on the loaded OpenMM ``Topology`` and
``positions``; no force-field parametrization is required (which
keeps ``diagnose`` importable in the CI fast lane).
"""

from __future__ import annotations

import math
from typing import Any

from dvbfixer.diagnose.report import Finding, Severity

# Bond-count ceiling per element (matches dvbfixer.pull.MAX_BONDS).
MAX_BONDS: dict[str, int] = {
    "C": 4,
    "N": 4,
    "O": 2,
    "S": 2,
    "H": 1,
    "P": 5,
}

# Canonical bond lengths (nm), atom-name-aware where the same
# element pair takes very different values depending on hybridisation.
# Keyed first by sorted (atom_name_1, atom_name_2) — falls back to
# element pair if the name pair isn't listed.
_CANONICAL_BOND_NM_BY_NAMES: dict[tuple[str, str], float] = {
    # Backbone carbonyl C=O — partial double bond, ~1.23 Å (not 1.43)
    ("C", "O"): 0.123,
    # Backbone peptide C-N (amide) — ~1.33 Å (not 1.47)
    ("C", "N"): 0.133,
    # C-terminal carboxylate
    ("C", "OXT"): 0.125,
    # SER / THR / TYR hydroxyl O-H
    ("HG", "OG"): 0.096,
    ("HG1", "OG1"): 0.096,
    ("HH", "OH"): 0.096,
    # CYS thiol
    ("HG", "SG"): 0.134,
    # Disulfide
    ("SG", "SG"): 0.205,
}

# Element-pair fallback (nm). Only used when the name pair isn't in
# _CANONICAL_BOND_NM_BY_NAMES. Keyed by sorted (element1, element2).
_CANONICAL_BOND_NM: dict[tuple[str, str], float] = {
    ("C", "C"): 0.153,   # sp3 C-C
    ("C", "H"): 0.109,
    ("C", "N"): 0.147,   # sp3 C-N (fallback; peptide C-N handled above)
    ("C", "O"): 0.143,   # sp3 C-O sidechain (backbone C=O handled above)
    ("H", "N"): 0.101,
    ("H", "O"): 0.096,
    ("H", "S"): 0.134,
    ("N", "N"): 0.145,
    ("O", "P"): 0.161,
    ("S", "S"): 0.204,
}

# Fractional tolerance for bond-length outlier detection. AMBER FF
# springs allow ~10% deviation before large forces kick in; > 30%
# indicates a broken input.
_BOND_LEN_WARN_FRAC = 0.10   # WARNING
_BOND_LEN_ERROR_FRAC = 0.30  # ERROR


def _resid_str(res_id: str, icode: str = "") -> str:
    icode = (icode or "").strip()
    return f"{res_id}{icode}" if icode else str(res_id)


def _res_loc(res: Any) -> tuple[str, str]:
    icode = getattr(res, "insertionCode", "").strip()
    return res.chain.id, _resid_str(res.id, icode)


def _bond_key(a: Any, b: Any) -> tuple[str, str]:
    return tuple(sorted([a.element.symbol, b.element.symbol]))  # type: ignore[return-value]


def _name_key(a: Any, b: Any) -> tuple[str, str]:
    return tuple(sorted([a.name, b.name]))  # type: ignore[return-value]


def _canonical_length_nm(a: Any, b: Any) -> float | None:
    """Return the canonical bond length in nm, preferring atom-name-
    specific values (which distinguish e.g. backbone C=O from sidechain
    C-O) over element-pair fallbacks.
    """
    name_specific = _CANONICAL_BOND_NM_BY_NAMES.get(_name_key(a, b))
    if name_specific is not None:
        return name_specific
    return _CANONICAL_BOND_NM.get(_bond_key(a, b))


def _pos(positions: Any, atom: Any) -> tuple[float, float, float]:
    from openmm.unit import nanometer
    p = positions[atom.index].value_in_unit(nanometer)
    return float(p[0]), float(p[1]), float(p[2])


def _dist(positions: Any, a: Any, b: Any) -> float:
    ax, ay, az = _pos(positions, a)
    bx, by, bz = _pos(positions, b)
    return math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)


def check_valences(topology: Any) -> list[Finding]:
    """Every atom's bond count must not exceed the per-element ceiling.

    Excess bonds usually mean spurious CONECT records inferred by
    OpenBabel or an upstream tool.
    """
    counts: dict[int, int] = {}
    for b1, b2 in topology.bonds():
        counts[b1.index] = counts.get(b1.index, 0) + 1
        counts[b2.index] = counts.get(b2.index, 0) + 1

    findings: list[Finding] = []
    for atom in topology.atoms():
        if atom.element is None:
            continue
        sym = atom.element.symbol
        cap = MAX_BONDS.get(sym)
        if cap is None:
            continue
        n = counts.get(atom.index, 0)
        if n > cap:
            chain, resid = _res_loc(atom.residue)
            findings.append(Finding(
                severity=Severity.ERROR,
                category="valence",
                chain=chain,
                resid=resid,
                resname=atom.residue.name,
                atom=atom.name,
                message=f"valence {n} exceeds max {cap} for {sym} "
                        f"(spurious CONECT?)",
                fix_hint="manual: inspect CONECT records; "
                         "or `dvbfixer conect --no-infer-conect` to "
                         "regenerate cleanly",
            ))
    return findings


def check_bond_lengths(topology: Any, positions: Any) -> list[Finding]:
    """For every bond in the topology, check that its length is within
    ±10% of the canonical AMBER value (WARNING; > 30% → ERROR).
    """
    findings: list[Finding] = []
    for b1, b2 in topology.bonds():
        if b1.element is None or b2.element is None:
            continue
        canonical = _canonical_length_nm(b1, b2)
        if canonical is None:
            continue
        d = _dist(positions, b1, b2)
        # Skip virtual bonds (some inter-residue OpenMM refs)
        if d < 0.05:  # < 0.5 Å — coincident atoms, reported elsewhere
            continue
        frac = abs(d - canonical) / canonical
        if frac < _BOND_LEN_WARN_FRAC:
            continue
        severity = Severity.ERROR if frac >= _BOND_LEN_ERROR_FRAC else Severity.WARNING
        chain, resid = _res_loc(b1.residue)
        atom_pair = f"{b1.name}-{b2.name}"
        if b1.residue is not b2.residue:
            atom_pair = f"{b1.residue.name}{b1.residue.id}:{b1.name}-" \
                        f"{b2.residue.name}{b2.residue.id}:{b2.name}"
        findings.append(Finding(
            severity=severity,
            category="bond_length",
            chain=chain,
            resid=resid,
            resname=b1.residue.name,
            atom=atom_pair,
            message=f"bond {b1.element.symbol}-{b2.element.symbol}: "
                    f"{d * 10:.3f} Å (expected {canonical * 10:.3f} Å, "
                    f"deviation {frac * 100:.1f}%)",
            fix_hint="dvbfixer minimize (energy minimisation relaxes bond "
                     "lengths to their FF targets)",
        ))
    return findings


def _dihedral_deg(
    p1: tuple[float, float, float],
    p2: tuple[float, float, float],
    p3: tuple[float, float, float],
    p4: tuple[float, float, float],
) -> float:
    """Standard 4-atom dihedral, returns angle in degrees in (-180, 180]."""
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def cross(a, b):
        return (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    def norm(a):
        return math.sqrt(dot(a, a))

    b1 = sub(p2, p1)
    b2 = sub(p3, p2)
    b3 = sub(p4, p3)
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    m = cross(n1, sub((0.0, 0.0, 0.0), (-b2[0], -b2[1], -b2[2])))
    len_b2 = norm(b2)
    x = dot(n1, n2) * len_b2
    y = dot(m, n2)
    return math.degrees(math.atan2(y, x))


def _find_atom(res: Any, name: str) -> Any | None:
    for a in res.atoms():
        if a.name == name:
            return a
    return None


def check_peptide_omegas(topology: Any, positions: Any) -> list[Finding]:
    """Report cis peptides (|ω| close to 0°) and non-planar amides.

    ω = dihedral Cα(i) – C(i) – N(i+1) – Cα(i+1). Normal trans amide
    ω ≈ 180°. cis (~0°) is rare except at PRO where it occurs in ~5%
    of residues naturally.
    """
    findings: list[Finding] = []
    # Walk consecutive residues on each chain
    for chain in topology.chains():
        residues = list(chain.residues())
        for i in range(len(residues) - 1):
            r1 = residues[i]
            r2 = residues[i + 1]
            ca1 = _find_atom(r1, "CA")
            c1 = _find_atom(r1, "C")
            n2 = _find_atom(r2, "N")
            ca2 = _find_atom(r2, "CA")
            if not (ca1 and c1 and n2 and ca2):
                continue
            omega = _dihedral_deg(
                _pos(positions, ca1),
                _pos(positions, c1),
                _pos(positions, n2),
                _pos(positions, ca2),
            )
            abs_omega = abs(omega)
            is_pro = r2.name in ("PRO", "HYP")
            chain_id, resid = _res_loc(r2)

            # cis peptide: |omega| within 30° of 0°.
            if abs_omega < 30:
                sev = Severity.INFO if is_pro else Severity.WARNING
                findings.append(Finding(
                    severity=sev,
                    category="cis_peptide",
                    chain=chain_id,
                    resid=resid,
                    resname=r2.name,
                    message=f"cis peptide (ω={omega:.1f}°) "
                            f"{'— rare but natural at PRO' if is_pro else '— unusual'}",
                    fix_hint="manual: verify against a reference structure; "
                             "MolProbity flags all cis non-PRO peptides",
                ))
            # non-planar: ω in [30°, 150°] range (neither cis nor trans)
            elif abs_omega < 150:
                findings.append(Finding(
                    severity=Severity.WARNING,
                    category="non_planar_amide",
                    chain=chain_id,
                    resid=resid,
                    resname=r2.name,
                    message=f"non-planar amide (ω={omega:.1f}°); trans is "
                            f"~180°, cis is ~0°",
                    fix_hint="dvbfixer minimize (amide planarity is one of "
                             "the strongest FF constraints)",
                ))
    return findings


def check_ca_chirality(topology: Any, positions: Any) -> list[Finding]:
    """Report D-amino acids via triple-product sign check.

    For a Cα with substituents N, C, CB, H the L / D configuration
    corresponds to the sign of the scalar triple product

        (v_N × v_C) · v_CB

    where v_X = pos(X) - pos(CA). L-amino acids have a POSITIVE triple
    product (right-handed convention); D-amino acids give a negative
    value. GLY has no CB (skipped). Non-standard residues missing any
    of N/CA/C/CB are skipped.
    """
    findings: list[Finding] = []
    for res in topology.residues():
        if res.name in ("GLY", "HOH", "WAT", "TIP3"):
            continue
        n = _find_atom(res, "N")
        ca = _find_atom(res, "CA")
        c = _find_atom(res, "C")
        cb = _find_atom(res, "CB")
        if not (n and ca and c and cb):
            continue
        ca_p = _pos(positions, ca)
        vn = tuple(a - b for a, b in zip(_pos(positions, n), ca_p))
        vc = tuple(a - b for a, b in zip(_pos(positions, c), ca_p))
        vcb = tuple(a - b for a, b in zip(_pos(positions, cb), ca_p))
        cross = (
            vn[1] * vc[2] - vn[2] * vc[1],
            vn[2] * vc[0] - vn[0] * vc[2],
            vn[0] * vc[1] - vn[1] * vc[0],
        )
        triple = cross[0] * vcb[0] + cross[1] * vcb[1] + cross[2] * vcb[2]
        # L-amino acid: triple > 0. D-amino acid: triple < 0.
        # Use a small deadband around zero to avoid false positives on
        # marginal geometry.
        if triple < -1e-4:
            chain, resid = _res_loc(res)
            findings.append(Finding(
                severity=Severity.ERROR,
                category="chirality",
                chain=chain,
                resid=resid,
                resname=res.name,
                atom="CA",
                message=f"D-amino acid Cα stereochemistry "
                        f"(triple product (N×C)·CB = {triple:.4f} nm³; "
                        f"L requires positive)",
                fix_hint="manual: rebuild the residue with correct L stereo",
            ))
    return findings


def run_all(topology: Any, positions: Any) -> list[Finding]:
    """Execute every chemistry check and return the concatenated findings."""
    findings: list[Finding] = []
    findings.extend(check_valences(topology))
    findings.extend(check_bond_lengths(topology, positions))
    findings.extend(check_peptide_omegas(topology, positions))
    findings.extend(check_ca_chirality(topology, positions))
    return findings
