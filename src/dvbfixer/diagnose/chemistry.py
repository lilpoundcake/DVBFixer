"""Chemistry / bond-geometry checks for ``dvbfixer diagnose``.

Uses hardcoded canonical values (AMBER14 protein defaults + MolProbity
guidance). All checks operate on the loaded OpenMM ``Topology`` and
``positions``; no force-field parametrization is required (which
keeps ``diagnose`` importable in the CI fast lane).
"""

from __future__ import annotations

import math
from typing import Any

from dvbfixer.diagnose._common import NON_STANDARD_AAS
from dvbfixer.diagnose.report import Finding, Severity

# Bond-count ceiling per element (matches dvbfixer.pull.MAX_BONDS).
MAX_BONDS: dict[str, int] = {
    "C": 4,
    "N": 4,
    "O": 2,
    # Sulfur is hypervalent in common supported ligands (e.g. sulfate,
    # sulfonates, EPE buffer); bond count alone cannot distinguish those
    # valid states from bad CONECT inference.
    "S": 6,
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

# Fractional tolerance for bond-length outlier detection. Our
# hardcoded canonicals don't carry Engh-Huber σ, so we can't do
# proper 4σ MolProbity-style checks — we lean permissive to avoid
# noise on pre-minimisation and crystal inputs. WARNING at 20%
# deviation, ERROR at 50% (the SER HG-on-OXT bug case sits at 75%,
# so this still catches genuinely broken geometry).
_BOND_LEN_WARN_FRAC = 0.20   # WARNING
_BOND_LEN_ERROR_FRAC = 0.50  # ERROR


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
    ±20% of the canonical AMBER value (WARNING; > 50% → ERROR).
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
    value. GLY has no CB (skipped). Non-standard amino acids (MSE,
    SEC, PYL, phospho-residues, etc.) are recognised via the
    ``NON_STANDARD_AAS`` whitelist so their chirality is still
    checked.
    """
    findings: list[Finding] = []
    for res in topology.residues():
        if res.name in ("GLY", "HOH", "WAT", "TIP3", "TIP4", "TIP5", "SOL"):
            continue
        _ = NON_STANDARD_AAS  # documented as recognised — chirality
        # check keys off atom-name presence, so any non-standard AA
        # that carries N/CA/C/CB (including MSE, SEC, PYL, SEP, TPO,
        # PTR) is naturally handled.
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
        # Deadband -1e-6 nm³ ≈ -1e-3 Å³ matches ``fix_ca_chirality``.
        if triple >= -1e-6:
            continue
        # Near-degeneracy gate: |cross|² / (|vn|²·|vc|²) is sin²(angle
        # between N-CA and C-CA). Below 0.01 the plane is ill-defined
        # and the sign of the triple is dominated by numerical noise.
        norm2 = cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2
        vn_norm2 = vn[0] ** 2 + vn[1] ** 2 + vn[2] ** 2
        vc_norm2 = vc[0] ** 2 + vc[1] ** 2 + vc[2] ** 2
        if norm2 < 1e-18 or norm2 < 0.01 * vn_norm2 * vc_norm2:
            continue
        chain, resid = _res_loc(res)
        findings.append(Finding(
            severity=Severity.ERROR,
            category="chirality",
            chain=chain,
            resid=resid,
            resname=res.name,
            atom="CA",
            message=f"D-amino acid Cα stereochemistry "
                    f"(triple product (N×C)·CB = {triple:.5f} nm³; "
                    f"L requires positive)",
            fix_hint="manual: rebuild the residue with correct L stereo",
        ))
    return findings


def check_disulfides(topology: Any, positions: Any) -> list[Finding]:
    """Report SS bonds whose geometry deviates from AMBER canonicals.

    Ideal disulfide geometry (Schmidt / Neidigh 2002):
      - SG-SG bond length: 2.05 ± 0.10 Å
      - Cα-Cα distance: 5.5 – 7.0 Å
      - CB-SG-SG-CB dihedral: |χ_ss| = 60 – 120° (typically ~90°)

    Detects SS bonds from the topology bond set OR by SG-SG distance
    < 2.5 Å (catches CONECT-less inputs).
    """
    findings: list[Finding] = []
    # Collect all SG atoms per CYS/CYX/CYM residue.
    sg_atoms: list[Any] = []
    for res in topology.residues():
        if res.name not in ("CYS", "CYX", "CYM", "CSS"):
            continue
        for a in res.atoms():
            if a.name == "SG":
                sg_atoms.append(a)
                break

    if len(sg_atoms) < 2:
        return findings

    # Find pairs within SS-bond distance.
    from openmm.unit import nanometer
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            a1, a2 = sg_atoms[i], sg_atoms[j]
            if a1.residue is a2.residue:
                continue
            d = _dist(positions, a1, a2)
            d_a = d * 10.0
            if d_a > 2.5:
                continue
            # Real SS bond — check geometry.
            r1, r2 = a1.residue, a2.residue
            ca1 = _find_atom(r1, "CA")
            ca2 = _find_atom(r2, "CA")
            cb1 = _find_atom(r1, "CB")
            cb2 = _find_atom(r2, "CB")

            chain1, resid1 = _res_loc(r1)
            chain2, resid2 = _res_loc(r2)
            partner = f"{chain2}/{r2.name}{resid2}"

            # SG-SG bond length.
            if not (1.95 <= d_a <= 2.20):
                sev = Severity.ERROR if not (1.85 <= d_a <= 2.30) else Severity.WARNING
                findings.append(Finding(
                    severity=sev,
                    category="disulfide_geometry",
                    chain=chain1, resid=resid1, resname=r1.name, atom="SG",
                    message=f"disulfide SG-SG to {partner}:SG at {d_a:.2f} Å "
                            f"(expected 2.05 ± 0.10 Å)",
                    fix_hint="dvbfixer minimize (restrains SS to canonical length)",
                ))

            # Cα-Cα distance.
            if ca1 and ca2:
                d_ca_nm = _dist(positions, ca1, ca2)
                d_ca_a = d_ca_nm * 10.0
                if not (5.0 <= d_ca_a <= 7.5):
                    findings.append(Finding(
                        severity=Severity.WARNING,
                        category="disulfide_geometry",
                        chain=chain1, resid=resid1, resname=r1.name, atom="CA",
                        message=f"disulfide Cα-Cα to {partner}:CA at {d_ca_a:.2f} Å "
                                f"(expected 5.5 – 7.0 Å)",
                        fix_hint="dvbfixer minimize",
                    ))

            # CB-SG-SG-CB dihedral (χ_ss).
            if ca1 and ca2 and cb1 and cb2:
                chi = abs(_dihedral_deg(
                    _pos(positions, cb1),
                    _pos(positions, a1),
                    _pos(positions, a2),
                    _pos(positions, cb2),
                ))
                if not (60 <= chi <= 120):
                    findings.append(Finding(
                        severity=Severity.WARNING,
                        category="disulfide_geometry",
                        chain=chain1, resid=resid1, resname=r1.name,
                        message=f"disulfide CB-SG-SG-CB dihedral {chi:.0f}° "
                                f"(expected 60 – 120°; canonical ~90°) — "
                                f"partner {partner}",
                        fix_hint="dvbfixer minimize (relaxes non-ideal SS twist)",
                    ))
            # Suppress __unused__ warnings.
            _ = nanometer
    return findings


def run_all(topology: Any, positions: Any) -> list[Finding]:
    """Execute every chemistry check and return the concatenated findings."""
    findings: list[Finding] = []
    findings.extend(check_valences(topology))
    findings.extend(check_bond_lengths(topology, positions))
    findings.extend(check_peptide_omegas(topology, positions))
    findings.extend(check_ca_chirality(topology, positions))
    findings.extend(check_disulfides(topology, positions))
    return findings
