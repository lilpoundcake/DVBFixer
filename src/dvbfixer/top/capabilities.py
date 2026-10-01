"""Declarative capabilities for RTP-based GROMACS topology backends."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol

from dvbfixer.top.ff_data import ION_PARAMS, STANDARD_AA


@dataclass(frozen=True)
class WaterIonPair:
    water_file: str
    ion_file: str | None
    ion_set: str | None
    source: str


class _TopArgs(Protocol):
    water: str
    ion_set: str
    ff_dir: str | None
    merge: bool
    acpype: bool
    pdb: str | None


@dataclass(frozen=True)
class TopologyForceField:
    name: str
    bundle: str
    family: str
    provenance: str
    default_selectable: bool
    water_ion_pairs: dict[str, WaterIonPair]
    required_files: frozenset[str]
    optional_files: frozenset[str] = frozenset()
    extra_rtp_files: tuple[str, ...] = ()
    extra_r2b_files: tuple[str, ...] = ()
    cmap_file: str | None = None
    nbfix_file: str | None = None
    supports_custom_ff_dir: bool = True
    supports_merge: bool = True
    supports_acpype: bool = True
    supports_matched_pdb: bool = True
    fail_closed: bool = False
    accepted_molecule_classes: frozenset[str] = frozenset()
    rejected_molecule_classes: frozenset[str] = frozenset()
    accepted_residues: frozenset[str] = frozenset()
    accepted_input_ions: frozenset[str] = frozenset()
    protonation_variants: frozenset[str] = frozenset()
    caps: frozenset[str] = frozenset()
    strict_terminal_variants: bool = False
    legacy_ion_substitution: bool = False
    bundled_ions: bool = False
    exact_ion_pairing: bool = False
    expected_sha256: dict[str, str] = field(default_factory=dict)

    @property
    def has_arn(self) -> bool:
        return "aminoacids.arn" in self.required_files

    @property
    def has_r2b(self) -> bool:
        return "aminoacids.r2b" in self.required_files

    @property
    def has_tdb(self) -> bool:
        return bool({"aminoacids.n.tdb", "aminoacids.c.tdb"} & self.required_files)

    @property
    def supports_cmap(self) -> bool:
        return self.cmap_file is not None


_CORE = frozenset({
    "aminoacids.rtp", "aminoacids.r2b", "atomtypes.atp",
    "ffbonded.itp", "ffnonbonded.itp", "forcefield.itp",
})

_AMBER_WATERS = {
    name: WaterIonPair(f"{name}.itp", None, ion_set, "DVBfixer reviewed substitution")
    for name, ion_set in {
        "tip3p": "jc-tip3p", "spc": "jc-spce", "spce": "jc-spce",
        "tip4p": "jc-tip4pew", "tip4pew": "jc-tip4pew", "opc": "lm-hfe-opc",
    }.items()
}

_CHARMM_WATERS = {
    name: WaterIonPair(f"{name}.itp", "ions.itp", None, "bundled CHARMM36 files")
    for name in ("tip3p", "spc", "spce")
}

_AMBER19SB_WATERS = {
    name: WaterIonPair(
        f"{name}.itp", f"ions_{name}.itp", f"amber19sb-{name}",
        "GROMACS v2026.3 upstream matched files",
    )
    for name in ("opc", "opc3", "spc", "spce", "tip3p", "tip4pew")
}

_AMBER19SB_HASHES = {
    "aminoacids.c.tdb": "30dd9f50e9d411f1875ee9e62ce8b2872c9c6d5a2c33a26227309e78c453caff",
    "aminoacids.hdb": "9f483656a790a7aabc6b8bdb6182fc38ff12f7eb1f85c8172c4e8009215ff370",
    "aminoacids.n.tdb": "30dd9f50e9d411f1875ee9e62ce8b2872c9c6d5a2c33a26227309e78c453caff",
    "aminoacids.r2b": "a5f2801c60c13f90e17f7129e3c4e25e099f0f9d01c9908a7af4c7ad10081698",
    "aminoacids.rtp": "ece03322c6df812c1467bfdb67240e694230ea28ba675d86cccb8a4823e9a419",
    "atomtypes.atp": "e3cbcf3ee696da0e37196541c565d38a64ca60edafddb8a05c33a17448455866",
    "cmap.itp": "61ceb9bb7ba200cf3b006fcd99407aa2c08dbde1e66e2f68e1bfb8d863a7f015",
    "ffbonded.itp": "b6fa511c5bd3f11a19484a418d10e32972e2c19788d464466926c7d8fbde1034",
    "ffnonbonded.itp": "608949f58a1bb41760c009a10f7c527a5d3e00c5608899836abef45402241da7",
    "forcefield.doc": "fa33def5ddc33648224a85ae2fcc510b7cb9d51f618686636367200e46764f28",
    "forcefield.itp": "4d1dba4fe3d70af435c0cfc2fa9eb620bfd3a88a7fbf78f4687d0e5f3c6463b7",
    "ions_opc.itp": "f34241738b5928ff64394a7d74ed6cf89cdbeeb46b9a1a0757e5d947461a8e22",
    "ions_opc3.itp": "966c0be6a46401bc1e0f6cda6ff572b51a3028a612a452dc21d0fca440244191",
    "ions_spc.itp": "6e4f457a3431fe61d26c013c9aac21b0a3c81fcaebf9079c777ddc78c2fbfeaf",
    "ions_spce.itp": "6e4f457a3431fe61d26c013c9aac21b0a3c81fcaebf9079c777ddc78c2fbfeaf",
    "ions_tip3p.itp": "a40654c00fc42fbf0250e4439be68e7ac3c682c7c14173df5d35fc453e1e9db9",
    "ions_tip4pew.itp": "8dff0830473179f499e99f4f31125a396c25a9d8f7473b9ed57be37efb16507c",
    "opc.itp": "6c757c00bbd28abaad5558011b8b2e84f9fdd9f1636e91f7872777029a5520c1",
    "opc3.itp": "eacaba41ced99ddcc0fafc05edd5c7a83a43d84d09eee07a540983f398e2fbc5",
    "spc.itp": "c15cfab95d07de31b06b8a92a160ffdba1e9ade41e966bd483ff34719852b455",
    "spce.itp": "8116a29837d4da0a9b1e72ccb8b32e6f4b007343af1ccf80f58988d054f8e438",
    "tip3p.itp": "814342d16d1df70885aa99480a68ed69fa4fc41f7b438c172227d918a3c8a985",
    "tip4pew.itp": "63b74979e1f5ae41c9e636901f1dfeafd5b2483d88c3e5b4870fa4988da8faf2",
    "watermodels.dat": "7df5f3bec651d63ca81a51d5d8ac9d7c7bb0f4f29c0f49af89d9e99bdbca9363",
}

_AMBER19SB_RESIDUES = frozenset(STANDARD_AA | {
    "ACE", "NME", "ASH", "GLH", "HID", "HIE", "HIP", "CYX", "CYM", "LYN", "HYP",
})

TOPOLOGY_FORCE_FIELDS = {
    "amber": TopologyForceField(
        name="amber", bundle="amber99sb-ildn-lipid21.ff", family="amber",
        provenance="DVBfixer legacy AMBER99SB-ILDN + Lipid21 bundle",
        default_selectable=True, water_ion_pairs=_AMBER_WATERS,
        required_files=_CORE | {"aminoacids.arn"},
        legacy_ion_substitution=True,
    ),
    "amber19sb": TopologyForceField(
        name="amber19sb", bundle="amber19sb.ff", family="amber19sb",
        provenance="GROMACS v2026.3 (121090014570a53a17ea391bcddae45e5ea05eb4)",
        default_selectable=False, water_ion_pairs=_AMBER19SB_WATERS,
        required_files=frozenset(_AMBER19SB_HASHES),
        optional_files=frozenset({"aminoacids.arn"}), cmap_file="cmap.itp",
        supports_custom_ff_dir=False, supports_acpype=False, fail_closed=True,
        accepted_molecule_classes=frozenset({"protein", "protein-cap", "water", "simple-ion"}),
        rejected_molecule_classes=frozenset({
            "dna", "rna", "carbohydrate", "lipid", "ligand", "ptm", "coordinated-metal",
        }),
        accepted_residues=_AMBER19SB_RESIDUES, caps=frozenset({"ACE", "NME"}),
        accepted_input_ions=frozenset({"NA", "CL", "K"}),
        protonation_variants=frozenset({"ASH", "GLH", "HID", "HIE", "HIP", "CYX", "CYM", "LYN"}),
        strict_terminal_variants=True, bundled_ions=True, exact_ion_pairing=True,
        expected_sha256=_AMBER19SB_HASHES,
    ),
    "charmm": TopologyForceField(
        name="charmm", bundle="charmm36_ljpme-jul2022.ff", family="charmm",
        provenance="bundled CHARMM36 LJ-PME July 2022",
        default_selectable=False, water_ion_pairs=_CHARMM_WATERS,
        required_files=_CORE | {"aminoacids.arn", "aminoacids.n.tdb", "aminoacids.c.tdb"},
        extra_rtp_files=("merged.rtp", "carb.rtp", "lipid.rtp", "na.rtp", "cgenff.rtp",
                         "ethers.rtp", "metals.rtp", "silicates.rtp", "solvent.rtp"),
        extra_r2b_files=("carb.r2b", "lipid.r2b", "na.r2b", "cgenff.r2b", "ethers.r2b",
                         "metals.r2b", "silicates.r2b", "solvent.r2b"),
        cmap_file="cmap.itp", nbfix_file="nbfix.itp", bundled_ions=True,
    ),
}

DEFAULT_TOPOLOGY_FORCE_FIELD = next(
    name for name, descriptor in TOPOLOGY_FORCE_FIELDS.items() if descriptor.default_selectable
)
WATER_CHOICES = tuple(sorted({w for ff in TOPOLOGY_FORCE_FIELDS.values() for w in ff.water_ion_pairs}))
ION_SET_CHOICES = ("auto", *ION_PARAMS, *(
    pair.ion_set
    for pair in _AMBER19SB_WATERS.values()
    if pair.ion_set is not None
))


def resolve_force_field(name: str) -> TopologyForceField:
    return TOPOLOGY_FORCE_FIELDS[name]


def resolve_bundle(descriptor: TopologyForceField, root: Path, custom: str | None = None) -> Path:
    if custom:
        if not descriptor.supports_custom_ff_dir:
            raise ValueError(f"--ff {descriptor.name} requires its verified bundled force field; --ff-dir is not supported")
        bundle = Path(custom)
    else:
        bundle = root / descriptor.bundle
    if not bundle.is_dir():
        raise ValueError(f"Force field directory not found: {bundle}")
    missing = sorted(name for name in descriptor.required_files if not (bundle / name).is_file())
    if missing:
        raise ValueError(f"Force field {descriptor.name} is incomplete; missing: {', '.join(missing)}")
    if descriptor.expected_sha256 and not custom:
        actual_names = {path.name for path in bundle.iterdir() if path.is_file()}
        expected_names = set(descriptor.expected_sha256)
        if actual_names != expected_names:
            raise ValueError(
                f"Force field {descriptor.name} inventory mismatch; expected {len(expected_names)} files, "
                f"found {len(actual_names)}"
            )
        for name, expected in descriptor.expected_sha256.items():
            actual = hashlib.sha256((bundle / name).read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError(f"Force field {descriptor.name} checksum mismatch: {name}")
    return bundle


def validate_options(descriptor: TopologyForceField, args: _TopArgs) -> str | None:
    water = args.water
    ion_set = args.ion_set
    if water not in descriptor.water_ion_pairs:
        allowed = "|".join(descriptor.water_ion_pairs)
        raise ValueError(f"--water {water} is not supported by --ff {descriptor.name}; use {allowed}")
    if args.ff_dir and not descriptor.supports_custom_ff_dir:
        raise ValueError(f"--ff-dir is not supported by --ff {descriptor.name}")
    if args.merge and not descriptor.supports_merge:
        raise ValueError(f"--merge is not supported by --ff {descriptor.name}")
    if args.acpype and not descriptor.supports_acpype:
        raise ValueError("--acpype is a separate AMBER14+GLYCAM route and cannot be combined with --ff amber19sb")
    if args.pdb and not descriptor.supports_matched_pdb:
        raise ValueError(f"--pdb is not supported by --ff {descriptor.name}")

    pair = descriptor.water_ion_pairs[water]
    if descriptor.legacy_ion_substitution:
        selected = pair.ion_set if ion_set == "auto" else ion_set
        if selected not in ION_PARAMS:
            raise ValueError(f"--ion-set {selected} is not supported by --ff {descriptor.name}")
        return selected
    if descriptor.exact_ion_pairing:
        if ion_set == "auto":
            return pair.ion_set
        if ion_set != pair.ion_set:
            raise ValueError(
                f"Unsupported Amber19SB water/ion pair: {water} with {ion_set}; "
                f"use {pair.ion_set} (or --ion-set auto)"
            )
        return ion_set
    return None
