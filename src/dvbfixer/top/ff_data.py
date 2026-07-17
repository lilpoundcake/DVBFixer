"""Data constants for ``dvbfixer top``.

Split out of the flat ``top.py`` in the Phase 2.4 follow-up work.
Pure data — no functions, no side effects on import. Owns the sugar/
lipid PDB↔GROMACS lookup tables (``PDB_TO_GMX``, ``PDB_TO_CARB``,
``PDB_TO_LIPID``, ``CERAMIDE_RTP``), the multi-water ion Lennard-Jones
parameter matrix (``ION_PARAMS``), and the CHARMM glycosidic-linkage
parameter blob (``_GLYCAN_LINKAGE_PARAMS``).

Everything here is either a Python literal or a triple-quoted GROMACS
FF snippet — no imports from anywhere in dvbfixer. Downstream code
does ``from dvbfixer.top.ff_data import ION_PARAMS`` etc.
"""

from __future__ import annotations

PDB_TO_GMX = {
    'HIS': 'HISE',    # default HIS → HIE/HSE
    'CYS': 'CYS',
    'ASP': 'ASP',
    'GLU': 'GLU',
    'LYS': 'LYS',
    # AMBER protonation names
    'HIE': 'HISE', 'HID': 'HISD', 'HIP': 'HISH',
    'ASH': 'ASPH', 'GLH': 'GLUH',
    'CYX': 'CYS2', 'CYM': 'CYS',
    'LYN': 'LYSN',
    # CHARMM protonation names
    'HSE': 'HISE', 'HSD': 'HISD', 'HSP': 'HISH',
    'ASPP': 'ASPH', 'GLUP': 'GLUH',
    'LSN': 'LYSN',
    # Non-canonical → canonical
    'MSE': 'MET',
    # GLYCAM glycosylated protein residues (no NLN/OLS/OLT in AMBER99SB-ILDN
    # or CHARMM36 RTP — map to standard parent; protein-glycan bond detection
    # adds the cross-residue bond to interchain_ss.itp)
    'NLN': 'ASN', 'OLS': 'SER', 'OLT': 'THR',
}

# Standard amino acids that can appear in PDB
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    'TYR', 'VAL',
}

# PDB→CHARMM atom name mapping for sugars
# PDB uses different naming for NAG acetyl group and some other atoms
CARB_ATOM_MAP = {
    'NAG': {
        'N2': 'N', 'C7': 'C', 'O7': 'O', 'C8': 'CT',
        'H81': 'HT1', 'H82': 'HT2', 'H83': 'HT3',
        'HN2': 'HN',
    },
    'NDG': {
        'N2': 'N', 'C7': 'C', 'O7': 'O', 'C8': 'CT',
        'H81': 'HT1', 'H82': 'HT2', 'H83': 'HT3',
        'HN2': 'HN',
    },
    'NGA': {
        'N2': 'N', 'C7': 'C', 'O7': 'O', 'C8': 'CT',
        'H81': 'HT1', 'H82': 'HT2', 'H83': 'HT3',
        'HN2': 'HN',
    },
    'A2G': {
        'N2': 'N', 'C7': 'C', 'O7': 'O', 'C8': 'CT',
        'H81': 'HT1', 'H82': 'HT2', 'H83': 'HT3',
        'HN2': 'HN',
    },
    # Mannose: PDB may use H11/H12 for C1 hydrogens
    'MAN': {'H11': 'H1'},
    'BMA': {'H11': 'H1'},
    # Galactose/Glucose: similar
    'GAL': {'H11': 'H1'},
    'GLC': {'H11': 'H1'},
}

# PDB sugar names -> CHARMM carb.rtp names
PDB_TO_CARB = {
    'NAG': 'BGLCNA',   # N-acetylglucosamine (beta)
    'NDG': 'BGLCNA',   # 2-(acetylamino)-2-deoxy-alpha-D-glucopyranose
    'BMA': 'BMAN',     # beta-mannose
    'MAN': 'AMAN',     # alpha-mannose
    'GAL': 'BGAL',     # beta-galactose
    'GLC': 'BGLC',     # beta-glucose
    'FUC': 'AFUC',     # alpha-fucose
    'FUL': 'BFUC',     # beta-fucose
    'SIA': 'ANE5AC',   # sialic acid (Neu5Ac)
    'NGA': 'BGALNA',   # N-acetylgalactosamine (beta)
    'A2G': 'AGALNA',   # N-acetylgalactosamine (alpha)
    'BGC': 'BGLC',     # beta-glucose
    'XYS': 'BXYL',     # beta-xylose
    'AFU': 'AFUC',     # alpha-L-fucose (alternate PDB code)
    'AMA': 'AMAN',     # alpha-mannose (alternate PDB code)
    'BGA': 'BGAL',     # beta-galactose (alternate PDB code)
    'BGL': 'BGLCNA',   # beta-N-acetylglucosamine (alternate PDB code)
    # CHARMM-GUI native names (already CHARMM RTP names, pass through)
    'BGLC': 'BGLC',
    'BGAL': 'BGAL',
    'AFUC': 'AFUC',
    'AMAN': 'AMAN',
    'BMAN': 'BMAN',
    'BGLCNA': 'BGLCNA',
    'BGALNA': 'BGALNA',
    'AGALNA': 'AGALNA',
    'ANE5AC': 'ANE5AC',
    'ANE5': 'ANE5AC',
    # Non-standard PDB names (from transplant/GLYCAM workflows)
    'AGL': 'AGAL',     # alpha-galactose (_resolve_sugar_rtp auto-detects AGALNA if N-acetyl)  # CHARMM-GUI short name for sialic acid
}

# PDB lipid names -> CHARMM lipid.rtp names (ceramides)
PDB_TO_LIPID = {
    'CER1': 'CER160',   # CHARMM-GUI: CER1 = d18:1/16:0
    'CER2': 'CER2',     # generic ceramide
    'CER3': 'CER3E',    # ceramide variant
}

# CHARMM lipid.rtp ceramide residue names
CERAMIDE_RTP = {
    'CER160', 'CER180', 'CER181', 'CER2', 'CER200',
    'CER220', 'CER240', 'CER241', 'CER3E',
}

# 4-char resnames from CGenFF/solvent that GROMACS writes into PDB cols 17-20
# (standard PDB uses 3-char in cols 18-20, but GROMACS left-justifies from col 17)
_KNOWN_4CHAR_RESNAMES = {
    'ACET', 'ACEH', 'ACEM',  # acetate/acetic acid/acetamide (CGenFF)
    'TIP3', 'SPC', 'SPCE',   # water models
}

# Water residue names (for counting SOL molecules in PDB)
_WATER_RESNAMES = {'SOL', 'HOH', 'WAT', 'TIP3', 'SPC', 'SPCE', 'TIP4', 'TIP5'}

# ---------------------------------------------------------------------------
# Water-model-matched ion Lennard-Jones parameters (AMBER side)
#
# Sources:
#   JC = Joung & Cheatham, J Phys Chem B 112, 9020 (2008) — monovalents for
#        TIP3P, SPC/E, TIP4P-Ew.
#   LM = Li, Roberts, Chakravorty, Merz, JCTC 9, 2733 (2013) — 12-6 divalents
#        (Ca/Mg/Zn) for the same three water models.
#   LSM = Li, Song, Merz, JCTC 16, 4429 (2020) [PMC8173364] — 12-6 divalents
#        for OPC/OPC3/TIP3P-FB/TIP4P-FB.
#   SLM = Sengupta, Li, Wynn, Merz, JCIM 61, 869 (2021) [PMC8173365] — 12-6
#        monovalents for the same OPC family.
#
# σ_nm = R*_Å × 0.17818; ε_kJ = ε_kcal × 4.184.
# Atom-type names match the bundled FF/amber99sb-ildn-lipid21.ff/ions.itp
# moleculetype references: Na, K, Cl, C0 (Ca²⁺), MG, Zn.
#
# Format: {ion_set: {atomtype: (atnum, mass, sigma_nm, eps_kj, charge,
#                                  moleculetype_resname)}}
# ---------------------------------------------------------------------------
ION_PARAMS = {
    'jc-tip3p': {
        # JC TIP3P monovalent: Na 1.369/0.0874, K 1.705/0.19368, Cl 2.513/0.03559
        'Na': (11, 22.990, 0.24393, 0.36586, +1.0, 'NA'),
        'K':  (19, 39.098, 0.30380, 0.81036, +1.0, 'K'),
        'Cl': (17, 35.453, 0.44786, 0.14887, -1.0, 'CL'),
        # LM 12-6 HFE TIP3P divalent: Ca 1.649/0.10593, Mg 1.360/0.01020, Zn 1.271/0.00330
        'C0': (20, 40.078, 0.29381, 0.44317, +2.0, 'CA'),
        'MG': (12, 24.305, 0.24232, 0.04269, +2.0, 'MG'),
        'Zn': (30, 65.380, 0.22647, 0.01382, +2.0, 'ZN'),
    },
    'jc-spce': {
        # JC SPC/E monovalent: Na 1.212/0.35264, K 1.593/0.42971, Cl 2.711/0.01279
        'Na': (11, 22.990, 0.21595, 1.47545, +1.0, 'NA'),
        'K':  (19, 39.098, 0.28384, 1.79789, +1.0, 'K'),
        'Cl': (17, 35.453, 0.48309, 0.05349, -1.0, 'CL'),
        # LM 12-6 HFE SPC/E divalent: Ca 1.635/0.09788, Mg 1.360/0.01020, Zn 1.276/0.00354
        'C0': (20, 40.078, 0.29132, 0.40955, +2.0, 'CA'),
        'MG': (12, 24.305, 0.24232, 0.04269, +2.0, 'MG'),
        'Zn': (30, 65.380, 0.22736, 0.01482, +2.0, 'ZN'),
    },
    'jc-tip4pew': {
        # JC TIP4P-Ew monovalent: Na 1.226/0.16844, K 1.590/0.27947, Cl 2.760/0.01166
        'Na': (11, 22.990, 0.21845, 0.70474, +1.0, 'NA'),
        'K':  (19, 39.098, 0.28332, 1.16928, +1.0, 'K'),
        'Cl': (17, 35.453, 0.49179, 0.04880, -1.0, 'CL'),
        # LM 12-6 HFE TIP4P-Ew divalent: Ca 1.657/0.11069, Mg 1.353/0.00942, Zn 1.252/0.00251
        'C0': (20, 40.078, 0.29524, 0.46312, +2.0, 'CA'),
        'MG': (12, 24.305, 0.24108, 0.03941, +2.0, 'MG'),
        'Zn': (30, 65.380, 0.22308, 0.01051, +2.0, 'ZN'),
    },
    'lm-hfe-opc': {
        # SLM 2021 HFE-OPC monovalent: Na 1.4670/0.02960, K 1.7020/0.13954, Cl 2.3600/0.67879
        'Na': (11, 22.990, 0.26139, 0.12385, +1.0, 'NA'),
        'K':  (19, 39.098, 0.30336, 0.58385, +1.0, 'K'),
        'Cl': (17, 35.453, 0.42050, 2.84010, -1.0, 'CL'),
        # LSM 2020 HFE-OPC divalent: Ca 1.4930/0.03685, Mg 1.2390/0.00206, Zn 1.1510/0.00046
        'C0': (20, 40.078, 0.26602, 0.15419, +2.0, 'CA'),
        'MG': (12, 24.305, 0.22076, 0.00864, +2.0, 'MG'),
        'Zn': (30, 65.380, 0.20509, 0.00192, +2.0, 'ZN'),
    },
    'lm-iod-opc': {
        # SLM 2021 IOD-OPC monovalent: Na 1.4400/0.02322, K 1.7380/0.16500, Cl 2.1500/0.52153
        'Na': (11, 22.990, 0.25658, 0.09715, +1.0, 'NA'),
        'K':  (19, 39.098, 0.30978, 0.69036, +1.0, 'K'),
        'Cl': (17, 35.453, 0.38308, 2.18207, -1.0, 'CL'),
        # LSM 2020 IOD-OPC divalent: Ca 1.5900/0.07447, Mg 1.3730/0.01179, Zn 1.3730/0.01179
        'C0': (20, 40.078, 0.28331, 0.31159, +2.0, 'CA'),
        'MG': (12, 24.305, 0.24464, 0.04935, +2.0, 'MG'),
        'Zn': (30, 65.380, 0.24464, 0.04935, +2.0, 'ZN'),
    },
    'dang-legacy': {
        # Bundled FF/amber99sb-ildn-lipid21.ff/ffnonbonded.itp (Aqvist Na, Dang Cl,
        # Aqvist K + Allnér Mg + Hoops Zn + Bradbrook Ca). Kept for backward
        # compatibility with topologies generated before water-matched ions.
        'Na': (11, 22.990, 0.33284, 0.01159, +1.0, 'NA'),
        'K':  (19, 39.098, 0.47360, 0.00137, +1.0, 'K'),
        'Cl': (17, 35.453, 0.44010, 0.41840, -1.0, 'CL'),
        'C0': (20, 40.078, 0.30524, 1.92376, +2.0, 'CA'),
        'MG': (12, 24.305, 0.14123, 3.74342, +2.0, 'MG'),
        'Zn': (30, 65.380, 0.19600, 0.05230, +2.0, 'ZN'),
    },
}

# Default ion set for each water model (used when --ion-set auto)
_WATER_DEFAULT_ION_SET = {
    'tip3p':   'jc-tip3p',
    'spc':     'jc-spce',     # plain SPC not parametrized by JC — use SPC/E
    'spce':    'jc-spce',
    'tip4p':   'jc-tip4pew',  # plain TIP4P not parametrized by JC — use TIP4P-Ew
    'tip4pew': 'jc-tip4pew',
    'opc':     'lm-hfe-opc',
}

# Water choices that trigger an alias warning (silent substitute)
_WATER_ION_ALIAS = {
    'spc':   'spce',
    'tip4p': 'tip4pew',
}

# Atom-type names in ffnonbonded.itp that ION_PARAMS replaces
_ION_ATOMTYPE_NAMES = {'Na', 'K', 'Cl', 'C0', 'MG', 'Zn'}

# Known atom name differences between AMBER RTP and PDB/IUPAC naming.
# Key = RTP name, value = PDB name.
# Applied per-residue after detecting which convention the PDB uses.
_EXPLICIT_RENAMES = {
    # ILE: AMBER uses CD, PDB/IUPAC uses CD1
    'CD': 'CD1',
    'HD1': 'HD11', 'HD2': 'HD12', 'HD3': 'HD13',
    # C-terminal: AMBER RTP uses OC1/OC2, PDB uses OXT/O
    'OC1': 'OXT', 'OC2': 'O',
}



# ---------------------------------------------------------------------------
# CHARMM glycosidic-linkage parameter extras (by analogy with
# CC321D/CC321C variants). When OC311→OC3C61 at linkages, some atom-type
# combos (CC321-OC3C61, CC3161-OC3C61-CC3162, etc.) aren't in the
# standard FF distribution.
# ---------------------------------------------------------------------------
_GLYCAN_LINKAGE_PARAMS = """\
; ======================================================================
; Extra parameters for glycosidic linkage sites (by analogy)
; ======================================================================

; --- Extra bondtypes ---
[ bondtypes ]
; i       j     func    b0          kb
  CC321   OC3C61     1   0.14150000    301248.00 ; from CC321D OC3C61

; --- Extra angletypes ---
[ angletypes ]
; i       j       k     func    theta0      ktheta      rub         kub
  HCA2    CC321   OC3C61     5   109.500000   376.560000   0.00000000         0.00 ; from HCA2 CC321D OC3C61
  CC3163  CC321   OC3C61     5   111.500000   376.560000   0.00000000         0.00 ; from OC3C61 CC321D CC311C
  CC3161  OC3C61  CC3162     5   109.700000   794.960000   0.00000000         0.00 ; from CC3163 OC3C61 CC3162
  CC3162  CC3161  OC3C61     5   106.000000   376.560000   0.00000000         0.00 ; from CC3161 CC3162 OC3C61
  CC321   OC3C61  CC3162     5   109.700000   794.960000   0.00000000         0.00 ; from CC321D OC3C61 CC321C
  OC3C61  CC3162  OC3C61     5   112.000000   753.120000   0.00000000         0.00 ; from OC301 CC3162 OC3C61

; --- Extra dihedraltypes ---
[ dihedraltypes ]
; i       j       k       l     func    phi0        kphi        mult
; C-C-O-C glycosidic torsions (from CC3161 CC3162 OC3C61 CC3163)
  CC3161  CC3163  CC321   OC3C61     9     0.000000     0.836800     3 ; from par27 X CT1 CT2 X
  CC3161  OC3C61  CC3162  HCA1       9     0.000000     0.836800     3 ; from HCA1 CC3161 CC3162 OC3C61
  CC3161  OC3C61  CC3162  OC3C61     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  CC3162  CC3161  CC3161  OC3C61     9   180.000000     1.297040     3 ; from CC3161 CC3161 CC3162 OC3C61
  CC3162  CC3161  OC3C61  CC3162     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  CC3163  CC321   OC3C61  CC3162     9     0.000000     0.836800     3 ; from CC321 CC3163 OC3C61 CC3162
  CC321   CC3163  CC3161  OC3C61     9     0.000000     0.836800     3 ; from par27 X CT1 CT2 X
  CC321   OC3C61  CC3162  CC3161     9     0.000000     0.836800     3 ; from CC3161 CC3163 OC3C61 CC3162
  CC321   OC3C61  CC3162  HCA1       9     0.000000     0.836800     3 ; from HCA1 CC3161 CC3162 OC3C61
  CC321   OC3C61  CC3162  OC3C61     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  HCA1    CC3161  OC3C61  CC3162     9     0.000000     0.836800     3 ; from HCA1 CC3161 CC3162 OC3C61
  HCA1    CC3162  CC3161  OC3C61     9     0.000000     0.836800     3 ; from HCA1 CC3161 CC3162 OC3C61
  HCA1    CC3163  CC321   OC3C61     9     0.000000     0.836800     3 ; from HCA2 CC321 CC3163 OC3C61
  HCA2    CC321   OC3C61  CC3162     9     0.000000     0.836800     3 ; from par27 X CT2 OC30A X
  OC3C61  CC3161  CC3161  CC3163     9   180.000000     1.297040     3 ; from CC3161 CC3161 CC3162 OC3C61
  OC3C61  CC3162  CC3161  OC3C61     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  OC3C61  CC3162  OC3C61  CC3163     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  OC3C61  CC3163  CC321   OC3C61     9     0.000000     0.836800     3 ; from par27 X CT2 CT1 X
; Additional C-C-O-C and C-O-C-C linkage torsions
  CC3161  CC3161  OC3C61  CC3162     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  CC3161  OC3C61  CC3162  CC3161     9     0.000000     0.836800     3 ; from CC3161 CC3163 OC3C61 CC3162
  CC3163  CC3161  OC3C61  CC3162     9     0.000000     0.836800     3 ; from CC3161 CC3163 OC3C61 CC3162

; ======================================================================
; Extra parameters for sialic acid (C2) glycosidic linkage sites
; ANE5AC C2 (type CC3062) links to parent sugar O (type OC3C61)
; ======================================================================

; --- Extra angletypes (sialic acid linkage) ---
[ angletypes ]
; i       j       k     func    theta0      ktheta      rub         kub
  CC3161  OC3C61  CC3062     5   109.700000   794.960000   0.00000000         0.00 ; from CC3163 OC3C61 CC3062 (sugar C-O-sialic C2)
  OC3C61  CC3062  OC3C61     5   112.000000   753.120000   0.00000000         0.00 ; from OC301 CC3062 OC3C61 (ring O - C2 - linked O)

; --- Extra dihedraltypes (sialic acid linkage) ---
[ dihedraltypes ]
; i       j       k       l     func    phi0        kphi        mult
; Torsions around CC3161-OC3C61-CC3062 sialic acid linkage
  CC3161  CC3161  OC3C61  CC3062     9     0.000000     0.836800     3 ; from CC3161 CC3163 OC3C61 CC3062
  CC3161  OC3C61  CC3062  CC2O2      9     0.000000     0.836800     3 ; from CC3163 OC3C61 CC3062 CC2O2
  CC3161  OC3C61  CC3062  OC3C61     9     0.000000     0.836800     3 ; from CC3163 OC3C61 CC3062 OC3C61
  CC3161  OC3C61  CC3062  CC3261     9     0.000000     0.836800     3 ; from CC3163 OC3C61 CC3062 CC3261
  HCA1    CC3161  OC3C61  CC3062     9     0.000000     1.188256     3 ; from HCA1 CC3163 OC3C61 CC3062
  OC3C61  CC3161  CC3161  OC3C61     9     0.000000     0.836800     3 ; from OC3C61 CC3162 CC3161 OC3C61
  OC3C61  CC3062  OC3C61  CC3163     9     0.000000     0.836800     3 ; from OC3C61 CC3162 OC3C61 CC3163
  OC3C61  CC3161  OC3C61  CC3062     9     0.000000     0.836800     3 ; from OC3C61 CC3161 OC3C61 CC3162
  NC2D1   CC3161  CC3161  OC3C61     9     0.000000     0.836800     3 ; from OC3C61 CC3261 CC3161 NC2D1

"""

