"""Shared data structures for `dvbfixer.top` — parsed-PDB representation
and GROMACS topology accumulation. Split out of `top/pipeline.py` so both
it and `top/topology_builder.py` can depend on these without a circular
import between the two.
"""
from dataclasses import dataclass, field


@dataclass
class PDBResidue:
    chain_id: str
    resname: str
    resseq: int
    icode: str
    atoms: list = field(default_factory=list)  # [(name, x, y, z)]


@dataclass
class PDBChain:
    chain_id: str
    residues: list = field(default_factory=list)  # [PDBResidue]


@dataclass
class AtomEntry:
    index: int          # 1-based
    atom_type: str
    resnr: int          # 1-based residue number in chain
    resname: str
    atomname: str
    cgnr: int
    charge: float
    mass: float
    x: float = 0.0     # coordinates for PDB output
    y: float = 0.0
    z: float = 0.0
    chain_id: str = ' '
    orig_resseq: int = 0
    orig_resname: str = ''


@dataclass
class ChainTopology:
    name: str
    nrexcl: int
    atoms: list = field(default_factory=list)       # [AtomEntry]
    bonds: list = field(default_factory=list)        # [(i, j)]
    pairs: list = field(default_factory=list)        # [(i, j)]
    angles: list = field(default_factory=list)       # [(i, j, k)] or [(i, j, k, ftype)]
    dihedrals: list = field(default_factory=list)    # [(i,j,k,l)] or [(i,j,k,l,type_name)]
    impropers: list = field(default_factory=list)    # [(i, j, k, l)]
    cmap: list = field(default_factory=list)         # [(i, j, k, l, m)]
