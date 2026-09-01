"""dvbfixer puppet — Strip PDB to backbone-only polyglycine ("puppet" model).

Removes all non-ATOM lines, strips sidechains, renames all residues to GLY.
Residues can be kept intact (all atoms, original name) via --keep.
Useful for creating minimal backbone scaffolds for modeling or alignment.
"""

import argparse
import sys
from pathlib import Path

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}


def _parse_keep(specs):
    """Parse --keep specs into a set of (chain, resseq) tuples.

    Formats:
        A:100       single residue
        A:100-110   range
        A:100,105   list
        A:100-110,120,130-135   mixed
    """
    kept = set()
    for spec in specs:
        chain, rest = spec.split(':', 1)
        for part in rest.split(','):
            if '-' in part:
                start, end = part.split('-', 1)
                for r in range(int(start), int(end) + 1):
                    kept.add((chain, r))
            else:
                kept.add((chain, int(part)))
    return kept


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog='dvbfixer puppet',
        description='Strip PDB to backbone-only polyglycine model.',
    )
    io = p.add_argument_group('Input / output')
    io.add_argument('input', help='Input PDB or PDBx/mmCIF file')
    io.add_argument('-o', '--output', help='Output PDB (default: <input>_puppet.pdb)')

    content = p.add_argument_group('Content selection')
    content.add_argument('--keep', action='append', default=[],
                         help='Keep residue(s) intact (all atoms, original name). '
                              'Format: CHAIN:NUM, CHAIN:START-END, or '
                              'CHAIN:NUM1,NUM2,START-END (repeatable)')

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(p, batch=True)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)

    if not input_path.exists():
        print(f"Error: {input_path} not found", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else \
        input_path.with_stem(input_path.stem + '_puppet')

    keep_residues = _parse_keep(args.keep) if args.keep else set()

    kept = 0
    removed = 0
    kept_intact = 0
    with open(input_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if not line.startswith('ATOM'):
                continue

            chain = line[21]
            resseq = int(line[22:26])

            if (chain, resseq) in keep_residues:
                fout.write(line)
                kept_intact += 1
                continue

            atomname = line[12:16].strip()
            if atomname not in BACKBONE_ATOMS:
                removed += 1
                continue
            # Rename residue to GLY
            line = line[:17] + 'GLY' + line[20:]
            fout.write(line)
            kept += 1

    msg = f"Wrote {output_path.name}: {kept} backbone atoms"
    if kept_intact:
        msg += f", {kept_intact} kept intact"
    msg += f" ({removed} sidechain/H atoms removed)"
    print(msg)


if __name__ == '__main__':
    main()
