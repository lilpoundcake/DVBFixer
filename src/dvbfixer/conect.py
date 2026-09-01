"""dvbfixer conect — add inferred CONECT records to a PDB file.

Reads an input PDB, infers connectivity via OpenBabel + domain overrides
(SS bonds, sugar glycosidic linkages, ASN/SER/THR glycosylation), merges
with any existing CONECT records, and writes the result to a new PDB.

The same inference engine runs automatically inside `prepare`, `top`,
`minimize`, `transplant`, and `convert` (gated by `--no-infer-conect`).
Run this subcommand explicitly when you want CONECT as a discrete
preprocessing step or to inspect what bonds dvbfixer perceives.
"""

import argparse
import sys
from pathlib import Path

from dvbfixer.pdbutils import (
    _has_any_conect,
    _strip_existing_conect,
    infer_conect_records,
    write_conect_block,
)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='dvbfixer conect',
        description='Infer and write CONECT records for a PDB file.',
    )
    io = parser.add_argument_group('Input / output')
    io.add_argument('input', help='Input PDB, PDBx/mmCIF, or crystallographic CIF file')
    io.add_argument('-o', '--output',
                    help='Output PDB file (default: <input>_conect.pdb)')
    io.add_argument('--force', action='store_true',
                    help='Allow in-place overwrite (when --output equals input)')

    content = parser.add_argument_group('Content selection')
    content.add_argument('--include-protein-backbone', action='store_true',
                         help='Also emit CONECT for standard amino-acid backbone '
                              'bonds. Off by default (FF templates own those).')

    diag = parser.add_argument_group('Diagnostics')
    diag.add_argument('-v', '--verbose', action='store_true',
                      help='Print bond counts and source breakdown')

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(parser, batch=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    in_path = Path(args.input)
    if not in_path.exists():
        print(f"Error: {in_path} not found", file=sys.stderr)
        sys.exit(1)

    if args.output:
        out_path = Path(args.output)
    else:
        out_path = in_path.with_name(in_path.stem + '_conect.pdb')

    if out_path.resolve() == in_path.resolve() and not args.force:
        print("Error: --output equals input. Pass --force to overwrite.",
              file=sys.stderr)
        sys.exit(1)

    has_existing = _has_any_conect(in_path)
    if args.verbose:
        print(f"Input: {in_path}")
        print(f"  has existing CONECT: {has_existing}")

    bonds = infer_conect_records(
        in_path,
        preserve_existing=True,
        include_protein_backbone=args.include_protein_backbone,
        verbose=args.verbose,
    )

    with open(in_path) as f:
        lines = f.readlines()
    stripped = _strip_existing_conect(lines)
    new_conect = write_conect_block(bonds)

    out_lines = []
    inserted = False
    for ln in stripped:
        if not inserted and ln.startswith(('END', 'ENDMDL')):
            out_lines.extend(new_conect)
            inserted = True
        out_lines.append(ln)
    if not inserted:
        out_lines.extend(new_conect)

    with open(out_path, 'w') as f:
        f.writelines(out_lines)

    print(f"Wrote {out_path} ({len(bonds)} CONECT bonds)")


if __name__ == '__main__':
    main()
