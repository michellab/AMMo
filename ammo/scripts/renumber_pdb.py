from argparse import ArgumentParser
from ammo.utils import renumber_pdb


def __main__():
    parser = ArgumentParser(description='Renumber the atom indices in a PDB file based on a reference')
    parser.add_argument('--input', type=str, required=True, help='Input PDB file')
    parser.add_argument('--reference', type=str, required=True, help='Reference files')
    parser.add_argument('--output', type=str, required=True, help='Output path')
    parser.add_argument('--matching', type=str, default='warn', help='How to handle when not all atoms in "input" are found in "reference". Allowed values are "ignore", "warn", and "error". Default: "warn"')
    parser.add_argument('--offset', type=int, default=0, help='Residue difference between the input and reference. Default: 0')
    args = parser.parse_args()

    if ',' in args.reference:
        args.reference = args.reference.split(',')
    
    renumber_pdb(args.input, args.reference, args.output, args.matching, args.offset)


if __name__ == '__main__':
    __main__()