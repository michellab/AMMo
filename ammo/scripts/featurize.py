import numpy as np
from argparse import ArgumentParser
from ammo.analysis import featurize


def __main__():
    parser = ArgumentParser(description='Reduce trajectory data to a single feature. Run in the folder containing snapshot folders')
    parser.add_argument('--topology', required=True, help='system topology file')
    parser.add_argument('--trajectory', required=True, help='seeded MD trajectory file"')
    parser.add_argument('--feature', required=True, help='Type of feature to calculate. Allowed: "rmsd", "torsion", "distance"')
    parser.add_argument('--mask', required=True, help='AMBER selection mask for the feature')
    parser.add_argument('--output', required=True, help='output file for feature')
    parser.add_argument('--reference', help='RMSD reference PDB file')
    parser.add_argument('--shared', default='!@/H', help='mask for atoms used for aligning the RMSD reference. Default: "!@/H"')
    args = parser.parse_args()

    featurized = featurize(args.trajectory, args.topology, args.feature, args.mask, args.reference, args.shared)
    np.savetxt(args.output, featurized)
    
    return None


if __name__ == '__main__':
    __main__()