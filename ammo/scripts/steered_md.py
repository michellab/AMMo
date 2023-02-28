from argparse import ArgumentParser
from ammo.steering import run_smd


def __main__():
    parser = ArgumentParser(description='Run a steered MD simulation')
    parser.add_argument('--topology', type=str, required=True, help='system topology')
    parser.add_argument('--coordinates', type=str, required=True, help='equilibrated coordinates')
    parser.add_argument('--input', type=str, required=True, help='path to pseudo PLUMED file')
    parser.add_argument('--reference', type=str, help='path to reference PDB file(s) if using RMSD as a CV. All of the atoms in the reference have to also appear in the system')
    parser.add_argument('--engine', type=str, default='AMBER', help='MD engine to run sMD with. Default: AMBER')
    parser.add_argument('--workdir', type=str, default='.', help='Working directory. If None, sMD will be run in a temporary folder and copied over. Default: "."')
    args = parser.parse_args()

    if args.reference is not None:
        args.reference = args.reference.split(',')
        if len(args.reference) == 1:
            args.reference = args.reference[0]
    if args.workdir == 'None':
        args.workdir = None

    run_smd(args.topology, args.coordinates, args.input, args.reference, args.engine, args.workdir)

    return None


if __name__ == '__main__':
    __main__()
