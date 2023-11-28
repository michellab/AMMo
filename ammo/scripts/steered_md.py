from argparse import ArgumentParser
from ammo.steering import run_smd


def __main__():
    parser = ArgumentParser(description='Run a steered MD simulation')
    parser.add_argument('--topology', type=str, required=True, help='system topology')
    parser.add_argument('--coordinates', type=str, required=True, help='equilibrated coordinates')
    parser.add_argument('--input', type=str, required=True, help='path to pseudo PLUMED file')
    parser.add_argument('--engine', type=str, default='AMBER', help='MD engine to run sMD with. Default: AMBER')
    parser.add_argument('--workdir', type=str, default='.', help='Working directory. If None, sMD will be run in a temporary folder and copied over. Default: "."')
    parser.add_argument('--suffix', type=int, help='A suffix to add to the output of the simulation, e.g. "steering_1.nc" or "steering_1.dat". For cases when the steering is done in more than one step. If None, nothing will be added')
    parser.add_argument('--restraint', type=str, help='A pseudo flat bottom restraint file that will be used during the steering (currently only available for AMBER). Instead of atom indices, AMBER masks are used')
    args = parser.parse_args()

    if args.workdir == 'None':
        args.workdir = None

    run_smd(args.topology, args.coordinates, args.input, args.engine, args.workdir, args.suffix, args.restraint)

    return None


if __name__ == '__main__':
    __main__()
