from allostery.equilibrium import run_eq_md
from argparse import ArgumentParser


def __main__():
    parser = ArgumentParser(description='Run equilibrium MD with BioSimSpace')
    parser.add_argument('--duration', type=float, required=True, help='MD simulation duration in ns')
    parser.add_argument('--topology', type=str, help='parameter file')
    parser.add_argument('--coordinates', type=str, help='cooridnate file')
    parser.add_argument('--output', type=str, help='output location, without extension (e.g. "production-1"). If None, the next available name will be used (e.g. "production-3" if "production-1" and "production-2" already exist)')
    args = parser.parse_args()

    run_eq_md(args.duration, args.topology, args.coordinates, args.output)
    return None


if __name__ == '__main__':
    __main__()