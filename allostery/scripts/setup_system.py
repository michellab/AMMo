from ast import literal_eval
from argparse import ArgumentParser
from allostery.setup import setup_system


def __clean_args(args):
    # change protocol input into a list
    args.protocol = literal_eval(args.protocol)

    # fix charges is given
    if args.charges is not None:
        args.charges = literal_eval(args.charges)

    # fix extra parameters is given
    if args.parameters is not None:
        args.parameters = args.parameters.split(',')

    return args


def __main__():
    parser = ArgumentParser(description='General system setup: parameterisation, solvation, minimisation, heating, and equilibration.')
    parser.add_argument('--input', required=True, type=str, help='System PDB file')
    parser.add_argument('--protocol', required=True, type=str, help='Comma separated list of minimisation steps, heating duration in ps, '
                                                     'and equilibration duration in ps')
    parser.add_argument('--engine', type=str, default='GROMACS', help='Simulation engine supported by BioSimSpace. Default: "GROMACS"')
    parser.add_argument('--charges', type=str, help='Ligand charges in the order they appear in the input PDB, '
                                                    'comma separated')
    parser.add_argument('--parameters', type=str, help='any additional parameter arguments to give LeAP separated by comma')
    parser.add_argument('--topology', type=str, help='Dry topology of system. If provided will be used instead of '
                                                     're-parameterising')
    parser.add_argument('--solvated', action='store_true', help='Indicate that the system is already solvated')
    args = parser.parse_args()
    args = __clean_args(args)

    setup_system(args.input, args.protocol, args.engine, args.charges, args.parameters, args.topology, args.solvated)

    return None


if __name__ == '__main__':
    __main__()
