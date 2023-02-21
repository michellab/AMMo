from argparse import ArgumentParser
from ast import literal_eval
from ammo.steering import run_smd


def __main__():
    parser = ArgumentParser(description='Run a steered MD simulation')
    parser.add_argument('--topology', type=str, required=True, help='system topology')
    parser.add_argument('--coordinates', type=str, required=True, help='equilibrated coordinates')
    parser.add_argument('--masks', type=str, required=True, help='a list of AMBER masks for each CV')
    parser.add_argument('--types', type=str, required=True, help='CV types')
    parser.add_argument('--timings', type=str, required=True, help='steering schedule in ns')
    parser.add_argument('--values', type=str, required=True, help='CV values in default PLUMED units at each point in "timings", with each list corresponding to an individual CV (unless there is only one point in "timings", in which case a single list of values is to be provided). "initial" will be replaced by a computed initial value. In additional to that, addition, subtraction, multiplication and division are allowed between the initial value and a number (e.g. "initial/2" or "5+initial") to help with multiple step steering protocols')
    parser.add_argument('--forces', type=str, required=True, help='forces to be applied to each CV in kJ/mol, in the same format as "--values"')
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
    
    args.masks = literal_eval(args.masks)
    args.types = literal_eval(args.types)
    args.timings = literal_eval(args.timings)
    args.values = literal_eval(args.values)
    args.forces = literal_eval(args.forces)

    run_smd(args.topology, args.coordinates, args.masks, args.types, args.timings, args.values, args.forces, args.reference, args.engine, args.workdir)

    return None


if __name__ == '__main__':
    __main__()
