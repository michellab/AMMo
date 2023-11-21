from argparse import ArgumentParser
from ammo.utils._utils import __parse_seeds as _parse_seeds
from ammo.analysis import compare_trajectories


def __main__():
    parser = ArgumentParser(description='Compare two trajectories in terms of their Hbond frequency and torsional angle divergences')
    parser.add_argument('--traj', required=True, help='trajectory file to analyse')
    parser.add_argument('--top', required=True, help='trajectory topology file')
    parser.add_argument('--ref_traj', required=True, help='reference trajectory file')
    parser.add_argument('--ref_top', required=True, help='reference trajectory topology file')
    parser.add_argument('--n_bins', type=int, default=100, help='number of bins for divergence calculation')
    parser.add_argument('--div_type', default='KL', help='divergence type. Allowed values: "KL" or "JS"')
    parser.add_argument('--filter', default=[20,20], help='top number (if int) or fraction (if float) of divergence/hbond frequency difference values to filter, in order of (hbonds, dihedrals); or the divergence/frequency cutoff(s). If None, all will be plotted/visualized')
    parser.add_argument('--filter_type', default='top', help='way of filtering. Allowed values are "top" (for top n or fraction values) and "cutoff" (for minimum divergence/frequency value(s))')
    parser.add_argument('--plot', default='plots', help='base file name to save the plots of filtered distributions in the working directory')
    parser.add_argument('--colors', default=['seagreen', 'indigo'], help='colors for plotting he trajectory and reference values respectively')
    parser.add_argument('--pymol', default='analysis', help='base filename to save a pymol session in the working directory. If None, no session will be saved')
    parser.add_argument('--residues', default=None, help='comma separated list or dash separated range of residues to analyse. If None, all residues will be included')
    parser.add_argument('--workdir', default=None, help='working directory for input/output files. If None, a directory will be created')
    parser.add_argument('--overwrite', action='store_true', help='overwrite files found in the working trajectory')
    args = parser.parse_args()

    if args.filter == 'None':
        args.filter = None
    elif ',' in args.filter:
        filter = []
        for val in args.filter.replace('[','').replace(']', '').split(','):
            if '.' in val:
                filter.append(float(val))
            else:
                filter.append(int(val))
        args.filter = filter

    if isinstance(args.colors, str):
        args.colors = args.colors.split(',')

    if args.pymol == 'None':
        args.pymol = None

    if args.residues == 'None':
        args.residues = None
    elif args.residues is not None:
        args.residues = _parse_seeds(args.residues)

    if args.workdir == 'None':
        args.workdir = None

    output = compare_trajectories(args.traj, args.ref_traj, args.top, args.ref_top, args.n_bins, args.div_type, args.filter, args.filter_type, args.plot, args.colors, args.pymol, args.residues, args.workdir, args.overwrite)

    return None


if __name__ == '__main__':
    __main__()