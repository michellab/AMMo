import os
from argparse import ArgumentParser
from allostery.equilibrium import run_eq_md


def run_seeded_md(folder, snapshot, duration, report=5000, workdir=None, clean=False):
    os.chdir(folder)
    if not os.path.exists(f'snapshot_{snapshot}'):
        os.mkdir(f'snapshot_{snapshot}')
    
    # find topology
    if os.path.exists('system.prm7'):
        topology = 'system.prm7'
    else:
        topology = '../system-setup/system.prm7'

    #run seeded MD (i.e. equilibrium md using the seed input)
    run_eq_md(duration, topology, f'snapshots/snapshot_{snapshot}.rst7', f'snapshot_{snapshot}/production', report, os.path.abspath(workdir), clean)
    return None

def __main__():
    parser = ArgumentParser(description='Run a seeded MD simulation')
    parser.add_argument('--folder', type=str, required=True, help='Seeded MD folder, containing topology and a snapshots directory')
    parser.add_argument('--snapshot', type=str, required=True, help='Seed index')
    parser.add_argument('--duration', type=float, required=True, help='Seeded MD duration in ns')
    parser.add_argument('--report', type=int, default=5000, help='Report interval in steps. Default: 2500')
    parser.add_argument('--workdir', type=str, help='Working directory')
    parser.add_argument('--clean', action='store_true', help='Remove unneeded process files')
    args = parser.parse_args()

    run_seeded_md(args.folder, args.snapshot, args.duration, args.report, args.workdir, args.clean)

    return None


if __name__ == '__main__':
    __main__()