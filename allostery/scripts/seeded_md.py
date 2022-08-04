import os
from argparse import ArgumentParser
from allostery.equilibrium import run_eq_md


def run_seeded_md(folder, snapshot, duration):
    os.chdir(folder)
    if not os.path.exists(f'snapshot_{snapshot}'):
        os.mkdir(f'snapshot_{snapshot}')
    
    # find topology
    if os.path.exists('system.prm7'):
        topology = 'system.prm7'
    else:
        topology = '../system-setup/system.prm7'

    #run seeded MD (i.e. equilibrium md using the seed input)
    run_eq_md(duration, topology, f'snapshots/snapshot_{snapshot}.rst7', f'snapshot_{snapshot}/production')
    return None

def __main__():
    parser = ArgumentParser(description='Run a seeded MD simulation')
    parser.add_argument('--folder', type=str, required=True, help='Seeded MD folder, containing topology and a snapshots directory')
    parser.add_argument('--snapshot', type=str, required=True, help='seed index')
    parser.add_argument('--duration', type=float, required=True, help='seeded MD duration in ns')
    args = parser.parse_args()

    run_seeded_md(args.folder, args.snapshot, args.duration)

    return None


if __name__ == '__main__':
    __main__()