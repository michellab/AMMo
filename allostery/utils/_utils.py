"""General useful functions"""

import os
from random import choice
from string import ascii_lowercase, digits
import subprocess
import numpy as np


def get_dry_trajectory(topology, trajectory, output):
    """Dry a trajectory using cpptraj

    Parameters
    ----------
    topology : str
        path to topology file

    trajectory : str
        path to trajectory file

    output : str
        where to save the dry output

    Return
    ------
    None
    """
    to_strip = ':WAT,:SOL,:HOH,:NA,:CL,:Na+,:Cl-,:CLA,:POT,:SOD'
    cpptraj_file = f'cpptraj_{"".join(choice(ascii_lowercase + digits) for _ in range(10))}.in'

    with open(cpptraj_file, 'w') as file:
        file.writelines([f'parm {topology}\n',
                        f'trajin {trajectory}\n',
                        f'strip {to_strip}\n',
                        f'trajout {output}\n',
                         'go\nquit\n'])

    subprocess.run(['cpptraj', '-i', cpptraj_file])
    os.remove(cpptraj_file)

    return None

def __parse_seeds(seeds):
    if seeds == 'all':
        seeds = [int(folder.split('_')[1]) for folder in os.listdir() if folder.startswith('snapshot_') and '.pdb' not in folder]
        seeds.sort()
    elif isinstance(seeds, (list, tuple)):
        if not all(isinstance(idx, int) for idx in seeds):
            raise TypeError('Seed indices must all be of type int')
        else:
            return seeds
    elif isinstance(seeds, str):
        if '-' in seeds:
            idx_range = [int(idx) for idx in seeds.split('-')]
            return np.arange(idx_range[0], idx_range[1]+1).tolist()
        elif ',' in seeds:
            return [int(idx) for idx in seeds.split(',')]
        else:
            return [int(seeds)]
    elif isinstance(seeds, int):
        return [seeds]
    else:
        raise TypeError('Seeds must be of type str, int, list or tuple')