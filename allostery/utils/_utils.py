"""General useful functions"""

import os
from random import choice
from string import ascii_lowercase, digits
import subprocess


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

