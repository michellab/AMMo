"""General useful functions"""

import os
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
    to_strip = ':WAT,:SOL,:NA,:CL,:Na+,:Cl-'
    with open('cpptraj_drying.in', 'w') as file:
        file.writelines([f'parm {topology}\n',
                        f'trajin {trajectory}\n',
                        f'strip {to_strip}\n',
                        f'trajout {output}\n',
                         'go\nquit\n'])

    subprocess.run(['cpptraj', '-i', 'cpptraj_drying.in'])
    os.remove('cpptraj_drying.in')

    return None

