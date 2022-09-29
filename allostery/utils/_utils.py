"""General useful functions"""

import os
from warnings import warn
from random import choice
from string import ascii_lowercase, digits
import subprocess
import numpy as np
from pytraj import load, Trajectory


def renumber_pdb(input, reference, output=None, matching='warn', offset=0):
    """Renumber a PDB file based on a reference structure.

    Parameters
    ----------
    input : str
        PDB input to renumber

    reference : pytraj.Trajectory, str, [str]
        reference system files

    output : str
        output path. If None, no output will be saved

    matching : str
        how to handle when not all atoms in "input" are found in "reference". Allowed values are
        "ignore", "warn", and "error"

    offset : int
        residue difference between the input and reference

    Returns
    -------
    new_pdb : [str]
        renumbered PDB file
    """
    # load reference
    if isinstance(reference, Trajectory):
        return reference
    elif isinstance(reference, str):
        if not reference.endswith('.pdb'):
            raise AssertionError('If only a single file provided for reference, it must be of PDB format')
        else:
            reference = load(reference)
    elif isinstance(reference, list):
        for file in reference:
            if file.endswith('.prm7') or file.endswith('.parm7') or file.endswith('.top'):
                topology = file
            else:
                coordinates = file
        reference = load(coordinates, top=topology)
    else:
        raise TypeError('Unsupported reference files')

    # read pdb
    with open(input, 'r') as file:
        pdb = file.readlines()

    # look for PDB atoms in reference
    new_pdb_unsorted = {}
    unmatched = []
    for line in pdb:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            res_idx = int(line[22:26].strip())
            at_name = line[12:16]
            at_idx = int(line[6:11].strip())
            search = reference.top.select(f':{res_idx+offset}@{at_name}')
            if len(search) == 0:
                unmatched.append(f':{res_idx}@{at_name}')
                new_pdb_unsorted[at_idx] = line
                last_idx = at_idx
            else:
                new_pdb_unsorted[search[0]+1] = line[:6] + '%5s'%(search[0]+1) + line[11:]
                last_idx = search[0]+1
        else:
            new_pdb_unsorted[last_idx] = new_pdb_unsorted[last_idx]+line
    
    # deal with unmatched atoms
    if len(unmatched) > 0:
        if matching == 'warn':
            print(f'Warning: residues {",".join(unmatched)} in system were not found in reference')
        elif matching == 'error':
            raise ValueError(f'residues {",".join(unmatched)} in system were not found in reference')
    
    # sort atoms by new index
    new_pdb = [new_pdb_unsorted[idx] for idx in sorted(new_pdb_unsorted)]

    # write output
    if output is not None:
        with open(output, 'w') as file:
            file.writelines(new_pdb)

    return new_pdb


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