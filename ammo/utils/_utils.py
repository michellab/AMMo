"""General useful functions"""

import os
from random import choice
from string import ascii_lowercase, digits
import subprocess
import numpy as np
from pytraj import load, Trajectory, save
from BioSimSpace.Process import Amber


def __random_name(n_chars=10):
    return "".join(choice(ascii_lowercase + digits) for _ in range(n_chars))


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
    cpptraj_file = f'cpptraj_{__random_name()}.in'

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
    if isinstance(seeds, (list, tuple)):
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


def __parse_time(input, output_units, output_digits=3, output_type='string'):
    """Parse time given as a string and return either as string with units or numerical value

    Parameters
    ----------
    input : str
        time in format "value unit", e.g. "10 ps"

    output_units : str
        units to give output in. Allowed values are: "fs", "ps", "ns, "us", "ms", and "s"

    output_digits : int
        place to round output to

    output_type : str
        type of output to give. If "string" it will be "value unit" and if "number" it will be the numerical value 
    """
    conversion = {'fs': -15, 'ps': -12, 'ns': -9, 'us': -6, 'ms': -3, 's': 0}

    components = input.split(' ')

    try:
        value = int(components[0]) # get numerical value
        value = value * 10**(conversion[components[1]]) # convert to seconds
    except:
        raise ValueError(f'Time has to be in the format "value unit", e.g. "10 ps". Time provided: {input}')
    
    output_value = round(value * 10**(-1*conversion[output_units]), output_digits)

    if output_type == 'string':
        return f'{output_value} {output_units.replace("u", "μ")}'
    elif output_type == 'number':
        return output_value
    else:
        raise ValueError(f'"output_type" can be either "string" or "number". Was: {output_type}')
    

def __check_cpptraj():
    if 'AMBERHOME' not in os.environ:
        raise LookupError('AMBERHOME not defined. Cannot run cpptraj')
    else:
        return None


def __get_trajectory(trajectory, topology):
    if isinstance(trajectory, str) and topology is None:
        raise ValueError('If providing a trajectory file path, a topology file is also required')
    elif isinstance(trajectory, str):
        return os.path.abspath(trajectory), os.path.abspath(topology), False
    elif isinstance(trajectory, Trajectory):
        trajectory_file = '/tmp/' + __random_name() + '.nc'
        topology_file = '/tmp/' + __random_name() + '.prm7'
        save(trajectory_file, trajectory)
        save(topology_file, trajectory.top)
        return trajectory_file, topology_file, True
    else:
        raise TypeError(f'Unsupported trajectory type: {type(trajectory)}. Trajectory has to be str of a file path or a pytraj.Trajectory')


def __add_restraint(process, restraint, system):
    # check that AMBER process
    if not isinstance(process, Amber):
        raise TypeError(f'Restraints supported for AMBER simulations only')

    # replace masks in restraint file with atom indices
    # read in file if str to path
    if isinstance(restraint, str):
        with open(restraint, 'r') as file:
            restraint = file.readlines()
    
    # loop over each line
    for i, line in enumerate(restraint):
        if 'iat' in line:
            parts = line.replace('\n', '').split() # remove newline character and split
            for j, part in enumerate(parts):
                if part != 'iat=' and part != 'iat' and part != '=': # check if not reasonable pointers to atoms. This means the restraint file has to be tidy, with atoms on a separate line
                    atom = system.topology.select(part)+1 # add 1 to index from 1
                    if len(atom) == 0:
                        raise ValueError(f'Restraint atom {part} not found')
                    elif len(atom) > 1:
                        raise ValueError(f'More than 1 restraint atom {part} found')
                    else:
                        parts[j] = f'{atom[0]},'
            restraint[i] = '  ' + ' '.join(parts) + '\n' # rejoin and add a newline character

    # write the restraint file to the process directory
    with open(f'{process.workDir()}/RST', 'w') as file:
        file.writelines(restraint)

    # change the process config file
    config = process.getConfig()[:-1] # remove the "/" character
    config += ["  nmropt=1,",
               " /",
               " &wt type='REST', istep1=0,istep2=3000,value1=0.1,value2=1.0,  /",
               " &wt type='REST', istep1=3000,istep2=0,value1=1.0,value2=1.0  /",
               "&wt type='END'  /",
               "DISANG=RST"]
    process.setConfig(config)

    return None