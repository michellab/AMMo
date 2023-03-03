import os
from random import choice
from string import ascii_lowercase
import subprocess
from numpy import loadtxt


def featurize(trajectory, topology, feature, mask, reference=None, shared='!@/H'):
    """Reduce trajectory data to a single feature.
    
    Parameters
    ----------
    trajectory : str, pytraj.Trajectory
        trajectory to be featurized. If 'str', 'topology' must be specified
        and it will be loaded. 
       
    feature : str
        feature type. Allowed: 'rmsd', 'distance', 'torsion'
        
    mask : str
        AMBER selection mask for the feature
        
    reference : str, pytraj.Trajectory
        reference for RMSD calculations. If 'str' of a PDB file, it will be
        loaded
        
    shared : str
        AMBER selection mask for shared atoms between trajectory and reference.
        If 'None', all atoms will be used. Default is not H, meaning all heavy atoms
        
    topology : str
        topology of the system
        
    Returns
    -------
    featurized : numpy.array
        a time series of the feature
    """
    cpptraj_input = [f'parm {os.path.abspath(topology)}\n',
                     f'trajin {os.path.abspath(trajectory)}\n']

    file = f'/tmp/cpptraj_{"".join(choice(ascii_lowercase) for i in range(10))}.in'
    output = f'/tmp/cpptraj_{"".join(choice(ascii_lowercase) for i in range(10))}.out'
        
    # calculate the feature
    if feature == 'rmsd':
        cpptraj_input += [f'parm {os.path.abspath(reference)} name ref_parm\n',
                          f'reference {os.path.abspath(reference)} parm ref_parm name rmsd_ref\n',
                          f'rms {shared} ref rmsd_ref\n',
                          f'rmsd {mask} out {output} nofit ref rmsd_ref\n']
    else:
        cpptraj_input += [f'{feature} {mask} out {output}\n']
    cpptraj_input += ['go\nquit\n']

    # run cpptraj
    with open(file, 'w') as fl:
        fl.writelines(cpptraj_input)
    process_output = subprocess.run(['cpptraj', '-i', file], capture_output=True, text=True)
    os.remove(file)

    # check for errors
    error = 'Error: Error(s) occurred during execution.'
    if error in process_output.stderr:
        print(process_output.stderr.split('\n')[0])
        return None
    
    # load results
    data = loadtxt(output, comments='#')
    os.remove(output)
        
    return data[:,1]