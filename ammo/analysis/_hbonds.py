import os
from string import ascii_lowercase
from random import choice
from pandas import read_csv
from subprocess import run, DEVNULL
from pytraj import save, Trajectory


def hbonds(trajectory, mask="*", distance=3.5, angle=135, topology=None, filter=None, frequency=True, quiet=True):
    """Analyse the HBonds of a trajectory.

    Parameters
    ----------
    trajectory : str, pytraj.Trajectory
        the trajectory to analyse. If str, "topology" must also be set

    mask : str
        mask for a region to search, e.g. :100<:5 means hbonds between residues 5 A within residue 100

    distance : float
        hbond distance cutoff

    angle : float
        hbond angle cutoff

    topology : str
        topology file path

    filter : str
        filter detected HBond labels with a regex value

    frequency : bool
        return the HBonds as a df with HBond labels as index, and % frequency as the values

    quiet : bool
        whether to run cpptraj quietly. Consider turning this off if errors arise
    
    Returns
    -------
    df : pandas.DataFrame
        a data frame with the hbond information
    """
    # check that AMBERHOME is set
    __check_cpptraj()

    # get trajectory file
    trajectory, topology, remove = __get_trajectory(trajectory, topology)

    # write input file
    input_file = '/tmp/' + ''.join([choice(ascii_lowercase) for _ in range(10)]) + '.in'
    output_file = '/tmp/' + ''.join([choice(ascii_lowercase) for _ in range(10)]) + '.txt'
    with open(input_file, 'w') as file:
        file.writelines([f'parm {topology}\n',
                         f'trajin {trajectory}\n',
                         f'reference {trajectory} 1\n'
                         f'hbond mask {mask} angle {angle} dist {distance} series uuseries {output_file}\n',
                          'go\nquit\n'])
        
    if quiet:
        run(['cpptraj', '-i', input_file], stdout=DEVNULL)
    else:
        run(['cpptraj', '-i', input_file])
    df = read_csv(output_file, delim_whitespace=True)

    # remove tmp files
    tmp_files = [input_file, output_file]
    if remove:
        tmp_files += [topology, trajectory]
    for file in tmp_files:
        os.remove(file)

    # post process as required
    if filter is not None:
        df = df.filter(regex=filter, axis=1)
    if frequency:
        df = df.mean().iloc[1:]*100

    return df


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
        trajectory_file = '/tmp/' + ''.join([choice(ascii_lowercase) for _ in range(10)]) + '.nc'
        topology_file = '/tmp/' + ''.join([choice(ascii_lowercase) for _ in range(10)]) + '.prm7'
        save(trajectory_file, trajectory)
        save(topology_file, trajectory.top)
        return trajectory_file, topology_file, True
    else:
        raise TypeError(f'Unsupported trajectory type: {type(trajectory)}. Trajectory has to be str of a file path or a pytraj.Trajectory')