import os
from shutil import move
from re import match
import BioSimSpace as BSS
from ammo.utils import get_dry_trajectory


from time import sleep

def __get_output_location():
    """Find how many production runs are already in equilibrium MD
    location and return the next one. This is consistent with how I
    name files.
    Returns
    -------
    output : str
        the path of next production output name, e.g. production-3
        does not include extension for easier manipulation
    """
    files = [file for file in os.listdir() if match(f'production-[0-9].out',file)]
    start = 1
    for file in files:
        current = int(file.split('.')[0].split('-')[1].split('_')[0])
        if current>=start:
            start = current+1
    output = f'production-{start}'
    
    return output

def run_eq_md(duration, topology, coordinates, output=None, report=2500, workdir=None, clean=False):
    """Run equilibrium MD script using BioSimSpace.
    Parameters
    ----------
    duration : float
        MD simulation duration in ns
    topology : str
        system topology
    coordinates : str
        system coordinates
    output : str
        output location, without extension (e.g. 'production-1'). If None, the next available name will
        be used (e.g. 'production-3' if 'production-1' and 'production-2' already exist)
    report : int
        report interval
    workdir : str
        workding directory
    clean : bool
        whether to remove unneeded process files

    Returns
    -------
    None
    """
    print(topology,coordinates)   
    system = BSS.IO.readMolecules([topology, coordinates])

    if output is None:
        output = __get_output_location()
    
    # set up process
    protocol = BSS.Protocol.Production(runtime=duration*BSS.Units.Time.nanosecond, restart_interval=report, report_interval=report)
    process = BSS.Process.Amber(system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda', work_dir=workdir)

    # run process
    process.start()
    process.wait()

    # dry trajectory
    get_dry_trajectory(f'{process.workDir()}/amber.prm7', f'{process.workDir()}/amber.nc', f'{output}_dry.nc')
    
    # save results
    files = {'nc': 'nc', 'crd': 'rst7', 'out': 'out'}
    for  src, dest in files.items():
        move(f'{process.workDir()}/amber.{src}', f'{output}.{dest}')

    # clean files
    to_remove = ['amber.prm7', 'amber.rst7', 'README.txt', 'amber.err', 'amber.nrg', 'amber.cfg']
    if clean and workdir is not None:
        for file in to_remove:
            os.remove(f'{process.workDir()}/{file}')
