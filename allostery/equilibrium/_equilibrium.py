import os
from shutil import copyfile
from re import match
import BioSimSpace as BSS
from allostery.utils import get_dry_trajectory


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

def run_eq_md(duration, topology, coordinates, output=None):
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
    Returns
    -------
    None
    """   
    system = BSS.IO.readMolecules([topology, coordinates])

    if output is None:
        output = __get_output_location()
    
    #set up process
    protocol = BSS.Protocol.Production(runtime=duration*BSS.Units.Time.nanosecond, restart_interval=2500, report_interval=2500)
    process = BSS.Process.Amber(system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda')

    # run process
    process.start()
    process.wait()
    
    #save results
    files = ['nc', 'rst7', 'out']
    for  ext in files:
        copyfile(f'{process.workDir()}/amber.{ext}', f'{output}.{ext}')
    #dry trajectory
    get_dry_trajectory(f'{process.workDir()}/amber.prm7', f'{process.workDir()}/amber.nc', f'{output}_dry.nc')