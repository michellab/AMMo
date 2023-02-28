import os
import subprocess
import BioSimSpace as BSS
import pytraj as pt
from shutil import move
from ammo.utils import get_dry_trajectory
from time import sleep


def __load_input(input):
    # first check if file
    if isinstance(input, str):
        if os.path.exists(input):
            with open(input, 'r') as file:
                plumed = file.readlines()
        else:
            raise ValueError(f'If "plumed" is a string, it has to point to a valid file.')
    # otherwise must be input string
    else:
        plumed = input

    return plumed


def __get_duration(plumed):
    """Find total number of steps in the plumed file"""
    reading_steps = False
    duration = 0
    for line in plumed:
        if 'MOVINGRESTRAINT' in line and not reading_steps:
            reading_steps = True
        elif 'MOVINGRESTRAINT' in line and reading_steps:
            reading_steps = False
            break
        if 'STEP' in line and reading_steps:
            parts = line.split()
            for part in parts:
                if 'STEP' in part:
                    step_duration = int(part.split('=')[1])
                    duration = max(duration, step_duration)
    return duration


def __create_process(system, duration, engine, workdir):
    # create a simple distance CV
    cv = BSS.Metadynamics.CollectiveVariable.Distance(1, 2)
    schedule = [0*BSS.Units.Time.nanosecond, duration*0.002*BSS.Units.Time.picosecond] # timestep is 0.002 ps
    restraints = [BSS.Metadynamics.Restraint(1*BSS.Units.Length.angstrom), BSS.Metadynamics.Restraint(1*BSS.Units.Length.angstrom)]

    # create the protocol
    protocol = BSS.Protocol.Steering(cv, schedule, restraints, runtime=schedule[-1], report_interval=2500, restart_interval=2500)
    
    # create the process
    workdir = os.path.abspath(workdir)
    if engine == 'AMBER':
        process = BSS.Process.Amber(system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda', work_dir=workdir) # specify using pmemd.cuda
    else:
        process = BSS.Process.createProcess(system, protocol, engine, work_dir=workdir)
    
    # go to the process workdir
    os.chdir(workdir)

    return protocol, process


def __edit_masks(plumed, system, pytraj_system, reference):
    rmsd_count = 0
    # loop over plumed input and replace masks with atoms
    # and "initial" values with actual system values
    for i, line in enumerate(plumed):
        # replace mask with atom indices
        if 'ATOMS=' in line:
            new_line = __parse_masks(line, pytraj_system)
            plumed[i] = new_line
        # write reference
        elif 'REFERENCE=' in line:
            if isinstance(reference, list):
                ref = reference[rmsd_count]
            else:
                ref = reference
            rmsd_count += 1
            new_line = __write_reference(line, system, ref, rmsd_count)
            plumed[i] = new_line
        # reach steering part of the file
        elif 'MOVINGRESTRAINT' in line:
            # write plumed file with only the CVs
            with open('plumed.dat', 'w') as file:
                file.writelines(plumed[:i]+['PRINT ARG=* FILE=initial.dat\n'])
            # run plumed driver to find the initial values
            subprocess.run(['plumed', 'driver', '--mf_pdb', 'input.pdb'])
            break

    return plumed


def __parse_masks(line, system):
    """Parse masks to return atom indices for the system"""
    parts = line.replace('\n', '').split() # remove newline character

    for i, part in enumerate(parts):
        if 'ATOMS' in part:
            mask = part.split('=')[1]
            atoms = system.topology.select(mask)+1 # add 1 since pytraj.Topology.select() indexes from 0 but PLUMED does from 1
            parts[i] = f'ATOMS={",".join([str(idx) for idx in atoms])}'
    
    new_line = ' '.join(parts) + '\n'

    return new_line


def __write_reference(line, system, reference, count):
    """Write a correctly ordered RMSD reference file"""
    parts = line.replace('\n', '').split() # remove newline character

    for i, part in enumerate(parts):
        if 'REFERENCE' in part:
            mask = part.split('=')[1]
            parts[i] = f'REFERENCE=reference_{count}.pdb'
    
    reference_atoms = pt.load(reference).topology.select(mask).tolist() # do not need to add 1 because indices will be handled by BSS
    reference = BSS.IO.readMolecules(reference).getMolecule(0)

    cv = BSS.Metadynamics.CollectiveVariable.RMSD(system, reference, reference_atoms, 0)
    with open(f'reference_{count}.pdb', 'w') as file:
        file.writelines([line+'\n' for line in cv.getReferencePDB()])

    new_line = ' '.join(parts) + '\n'

    return new_line


def __edit_values(plumed):
    # read in the initial values
    with open('initial.dat', 'r') as file:
        initial_contents = file.readlines()
    all_labels = initial_contents[0].split()[3:]
    all_values = [float(val) for val in initial_contents[-1].split()[1:]]
    values = {}
    for i in range(len(all_labels)):
        values[all_labels[i]] = all_values[i]

    # go over plumed file
    reading_restraint = False
    for i, line in enumerate(plumed):
        # find MOVINGRESTRAINT part
        if 'MOVINGRESTRAINT' in line and not reading_restraint:
            reading_restraint = True
        # reach end of restraint
        elif 'MOVINRESTRAINT' in line and reading_restraint:
            reading_restraint = False
            break
        # find which CVs are used for steering
        if 'ARG=' in line and reading_restraint:
            parts = line.replace('\n', '').split()
            for part in parts:
                if 'ARG=' in part:
                    labels = part.replace('ARG=', '').split(',')
        # find step definitions
        # and check if there are "initial"
        # values that need to be replaces
        elif 'STEP' in line and 'initial' in line and reading_restraint:
            parts = line.replace('\n', '').split()
            for j, part in enumerate(parts):
                if 'AT' in part and 'initial' in part:
                    step_values = part.split('=')[1].split(',')
                    for k, val in enumerate(step_values):
                        if 'initial' in val:
                            step_values[k] = str(__parse_operation(val, values[labels[k]]))
                    parts[j] = f'{part.split("=")[0]}={",".join(step_values)}'
            plumed[i] = '  ' + ' '.join(parts) + '\n'
    
    return plumed


def __parse_operation(expression, initial):
    """Parse a simple expression involving the initial CV value"""
    expression = expression.replace('initial', str(initial))
    if '/' in expression:
        components = [float(component) for component in expression.split('/')]
        value = components[0]/components[1]
    elif '*' in expression:
        components = [float(component) for component in expression.split('*')]
        value = components[0]*components[1]
    elif '+' in expression:
        components = [float(component) for component in expression.split('+')]
        value = components[0]+components[1]
    elif '-' in expression:
        components = [float(component) for component in expression.split('-')]
        value = components[0]-components[1]
    else:
        value = initial
    return value


def run_smd(topology, coordinates, input, reference=None, engine='AMBER', workdir='.'):
    """
    Run a steered MD simulation with AMBER or GROMACS and PLUMED.

    Parameters
    ----------
    topology : str
        topology file
    coordinates : str
        equilibrated system coordinate file
    input : str, [str]
        pseudo PLUMED input, either as the input itself, or as a path to a file containing it. ATOMS are replaced by AMBER masks, and CV values can be "initial" and with simple arithmetic operations, e.g. "initial/2"
    reference : str
        path to reference PDB file if using RMSD as a CV. All of the atoms in the reference have to also appear in the system, but not vice versa.
    engine : str
        MD engine to run sMD with. Allowed 'AMBER' and 'GROMACS'
    workdir : str
        working directory for MD. If None, it will be run in a temporary directory and copied to the current directory

    Returns
    -------
    process : BioSimSpace.Process
        BioSimSpace process used for steering
    """
    # load system files
    system = BSS.IO.readMolecules([topology, coordinates])
    pytraj_system = pt.load(coordinates, top=topology) # using pytraj system for searching and some initial value computing
    if isinstance(reference, list):
        reference = [os.path.abspath(ref) for ref in reference]
    elif isinstance(reference, str):
        reference = os.path.abspath(reference)

    # load input
    plumed = __load_input(input)
    duration = __get_duration(plumed)

    # create protocol and process
    protocol, process = __create_process(system, duration, engine, workdir)
    # now in process workdir
    # write input as PDB
    BSS.IO.saveMolecules('input', system, 'pdb')

    # fix plumed
    plumed = __edit_masks(plumed, system, pytraj_system, reference)
    plumed = __edit_values(plumed)
    with open('plumed.dat', 'w') as file:
        file.writelines(plumed)

    # run process
    process.start()
    print(f'Process running in {process.workDir}')
    process.wait()

    to_copy = {'amber.nc': 'steering.nc', 'plumed.dat': 'plumed.dat', 'COLVAR': 'steering.dat', 'amber.out': 'steering.out'}
    if isinstance(reference, list):
        for i in range(len(reference)):
            to_copy[f'reference_{i+1}.pdb'] = f'reference_{i+1}.pdb'
    else:
        to_copy['reference_1.pdb'] = 'reference_1.pdb'
    for file in to_copy:
        move(f'{process.workDir()}/{file}', to_copy[file])

    # dry trajectory
    get_dry_trajectory(topology, 'steering.nc', 'steering_dry.nc')

    return None
