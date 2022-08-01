import os

import BioSimSpace as BSS
import pytraj as pt
from shutil import copyfile
from allostery.steering import PlumedFile, RMSDReference
from allostery.utils import get_dry_trajectory
from time import sleep


def __validate_input(timings, values, forces):
    if len(timings) == 1: # if a single steering step
        # check if values given as a single list
        if not isinstance(values[0], list):
            values = [[val] for val in values]
        # check the same for forces
        if not isinstance(forces[0], list):
            forces = [[force] for force in forces]
    
    # add the standard start and end steps to all
    timings = [0, 0.004] + timings + [timings[-1]+2]
    values = [['initial', 'initial'] + val + [val[-1]] for val in values]
    forces = [[0, force[0]] + force + [0] for force in forces]

    return timings, values, forces


def __parse_masks(masks, system, types, reference=None):
    """
    Parse masks to return atom indices for the system

    Parameters
    ----------
    masks : [str], str
        a list of AMBER masks for each CV in the order they appear in the plumed file
    system : pytraj.Trajectory
        the system to be steered
    types : [str]
        CV types
    reference : str, [str]
        RMSD reference. Required if rmsd is one of the CV types

    Returns
    -------
    atoms : [[int]]
        a list of atom lists for each CV
    """
    if isinstance(types, str):
        types = [types]
    if isinstance(masks, str):
        masks = [masks]
    if len(types) != len(masks):
        raise ValueError(f'Number of types ({len(types)}) has to match masks ({len(masks)}).')
    if isinstance(reference, str):
        reference = [reference]*types.count('rmsd')

    atoms = []
    rmsd_count = 0
    for mask, cv_type in zip(masks,types):
        if cv_type != 'rmsd':
            results = system.topology.select(mask)
        else:
            results = pt.load(reference[rmsd_count]).topology.select(mask)
            rmsd_count+=1
        atoms.append(results.tolist())

    return atoms


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
    else:
        components = [float(component) for component in expression.split('-')]
        value = components[0]-components[1]
    return value


def __create_cvs(system, types, atoms, reference=None):
    """Create CVs for sMD simulations

    Parameters
    ----------
    system : BioSimSpace._SireWrappers.System
        system to be steered
  
    types : str, [str]
        CV types. Allowed values are 'rmsd', 'torsion', and 'distance'

    atoms : str, [str]
        atom indices for each CV

    reference : BioSimSpace._SireWrappers.Molecule, [BioSimSpace._SireWrappers.Molecule]
        RMSD reference. Required if rmsd is one of the CV types

    Returns
    -------
    cvs : [BioSimSpace.Metadyanmics.CollectiveVariable]
        the CVs created
    """
    cvs = []

    # validate input
    if isinstance(reference, BSS._SireWrappers.Molecule):
        reference = [reference]*types.count('rmsd')

    # loop and create cvs
    rmsd_count = 0 # keep track of RMSD
    for cv_type, indices in zip(types, atoms):
        if cv_type == 'distance':
            cv = BSS.Metadynamics.CollectiveVariable.Distance(indices[0], indices[1])
        elif cv_type == 'torsion':
            cv = BSS.Metadynamics.CollectiveVariable.Torsion(indices)
        elif cv_type == 'rmsd':
            cv = BSS.Metadynamics.CollectiveVariable.RMSD(system, reference[rmsd_count], indices, 0)
            rmsd_count+=1
        else:
            raise ValueError(f'CV type {type} not supported. Available types are: distance, torsion, rmsd')
        cvs.append(cv)
    
    return cvs


def __compute_value(system, cv_type, atoms):
    """Compute CV values using pytraj"""
    if cv_type == 'distance':
        result = pt.calc_distance(system, atoms)[0][0]/10
    elif cv_type == 'torsion':
        result = pt.calc_dihedral(system, atoms)[0][0]*0.017453
    else:
        raise TypeError(f'Unsupported CV: {cv_type}')
    return result


def __compute_restraints(cvs, values, forces, system):
    """
    Convert the CV values into steering restraints

    Parameters
    ----------
    cvs : [BioSimSpace.Metadyanmics.CollectiveVariable]
        the CVs created
    values : [[float,int,str]]
        CV values in default PLUMED units. 'initial' will be replaced by a computed initial value
    forces : [[float,int]]
        force constant for each CV at each step
    system : pytraj.Trajectory
        pytraj system

    Returns
    -------
    restraints : [[BSS.Metadynamics.Restraint]]
        steering restrants
    """
    unord_restraints = []

    # loop through each CV
    for cv, cv_values, cv_forces in zip(cvs, values, forces):
        # compute initial values and get units
        if isinstance(cv, BSS.Metadynamics.CollectiveVariable.RMSD):
            initial = cv.getInitialValue().nanometers().value()
            units = BSS.Units.Length.nanometer
        elif isinstance(cv, BSS.Metadynamics.CollectiveVariable.Distance):
            initial = __compute_value(system, 'distance', [cv.getAtom0(), cv.getAtom1()])
            units = BSS.Units.Length.nanometer
        else:
            initial = __compute_value(system, 'torsion', cv.getAtoms())
            units = BSS.Units.Angle.radian

        # create the restraints
        cv_restraints = []
        for val, force in zip(cv_values, cv_forces):
            if isinstance(val, str):
                if val == 'initial':
                    val = initial
                elif 'initial' in val: # allow a simple mathematical operation
                    val = __parse_operation(val, initial)
            cv_restraints.append(BSS.Metadynamics.Restraint(val*units, force))

        unord_restraints.append(cv_restraints)
    
    restraints = []
    for i in range(len(unord_restraints[0])):
        restraints.append([restraint[i] for restraint in unord_restraints])

    return restraints


def __create_protocol(cvs, timings, restraints):
    """Create protocol

    Parameters
    ----------
    cvs : [BioSimSpace.Metadyanmics.CollectiveVariable]
        the CVs
    timings : [int, float]
        steering schedule in ns
    restraints : [[BSS.Metadynamics.Restraint]]
        steering restrants

    Returns
    -------
    protocol : BSS.Protocol.Steering
        steering protocol
    """
    ns = BSS.Units.Time.nanosecond
    schedule = [time*ns for time in timings]
    protocol = BSS.Protocol.Steering(cvs, schedule, restraints, runtime=schedule[-1], report_interval=2500, restart_interval=2500)

    return protocol


def run_smd(topology, coordinates, masks, types, timings, values, forces, reference=None, engine='AMBER'):
    """
    Run a steered MD simulation with AMBER or GROMACS and PLUMED.

    Parameters
    ----------
    topology : str
        topology file
    coordinates : str
        equilibrated system coordinate file
    masks : [str], str
        a list of AMBER masks for each CV in the order they appear in the plumed file
    types : [str], str
        CV types
    timings : [float, int]
        steering schedule in ns
    values : [[float,str]]
        CV values in default PLUMED units. 'initial' will be replaced by a computed initial value
    forces : [[float,int]]
        forces to be applied to each CV in kJ/mol
    reference : str
        path to reference PDB file if using RMSD as a CV. All of the atoms in the reference have to also appear in the system, but not vice versa.
    engine : str
        MD engine to run sMD with. Allowed 'AMBER' and 'GROMACS'

    Returns
    -------
    process : BioSimSpace.Process
        BioSimSpace process used for steering
    """
    # load system
    system = BSS.IO.readMolecules([topology, coordinates])
    pytraj_system = pt.load(coordinates, top=topology) # using pytraj system for searching and some initial value computing

    # parse masks to atoms
    cv_atoms = __parse_masks(masks, pytraj_system, types, reference)

    #load reference
    if reference is not None:
        if isinstance(reference, str):
            reference = BSS.IO.readMolecules(reference).getMolecule(0)
        else:
            reference = [BSS.IO.readMolecules(ref).getMolecule(0) for ref in reference]

    # create CVs
    timings, values, forces = __validate_input(timings, values, forces)
    cvs = __create_cvs(system, types, cv_atoms, reference)

    # append values with initial values
    # at the moment torsion and distance CVs do not have
    # a 'getInitialValue()' option
    restraints = __compute_restraints(cvs, values, forces, pytraj_system)

    # create the protocol
    protocol = __create_protocol(cvs, timings, restraints)

    # create the process
    if engine == 'AMBER':
        process = BSS.Process.Amber(system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda') # specify using pmemd.cuda
    else:
        process = BSS.Process.createProcess(system, protocol, engine)

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
        copyfile(f'{process.workDir()}/{file}', to_copy[file])

    # dry trajectory
    get_dry_trajectory(topology, 'steering.nc', 'steering_dry.nc')

    return None
