import os

import BioSimSpace as BSS
import pytraj as pt
from shutil import copyfile

from allostery.steering import PlumedFile, RMSDReference


def __parse_masks(masks, system):
    """
    Parse masks to return atom indices for the system

    Parameters
    ----------
    masks : [str], str
        a list of AMBER masks for each CV in the order they appear in the plumed file
    system : pytraj.Trajectory
        the system to be steered

    Returns
    -------
    atoms : [[int]]
        a list of atom lists for each CV
    """
    if type(masks) is str:
        masks = [masks]

    atoms = []
    for mask in masks:
        results = system.topology.select(mask)
        atoms.append(results)

    return atoms


def __compute_value(system, cv_type, atoms):
    """Compute CV values using pytraj"""
    if cv_type == 'DISTANCE':
        result = pt.calc_distance(system, atoms)[0][0]/10
    elif cv_type == 'TORSION':
        result = pt.calc_dihedral(system, atoms)[0][0]*0.017453
    else:
        raise TypeError(f'Unsupported CV: {cv_type}')
    return result


def __compute_initial_values(types, atoms, system, work_dir, reference=None):
    """
    Compute the initial values of CVs

    Parameters
    ----------
    types : [str]
        CV types
    atoms : [[int]]
        list of atom lists for each CV
    system : pytraj.Trajectory
        pytraj system
    work_dir : str
        path to the BioSimSpace Process working directory
    reference : BioSimSpace.Molecule
        RMSD reference

    Returns
    -------
    initial_values : [float]
        starting CV values
    """
    initial_values = []
    rmsd_index = 1
    for i in range(len(types)):
        # if RMSD, create BSS CV
        if types[i] == 'RMSD':
            rmsd_reference = RMSDReference(reference)
            rmsd_reference.renumber_atoms(system)
            initial_values.append(rmsd_reference.get_initial_value(system, atoms[i]+1, f'reference_{rmsd_index}.pdb', work_dir))
            rmsd_index += 1
        else:
            initial_values.append(__compute_value(system, types[i], atoms[i]))
    return initial_values


def __edit_plumed(plumed, atoms, initial_values):
    """
    Edit the PLUMED input file with correct atoms, RMSD references and
    Parameters
    ----------
    plumed : allostery.steering.PlumedFile
        provided PLUMED input file
    atoms : [[int]]
        a list of atom lists for each CV, in the order they appear in the PLUMED file
    initial_values : [float]
        starting CV values

    Returns
    -------
    None
    """
    # change atoms first
    for i in range(len(atoms)):
        entry = plumed.variables[i]
        # ignore RMSD entries as their reference was written before
        if entry.type != 'RMSD':
            entry_atoms = [str(atom+1) for atom in atoms[i]]
            entry.atoms = ','.join(entry_atoms)
    # change initial values in steps 0 and 1
    initial_values = [str(value) for value in initial_values]
    for i, step in enumerate(plumed.restraint):
        if 'initial' in step:
            parts = step.split()
            step_at = parts[1].split('=')[0]
            step_values = parts[1].split('=')[1].split(',')
            for j, value in enumerate(step_values):
                if value=='initial':
                    step_values[j] = initial_values[j]
            new_step = f'{parts[0]} {step_at}={",".join(step_values)} {parts[2]}\n'
            plumed.restraint[i] = new_step
    return None


def run_smd(topology, coordinates, masks, plumed, reference=None, engine='AMBER'):
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
    plumed : str
        plumed file with steering protocol, where "ATOMS=" parameter is left empty for all CVs
    reference : str
        path to reference PDB file if using RMSD as a CV. All of the atoms in the reference have to also appear in the system, but not vice versa
    engine : str
        MD engine to run sMD with. Allowed 'AMBER' and 'GROMACS'


    Returns
    -------
    process : BioSimSpace.Process
        BioSimSpace process used for steering
    """
    # load system
    system = BSS.IO.readMolecules([topology, coordinates])
    pytraj_system = pt.load(coordinates, top=topology) # using pytraj system for searching since the AMBER mask selection is much better

    # parse masks to atoms
    cv_atoms = __parse_masks(masks, pytraj_system)

    # load in the provided PLUMED file
    plumed_file = PlumedFile(plumed)
    # and get total simulation time
    simulation_time = plumed_file.get_total_time()*BSS.Units.Time.nanosecond

    # create a dummy CV, protocol and process so that files could be copied there
    tmp_cv = BSS.Metadynamics.CollectiveVariable.Distance(1, 2) # system will have at least 2 atoms
    protocol = BSS.Protocol.Steering(tmp_cv, [0 * BSS.Units.Time.nanosecond, simulation_time],
                                     [BSS.Metadynamics.Restraint(1*BSS.Units.Length.nanometer, 0)]*2,
                                     runtime=simulation_time,
                                     report_interval=2500, restart_interval=2500)
    if engine == 'AMBER':
        process = BSS.Process.Amber(system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda')
    else:
        process = BSS.Process.createProcess(system, protocol, engine)

    # compute initial values
    types = plumed_file.get_types()
    initial_values = __compute_initial_values(types, cv_atoms, pytraj_system, process.workDir(), reference)

    # edit the PLUMED file
    __edit_plumed(plumed_file, cv_atoms, initial_values)
    # change process PLUMED file
    plumed_file.get_file(f'{process.workDir()}/plumed.dat')

    # run process
    process.start()
    process.wait()
    to_copy = {'amber.nc': 'steering.nc', 'plumed.dat': 'plumed.in', 'steering.dat': 'steering.dat', 'reference_1.pdb': 'reference_1.pdb', 'reference_2.pdb': 'reference_2.pdb'}
    for file in to_copy:
        copyfile(f'{process.workDir()}/{file}', to_copy[file])

    return None
