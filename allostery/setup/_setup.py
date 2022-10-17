import BioSimSpace as BSS


def __get_molecule_types(molecules):
    """Get molecule types for parameterisation. It is assumed that the first molecule is the "main" protein and anything
     else with more than one residue is a "peptide". Single residue molecules are classed as ligands

    Parameters
    ----------
    molecules : [BioSimSpace._SireWrappers.Molecule]
        system molecules

    Returns
    -------
    types : [str]
        molecule types: "protein", "peptide", "ligand"
    """
    types = []
    for i, molecule in enumerate(molecules):
        if i==0:
            types.append('protein')
        elif molecule.nResidues()>1:
            types.append('peptide')
        else:
            types.append('ligand')

    return types


def __load_system(input_file, parameters=None, ligand_charges=None):
    """Loads a system from

    Parameters
    ----------
    input_file : str
        System PDB file

    parameters : [str]
        any additional parameter arguments to give LeAP

    ligand_charges : [int]
        list of ligand charges in the order they appear in the system

    Returns
    -------
    system : :class:`System <BioSimSpace._SireWrappers.System>`
        the combined molecular system
    """
    #load input
    input_system = BSS.IO.readMolecules(input_file)
    molecules = input_system.getMolecules()
    types = __get_molecule_types(molecules)

    #parameterise and combine into system
    system = BSS.Parameters.ff14SB(molecules[0], leap_commands=parameters).getMolecule().toSystem()
    ligand_count = 0
    for i in range(1,len(molecules)):
        if types[i] == 'peptide':
            parameterised = BSS.Parameters.ff14SB(molecules[i], leap_commands=parameters).getMolecule()
        elif types[i] == 'ligand':
            parameterised = BSS.Parameters.gaff2(molecules[i], net_charge=ligand_charges[ligand_count]).getMolecule()
            ligand_count+=1
        else:
            raise TypeError(f'Unsupported molecule type {type} for molecule with index {i}')
        system.addMolecules(parameterised)

    BSS.IO.saveMolecules('system_dry', system, 'prm7')

    return system


def __load_parameterised_system(coordinates, topology):
    """Load and saves a system that already has been parameterised.

    Parameters
    ----------
    coordinates : str
        System coordinate file

    topology : str
        system topology

    Returns
    -------
    system : :class:`System <BioSimSpace._SireWrappers.System>`
        the molecular system
    """
    system = BSS.IO.readMolecules([coordinates, topology])
    BSS.IO.saveMolecules('system_dry', system, 'prm7')

    return system


def __solvate(system, solvation='shell,10'):
    """Solvate in TIP3P water and save system.

    Parameters
    ----------
    system : :class:`System <BioSimSpace._SireWrappers.System>`
        the combined molecular system  

    solvation : bool
        way to solvate the system, e.g. 10 A shell is "shell,10", while a 15 A box in all dimensions is "box,15,15,15" (x, y, and z respectively)
        If None, the system will be treated as already solvated.

    Returns
    -------
    solvated_system : :class:`System <BioSimSpace._SireWrappers.System>`
        solvated system (TIP3P water)
    """
    if solvation is None:
        solvated_system = system
    else:
        solvation = solvation.split(',')
        # get solvation type, i.e. if box or shell
        box = solvation[0]
        # get dimensions
        dimensions = [float(solvation[i])*BSS.Units.Length.angstrom for i in range(1, len(solvation))]
        if box == 'shell':
            solvated_system = BSS.Solvent.tip3p(system, shell=dimensions[0], ion_conc=0.15, is_neutral=True)
        elif box == 'box':
            solvated_system = BSS.Solvent.tip3p(system, box=dimensions, ion_conc=0.15, is_neutral=True)
        else:
            raise ValueError(f'Solvation has to be either "shell" or "box", but is {box}') 

    BSS.IO.saveMolecules('system', solvated_system, ['prm7', 'rst7'])
    solvated_system = BSS.IO.readMolecules(['system.prm7', 'system.rst7'])

    return solvated_system


def __minimise(system, steps, engine):
    """Minimise PTP1B system.

    Parameters
    ----------
    system : :class:`System <BioSimSpace._SireWrappers.System>`
        the molecular system

    steps : int
        minimisation steps

    engine : str
        simulation engine.

    Returns
    -------
    minimised : :class:`System <BioSimSpace._SireWrappers.System>`
        the minimised system
    """

    protocol = BSS.Protocol.Minimisation(steps=steps)
    process = BSS.Process.createProcess(system, protocol, engine)

    process.start()
    process.wait()

    #save results
    minimised = process.getSystem()
    BSS.IO.saveMolecules('system_minimised', minimised, 'rst7')
    with open('min.out', 'w') as file:
        file.writelines([line+'\n' for line in process.getStdout()])

    return minimised


def __heat(system, runtime, engine):
    """Heat system from 0 K to 300 K over a specified time.

    Parameters
    ----------
    system : :class:`System <BioSimSpace._SireWrappers.System>`
        the molecular system

    runtime : :class:`Time <BioSimSpace.Types.Time>`
        time over which system is heated

    engine : str
        simulation engine.

    Returns
    -------
    heated : :class:`System <BioSimSpace._SireWrappers.System>`
        the heated system
    """
    K = BSS.Units.Temperature.kelvin
    protocol = BSS.Protocol.Equilibration(runtime=runtime, temperature_start=0*K, temperature_end=300*K)
    process = BSS.Process.createProcess(system, protocol, engine)

    process.start()
    process.wait()

    #save results
    heated = process.getSystem()
    BSS.IO.saveMolecules('system_heated', heated, 'rst7')
    with open('heat.out', 'w') as file:
        file.writelines([line+'\n' for line in process.getStdout()])

    return heated


def __equilibrate(system, runtime, engine):
    """Equilibrate system over a specified time.

    Parameters
    ----------
    system : :class:`System <BioSimSpace._SireWrappers.System>`
        the molecular system

    runtime : :class:`Time <BioSimSpace.Types.Time>`
        time over which system is equilibrated

    engine : str
        simulation engine.

    Returns
    -------
    None
    """
    protocol = BSS.Protocol.Equilibration(runtime=runtime, temperature=300*BSS.Units.Temperature.kelvin,
                                          pressure=1*BSS.Units.Pressure.atm)
    process = BSS.Process.createProcess(system, protocol, engine)

    process.start()
    process.wait()

    equilibrated = process.getSystem()
    BSS.IO.saveMolecules('system_equilibrated', equilibrated, 'rst7')
    with open('equil.out', 'w') as file:
        file.writelines([line+'\n' for line in process.getStdout()])

    return equilibrated


def setup_system(input_file, protocol, engine='GROMACS', ligand_charges=None, parameters=None, topology=None, solvation='shell,10'):
    """General system setup protocol. The input system will be minimised, heated in the NVT ensemble and then equilibrated in the NPT ensemble.
    The system can be provided as a PDB, in which case it will be parameterised with ff14SB. The first molecule is assumed to be the protein, and the following molecules that have more than one residue are treated as "peptide". Single residue molecules are treated as "ligand" and parameterised with gaff2.
    The system is also solvated in TIP3P water (0.15 mM NaCl conc) with a 10 Angstrom shell, unless specified as already solvated.

    Each step will produce a saved system coordinate file in the current directory, as well as the process output. If not "solvated=True" the dry system topology will also be saved.

    Parameters
    ----------
    input_file : str
        System coordinate file

    protocol : tuple
        A tuple consisting of number of minimisation steps, heating duration in ps and equilibration duration in ps

    engine : str
        simulation engine supported by BioSimSpace. Default: 'GROMACS'

    ligand_charges : int, [int]
        list of ligand charges in the order they appear in the system

    parameters : [str]
        any additional parameter arguments to give LeAP when parameterising proteins

    topology : str
        topology of the system. If provided, the system will not be parameterised

    solvation : str, None
        way to solvate the system, e.g. 10 A shell is "shell,10", while a 15 A box in all dimensions is "box,15,15,15" (x, y, and z respectively)
        If None, the system will be treated as already solvated.

    Returns
    -------
    equilibrated :
        equilibrated system
    """
    if isinstance(ligand_charges, int):
        ligand_charges = [ligand_charges]
    if topology is None:
        system = __load_system(input_file, parameters, ligand_charges)
        print('Parameterising system...', end = '')
    else:
        system = __load_parameterised_system(input_file, topology)
        print('Loading parameterised system...', end = '')
    
    solvated_system = __solvate(system, solvation)
    print('done.\n------------------------------\nMinimising system...', end='')

    minimised = __minimise(solvated_system, protocol[0], engine)
    print('done.\n------------------------------\nHeating system...', end='')
    heated = __heat(minimised, protocol[1]*BSS.Units.Time.picosecond, engine)
    print('done.\n------------------------------\nEquilibrating system...', end='')
    equilibrated = __equilibrate(heated, protocol[2]*BSS.Units.Time.picosecond, engine)
    print('done.\n------------------------------')
    print('System setup completed successfully.\n------------------------------')

    return equilibrated
