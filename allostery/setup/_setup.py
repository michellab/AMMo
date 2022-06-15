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
            parameterised = BSS.Parameters.gaff(molecules[i], net_charge=ligand_charges[ligand_count]).getMolecule()
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


def __solvate(system):
    """Solvate and save dry system.

    Parameters
    ----------
    system : :class:`System <BioSimSpace._SireWrappers.System>`
        the combined molecular system

    Returns
    -------
    solvated : :class:`System <BioSimSpace._SireWrappers.System>`
        solvated system (TIP3P water)
    """
    solvated = BSS.Solvent.tip3p(system, shell=10*BSS.Units.Length.angstrom, ion_conc=0.15, is_neutral=True)
    BSS.IO.saveMolecules('system', solvated, ['prm7', 'rst7'])
    solvated = BSS.IO.readMolecules(['system.prm7', 'system.rst7'])

    return solvated


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


def setup_system(input_file, protocol, engine='GROMACS', ligand_charges=None, parameters=None, topology=None):
    """General system setup protocol, intended for PTP1B.

    Parameters
    ----------
    input_file : str
        System PDB file

    protocol : tuple
        A tuple consisting of number of minimisation steps, heating duration in ps and equilibration duration in ps

    engine : str
        simulation engine. Can be either 'AMBER' or 'GROMACS'. Default = 'GROMACS'

    ligand_charges : int, [int]
        list of ligand charges in the order they appear in the system

    parameters : [str]
        any additional parameter arguments to give LeAP

    topology : str
        AMBER topology of the system. If provided, the dry system will not be parameterised

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
    solvated = __solvate(system)
    print('done.\n------------------------------\nMinimising system...', end='')

    minimised = __minimise(solvated, protocol[0], engine)
    print('done.\n------------------------------\nHeating system...', end='')
    heated = __heat(minimised, protocol[1]*BSS.Units.Time.picosecond, engine)
    print('done.\n------------------------------\nEquilibrating system...', end='')
    equilibrated = __equilibrate(heated, protocol[2]*BSS.Units.Time.picosecond, engine)
    print('done.\n------------------------------')
    print('System setup completed successfully.\n------------------------------')

    return equilibrated
