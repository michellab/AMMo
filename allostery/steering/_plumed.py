import os
import subprocess
import numpy as np
import pytraj as pt

plumed_functions = (
    'COMBINE', 'CUSTOM', 'ENSEMBLE', 'FUNCPATHMSD', 'FUNCSUMHILLS', 'LOCALENSEMBLE', 'MATHEVAL', 'PIECEWISE', 'SORT',
    'STATS')


class PlumedFile:
    """A class for managing PLUMED input files for MOVINGRESTRAINT"""

    def __init__(self, file=None):
        """
        Parameters
        ----------
        file : str
            existing PLUMED input file
        """
        self.variables = []  # variables such as distance, rmsd, with titles as
        self.functions = []  # any functions
        self.restraint = []  # steering protocol
        self.print = "PRINT ARG=* FILE=COLVAR\n"  # by default print everything to file COLVAR

        if file is not None:
            self.parse_file(file)

    def parse_file(self, file):
        """Parse an existing PLUMED input file"""
        with open(file, 'r') as fl:
            contents = fl.readlines()

        restraint = False
        for line in contents:
            # determine whether started restraint section
            if 'MOVINGRESTRAINT' in line:
                if restraint:
                    restraint = False
                else:
                    restraint = True
                continue
            # if reached print line, parse and stop
            elif line.startswith('PRINT'):
                self.print = line
                break

            # if not restraint yet
            if not restraint:
                entry = PlumedEntry(line)
                if entry.type in plumed_functions:
                    self.functions.append(entry)
                # else add to variables
                else:
                    self.variables.append(entry)
            # if restraint, parse restraint
            else:
                self.restraint.append(line)
        return None

    def change_variable(self, label, argument):
        """
        Change a variable in the PLUMED file
        Parameters
        ----------
        label : str
            variable label
        argument : (str, str)
            argument to change, e.g. ("ATOMS", "10-20") to change the atom entry to ATOMS=10-20

        Returns
        -------
        None
        """
        for entry in self.variables:
            if entry.label == label:
                if argument[0] == 'type':
                    entry.type = argument[1]
                elif argument[0] == 'ATOMS':
                    entry.atoms = argument[1]
                else:
                    entry.args[argument[0]] = argument[1]
        return None

    def get_total_time(self, timestep=0.000002):
        """Get the total number of steps in the steering protocol
        Parameters
        ----------
        timestep : float
            timestep in ns. Default : 2 fs

        Returns
        -------
        time : float
            total simulation time in ns
        """
        final_step = self.restraint[-1]
        steps = int(final_step.split()[0].split('=')[1])  # number of steps
        time = steps * timestep  # in ns
        return time

    def get_types(self):
        """Get all variable types"""
        types = []
        for variable in self.variables:
            types.append(variable.type)
        return types

    def get_file(self, path=None):
        """Get the PLUMED file contents. If path is specified, write out to path"""
        contents = []
        for entry in self.variables:
            contents.append(entry.get_entry())
        for entry in self.functions:
            contents.append(entry.get_entry())
        contents.append('MOVINGRESTRAINT ...\n')
        for step in self.restraint:
            contents.append(f'{step}')
        contents.append('... MOVINGRESTRAINT\n')
        contents.append(self.print)
        contents.append('\n')

        if path is not None:
            with open(path, 'w') as fl:
                fl.writelines(contents)
        return contents


class PlumedEntry:
    """A class for single PLUMED entries, such as CVs and functions"""

    def __init__(self, entry=None):
        """
        Parameters
        ----------
        entry : str
            an entry line from a PLUMED input file
        """
        self.label = None
        self.type = None
        self.atoms = None
        self.args = {}

        if entry is not None:
            self.parse_entry(entry)

    def parse_entry(self, entry):
        """
        Parse an entry line
        Parameters
        ----------
        entry : str
            an entry line from a PLUMED input file
        """
        # remove new line character
        if entry[-1] == '\n':
            entry = entry[:-1]

        # split into arguments
        args = entry.split()

        # get label if not an argument
        # then first entry in form label:
        if 'LABEL' not in entry:
            self.label = args[0][:-1]
            args = args[1:]

        # get type
        self.type = args[0]

        # sort the other arguments
        for arg in args[1:]:
            arg_parts = arg.split('=')
            if arg_parts[0] == 'LABEL':
                self.label = arg_parts[1]
            elif arg_parts[0] == 'ATOMS':
                self.atoms = arg_parts[1]
            else:
                self.args[arg_parts[0]] = arg_parts[1]

    def get_entry(self):
        """

        Returns
        -------
        entry : str
            entry put together into a line for a PLUMED file
        """
        entry = f'{self.label}: {self.type}'
        if self.atoms is not None:
            entry = f'{entry} ATOMS={self.atoms}'
        for key, value in self.args.items():
            entry = f'{entry} {key}={value}'
        entry = f'{entry}\n'

        return entry


class RMSDReference:
    """A class for managing PLUMED RMSD references"""

    def __init__(self, file):
        self.chains = []
        self.residues = []
        self.atoms = []
        with open(file, 'r') as fl:
            contents = fl.readlines()

        chain = RMSDChain(None)
        residue = RMSDResidue(None, None, None)
        for line in contents:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # parse chain
                chain_id = line[21]
                if chain_id != chain.id:
                    chain = RMSDChain(chain_id)
                    self.chains.append(chain)
                # parse residue
                res_idx = int(line[22:26].strip())
                res_name = line[17:21]
                if res_idx != residue.idx:
                    residue = RMSDResidue(res_idx, res_name, chain)
                    chain.residues.append(residue)
                    self.residues.append(residue)
                # parse atom
                at_type = line[:6].strip()
                at_idx = int(line[6:11].strip())
                at_name = line[12:16]
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                element = line[76:78].strip()
                atom = RMSDAtom(at_idx, at_name, x, y, z, element, at_type, residue)
                residue.atoms.append(atom)
                self.atoms.append(atom)

    def renumber_atoms(self, system):
        """Renumber reference atoms based on system"""
        for atom in self.atoms:
            atom_mask = f':{atom.residue.idx}@{atom.name}'
            system_idx = system.top.select(atom_mask)[0] + 1
            atom.idx = system_idx
        self.atoms.sort(key=lambda x: x.idx)

    def get_initial_value(self, system, rmsd_indices, output='reference.pdb', work_dir='.'):
        """
        Get initial RMSD value using PLUMED driver

        Parameters
        ----------
        system : str, pytraj.Trajectory
            path to system PDB or a pytraj trajectory
        rmsd_indices : [int]
            list of atom indices to be used for RMSD
        output : str
            output to save reference to
        work_dir : str
            path to working directory

        Returns
        -------
        initial_value : float
            initial RMSD value in nm
        """
        if type(system) == pt.Trajectory:
            trajectory = f'{work_dir}/system.pdb'
            pt.save(trajectory, system, overwrite=True)
        else:
            trajectory = f'{work_dir}/{system}'

        self.write_file(rmsd_indices, f'{work_dir}/{output}')

        plumed_file = f'{work_dir}/initial_value.dat'
        with open(plumed_file, 'w') as file:
            file.writelines([f'rmsd: RMSD reference={work_dir}/{output} TYPE=OPTIMAL\n',
                             f'PRINT ARG=rmsd FILE={work_dir}/OUTPUT\n'])
        subprocess.run(['plumed', 'driver', '--mf_pdb', trajectory, '--plumed', plumed_file])

        initial_value = np.loadtxt(f'{work_dir}/OUTPUT')[1]
        for file in [plumed_file, f'{work_dir}/OUTPUT']:
            os.remove(file)
        return initial_value

    def write_file(self, rmsd_indices=(), output=None):
        """
        Write out reference file

        Parameters
        ----------
        rmsd_indices : [int]
            atom indices to be used for RMSD calculation
        output : str
            path to reference output

        Returns
        -------
        contents : [str]
            contents of reference PDB file
        """
        contents = []
        for atom in self.atoms:
            if atom.idx in rmsd_indices:
                occbeta = '  0.00  1.00'
            else:
                occbeta = '  1.00  0.00'
            contents.append(
                f'{atom.type:<6}{atom.idx:>5} {atom.name} {atom.residue.name}{atom.residue.chain.id}{atom.residue.idx:>4}    {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{occbeta}          {atom.element:>2}\n')

        if output is not None:
            with open(output, 'w') as file:
                file.writelines(contents)
        return contents


class RMSDChain:
    """RMSD reference chain"""
    def __init__(self, chain_id):
        self.id = chain_id
        self.residues = []


class RMSDResidue:
    """RMSD reference residue"""
    def __init__(self, idx, name, chain):
        self.idx = idx
        self.name = name
        self.atoms = []
        self.chain = chain


class RMSDAtom:
    """RMSD reference atom"""
    def __init__(self, idx, name, x, y, z, element, at_type, residue):
        self.idx = idx
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.element = element
        self.type = at_type
        self.residue = residue
