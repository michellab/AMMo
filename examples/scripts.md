# Allostery example

## Contents

1. [System setup](#System-setup)
2. [Running steered MD](#Running-steered-MD)
3. [Analysing steered MD data](#Analysing-steered-MD-data)
4. [Seeded MD](#Seeded-MD)
5. [Trajectory featurization](#Trajectory-featurization)

## System setup
[top](#Allostery-example)

A common starting point in MD simulations is a PDB file. The script `setup_system.py` runs `allostery.setup.setup_system()` from command line. The inputs are the same as for the function:

```python
$ python setup_system.py -h
usage: setup_system.py [-h] --input INPUT --protocol PROTOCOL [--engine ENGINE] [--charges CHARGES] [--parameters PARAMETERS]
                       [--topology TOPOLOGY] [--solvated]

General system setup: parameterisation, solvation, minimisation, heating, and equilibration.

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         System PDB file
  --protocol PROTOCOL   Comma separated list of minimisation steps, heating duration in ps, and equilibration duration in ps
  --engine ENGINE       Simulation engine supported by BioSimSpace. Default: "GROMACS"
  --charges CHARGES     Ligand charges in the order they appear in the input PDB, comma separated
  --parameters PARAMETERS
                        any additional parameter arguments to give LeAP separated by comma
  --topology TOPOLOGY   Dry topology of system. If provided will be used instead of re-parameterising
  --solvated            Indicate that the system is already solvated
```

If a `--topology` file is provided, the system will not be re-parameterised. Additionally, if `--solvated`, no additional preparation will be done, and the script will go straight to minimisation.

Running the script as:
```python
$ python setup_system.py --input input_protein.pdb --protocol "7500,100,250" --parameters "source leaprc.phosaa10"
```
produces system topology files (dry and solvated), and coordinate files for the minimised, heated, and equilibrated system.

## Running steered MD
[top](#Allostery-example)

Once the system is prepared, the next step is to run steered MD simulations. This allows for better sampling of intermediate conformations which are unstable and therefore short-lived. They can be run using the `steered_md.py` script:
```bash
$ python steered_md.py -h
usage: steered_md.py [-h] --topology TOPOLOGY --coordinates COORDINATES --masks MASKS --types TYPES --timings TIMINGS --values VALUES --forces FORCES [--reference REFERENCE]
                     [--engine ENGINE]

Run a steered MD simulation

optional arguments:
  -h, --help            show this help message and exit
  --topology TOPOLOGY   system topology
  --coordinates COORDINATES
                        equilibrated coordinates
  --masks MASKS         a list of AMBER masks for each CV
  --types TYPES         CV types
  --timings TIMINGS     steering schedule in ns
  --values VALUES       CV values in default PLUMED units. "initial" will be replaced by a computed initial value
  --forces FORCES       forces to be applied to each CV in kJ/mol
  --reference REFERENCE
                        path to reference PDB file(s) if using RMSD as a CV. All of the atoms in the reference have to also appear in the system
  --engine ENGINE       MD engine to run sMD with. Default: AMBER
```

The `--topology` and `--coordinate` parameters are self-explanatory. The `--masks` are AMBER selection masks, corresponding to the atoms involved in each CV used for steering. For example, the distance between the C 𝛼  atoms of residues 100 and 200 would be ":100@CA :200@CA". More information can be found [here](https://amberhub.chpc.utah.edu/atom-mask-selection-syntax/). `--types` corresponds to CV types supported by BioSimSpace.

#### Single steering step

Below is an example of multi-CV steering. The CVs will be the  𝜒 1 angle of Tyr152, the distance between C 𝛾  atoms of residues 196 and 280, and the heavy atom RMSD of residues 178-184. Note that in the masks below, the residue numbers are offset by 1. The system includes an ACE cap at the start, and the mask selection indexes starting from 1. The steering will be carried out in 100 ns. In addition to the specified values, times, and forces, additional steps will be added to apply the force over 4 ps, keeping the CV values as initial. The target values and forces used are based on knowledge of the system.

```bash
$ python steered_md.py --topology system.prm7 \
                       --coordinates system_equilibrated.rst7 \
                       --masks "[':153&(@N,CA,CB,CG)', ':197@CG :281@CG', ':179-185&!(@/H)']" \
                       --types "['torsion', 'distance', 'rmsd']" \
                       --timings "[100]" \
                       --values "[-1.047, 0.7, 0]" \
                       --forces "[2500,2500,2500]" \
                       --reference reference.pdb
```

After the steering process is finished, the output files are copied over from the working directory, and additionally a dry copy of the trajectory (no waters or ions) is saved (in case of further analysis required).

#### Multiple step steering

In order to specify multiple steering steps, `--timings`, `--values` and `--forces` need to be provided. For example, if we wanted to steer the dihedral angle during the first 50 ns of the simulation and the distance during the second, while steering the RMSD throughout, the input would look like this:
```bash
$ python steered_md.py --topology system.prm7 \
                       --coordinates system_equilibrated.rst7 \
                       --masks "[':153&(@N,CA,CB,CG)', ':197@CG :281@CG', ':179-185&!(@/H)']" \
                       --types "['torsion', 'distance', 'rmsd']" \
                       --timings "[50,100]" \
                       --values "[[-1.047, -1.047,], ['initial', 0.7], ['initial/2', 0]]" \
                       --forces "[[2500, 2500], [2500, 2500], [2500, 2500]]" \
                       --reference reference.pdb
```

Note that simple mathematical operations are allowed for the initial value, and this way the RMSD steering is not affected.

In this case the dihedral angle CV was steered to its target value and kept constant by applying force, while the distance CV was kept at its initial value by applying force during the first half of the simulation. An alternative protocol where they are not steered at all beyond changing the CV value could be employed by simply changing the appropriate force constants to 0.

Because a steered MD simulation can be quite complex if there are a lot of CVs and steps involved, running the script as above can involve a lot of typing, remembering values and masks. Therefore the functionality of creating an allostery project and running sMD using pre-configured settings can be useful in the long run, while the above direct running of the script is recommended mainly for protocol testing and troubleshooting.

## Analysing steered MD data
[top](#Allostery-example)

Once a steered MD trajectory is produced, it has to be checked to ensure steering has been successful, and snapshots need to be saved for seeded MD simulations. This simple trajectory analysis can be done however the user choses. There is an example notebook in `$ALLOSTERYHOME/data/sMD_analysis.ipynb` which can be a good starting point.

## Seeded MD
[top](#Allostery-example)

With snapshot saved from the sMD trajectory, they can be used as "seeds" to run equilibrium MD simulations. Since they are indeed just equilibrium MD simulations, they can be run using two scripts here: `equilibrium_md.py` and `seeded_md.py`.

```bash
$ python equilibrium_md.py -h
usage: equilibrium_md.py [-h] --duration DURATION [--topology TOPOLOGY]
                         [--coordinates COORDINATES] [--output OUTPUT]
                         [--report REPORT] [--workdir WORKDIR] [--clean]

Run equilibrium MD with BioSimSpace

optional arguments:
  -h, --help            show this help message and exit
  --duration DURATION   MD simulation duration in ns
  --topology TOPOLOGY   parameter file
  --coordinates COORDINATES
                        cooridnate file
  --output OUTPUT       output location, without extension (e.g.
                        "production-1"). If None, the next available name will
                        be used (e.g. "production-3" if "production-1" and
                        "production-2" already exist)
  --report REPORT       Report interval in steps. Default: 2500
  --workdir WORKDIR     Working directory
  --clean               Remove unneeded process files
```

`equilibrium_md.py` runs a simple MD simulation with AMBER from the command line. If `--coordinates` is a saved snapshot from the previous sMD trajectory, it runs seeded MD. Alternatively, `seeded_md.py` makes some assumptions on how data is structured:

```bash
$ python seeded_md.py -h
usage: seeded_md.py [-h] --folder FOLDER --snapshot SNAPSHOT --duration
                    DURATION [--report REPORT] [--clean]

Run a seeded MD simulation

optional arguments:
  -h, --help           show this help message and exit
  --folder FOLDER      Seeded MD folder, containing topology and a snapshots
                       directory
  --snapshot SNAPSHOT  Seed index
  --duration DURATION  Seeded MD duration in ns
  --report REPORT      Report interval in steps. Default: 2500
  --clean              Remove unneeded process files
```

It needs to be given a directory containing a `snapshots` folder, where the seed coordinate files are named `snapshot_[seed_idx].rst7`. When running the MD simulation, it will create a directory called `snapshot_[seed_idx]` and run MD there. For example:

```bash
$ python seeded_md.py --folder seeded_md_folder --snapshot 45 --duration 100 --report 5000 --clean
```

Will run a 100 ns MD simulation in the directory `seeded_md_folder/snapshot_45`, reporting every 10 ps and using coordinates from file `seeded_md_folder/snapshots/snapshot_45.rst7`.

Running seeded MD via a script (whether it be `equilibirum_md.py` or `seeded_md.py`) makes it easier to execute this part of the process using a job scheduler. Alternatively, it is very easy to write a simple python script to run MD differently with BioSimSpace, e.g. using a different MD engine.

## Trajectory featurization
[top](#Allostery-example)

Once seeded MD simulations are finished, they can be used to build a Markov State Model. However, that requires dimensionality reduction, which starts by reducing trajectory data from all atom coordinates to select features.

```bash
$ python featurize.py
usage: featurize.py [-h] --topology TOPOLOGY --trajectory TRAJECTORY --feature
                    FEATURE --mask MASK --output OUTPUT
                    [--reference REFERENCE] [--shared SHARED]

Reduce trajectory data to a single feature. Run in the folder containing
snapshot folders

optional arguments:
  -h, --help            show this help message and exit
  --topology TOPOLOGY   system topology file 
  --trajectory TRAJECTORY
                        seeded MD trajectory file
  --feature FEATURE     Type of feature to calculate. Allowed: "rmsd",
                        "torsion", "distance"
  --mask MASK           AMBER selection mask for the feature
  --output OUTPUT       output file for feature
  --reference REFERENCE
                        RMSD reference PDB file
  --shared SHARED       mask for atoms used for aligning the RMSD reference.
                        Default: "!@/H"
```

`featurize()` computes a distance, dihedral or RMSD values for the trajectory specified, using an AMBER selection mask ([documentation](https://amberhub.chpc.utah.edu/atom-mask-selection-syntax/)). If RMSD is being calculated, `reference` has to be provided as well, and `shared` is the selection mask for atoms used for alignment. For example:

```bash
$ python featurize.py --topology 'system.prm7' --trajectory 'production.nc' --feature 'distance', --mask ':10@CA :20@CA' --output 'featurized.txt'
```

A good idea would be to include this featurization as part of the seeded MD job above, if the MSM features are known beforehand.

## In progress
This tool is currently a work in progress. Functionality still to come is:
* MSM building