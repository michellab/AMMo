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
  --solvation SOLVATION
                        Way to solvate the system, e.g. 10 A shell is "shell,10", while a 15 A box in all
                        dimensions is "box,15,15,15" (x, y, and z respectively). If None, the system will be
                        treated as already solvated.
```

If a `--topology` file is provided, the system will not be re-parameterised. Additionally, if `--solvation None`, no additional preparation will be done, and the script will go straight to minimisation. Otherwise, solvation can specify the size of box or shell to solvate the system in.

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
usage: steered_md.py [-h] --topology TOPOLOGY --coordinates COORDINATES --input INPUT [--engine ENGINE] [--workdir WORKDIR] [--suffix SUFFIX] [--restraint RESTRAINT]

Run a steered MD simulation

optional arguments:
  -h, --help            show this help message and exit
  --topology TOPOLOGY   system topology
  --coordinates COORDINATES
                        equilibrated coordinates
  --input INPUT         path to pseudo PLUMED file
  --engine ENGINE       MD engine to run sMD with. Default: AMBER
  --workdir WORKDIR     Working directory. If None, sMD will be run in a temporary folder and copied over. Default: "."
  --suffix SUFFIX       A suffix to add to the output of the simulation, e.g. "steering_1.nc" or "steering_1.dat". For cases when the steering is done in more than one step. If None, nothing will be added
  --restraint RESTRAINT
                        A pseudo flat bottom restraint file that will be used during the steering (currently only available for AMBER). Instead of atom indices, AMBER masks are used

```

The `--topology` and `--coordinate` parameters are self-explanatory. The `input` requires as pseudo PLUMED file, containing all the required information for steering, except the specific `ATOM` indices are replaced with [AMBER selection masks](https://amberhub.chpc.utah.edu/atom-mask-selection-syntax/), and in case of RMSD collective variables, an additional reference `FILE` parameter is added, which will be removed during PLUMED input preparation. An example pseudo PLUMED input file is given in `example_data/plumed_input.dat` (as well as in the specific use case examples), and more information can be found on the [PLUMED website](https://www.plumed.org/doc-v2.8/user-doc/html/_m_o_v_i_n_g_r_e_s_t_r_a_i_n_t.html).

#### Single steering step

Below is an example of multi-CV steering. The CVs will be the heavy atom RMSD of residues 178-184, the 𝜒1 angle of Tyr152, the stacking of residues 185 and 179, and the distance between C𝛾 atoms of residues 196 and 280. Note that in the masks below, the residue numbers are offset by 1. The system includes an ACE cap at the start, and the mask selection indices starting from 1. The steering will be carried out in 100 ns. The target values and forces used are based on knowledge of the system. 

```bash
$ cat example_data/plumed_input.dat

rmsd: RMSD REFERENCE=:179-185&(!@/H) TYPE=OPTIMAL FILE=example_data/reference.pdb
tyr: TORSION ATOMS=:153@N:153@CA:153@CB:153@CG
pro1: DISTANCE ATOMS=:180@CE2:186@CD
pro2: DISTANCE ATOMS=:180@CD1:185@CA
stacking: CUSTOM ARG=pro1,pro2 FUNC=abs(x-y) PERIODIC=NO
phe: DISTANCE ATOMS=:197@CG:281@CG
MOVINGRESTRAINT ...
  ARG=rmsd,tyr,stacking,phe
  STEP0=0    AT0=initial,initial,initial,initial    KAPPA0=0.0,0.0,0.0,0.0
  STEP1=2000    AT1=initial,initial,initial,initial    KAPPA1=3500.0,3500.0,3500.0,3500.0
  STEP2=75000000    AT2=0.0,1.047,0.0,0.45    KAPPA2=3500.0,3500.0,3500.0,3500.0
  STEP3=76000000    AT3=0.0,1.047,0.0,0.45    KAPPA3=0.0,0.0,0.0,0.0
... MOVINGRESTRAINT
PRINT STRIDE=2500 ARG=* FILE=steering.dat
```

The `"initial"` values for the CVs at steps 0 and 1 will be computed using PLUMED and filled in during final file setup. This, together with the use of AMBER atom masks, allows for easier steering preparation while still using the whole range of CVs in PLUMED.

```bash
$ python steered_md.py --topology system.prm7 \
                       --coordinates system_equilibrated.rst7 \
                       --input example_data/plumed_input.dat
```

After the steering process is finished, the output files are copied over from the working directory, and additionally a dry copy of the trajectory (no waters or ions) is saved (in case of further analysis required).

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

## MSM building

An example of Markov State Modelling is available in a [notebook](msm.ipynb), as the `MSMCollection` functionality is only available as part of the `ammo` python library.