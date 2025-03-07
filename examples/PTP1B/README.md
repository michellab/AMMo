# Example - PTP1B

An example application of AMMo to a protein system - Protein Tyrosine Phosphatase 1B (PTP1B).

Inhibition of PTP1B has been proposed as a therapeutic strategy for for type II diabetes treatments[[1]](#1). The protein enzymatic activity is regulated by the conformations adopted by the WPD loop that sits above the active site. The WPD loop adopts two major conformations - open and closed (Figure 1a). The closed conformation is catalytically active, as it positions the catalytic residues in range of the substrate[[2]](#2). The closing and opening of the loop occurs on multimicrosecond timescales[[3]](#3). Three experimentally characterised inhibitors (1-3) and a fragment binder (4) with unknown functional effect were used to validate the methodology. They are shown in Figure 1b.

## Contents



## Preparation
[top](#Example---PTP1B)

The first step is to create a project folder, where all of the data will be located. AMMo commands in the CLI start with `ammo`. To see all available commands and their brief decriptions, run:
```bash
$ ammo -h
```

To create a project folder, use `ammo project`:
```bash
$ ammo project --create PTP1B
```
This will open a settings file, where project defaults are set. The settings used in this case are given in .defaults/config (where they are stored for every project). To change these settings, simply run `ammo project --reconfigure`. Note that this will reconfigure the current active project. To see which project is active, use `ammo project --current`, and to activate a different project use `ammo project --activate [name]`. For the rest if this example, it is assumed that commands are issued from the project folder.

## Setup
[top](#Example---PTP1B)

Each variation on the protein system (different ligands, no peptide substrate, etc.) is a different system, all of which are found in the `systems` folder. To add a new system, use `ammo systems --create [name]`. Alternatively, the folder can be created manually, but the command will create all required subdirectories. For example:
```bash
$ ammo systems --create peptide-1
```

The starting point to use AMMo is a clean, prepared PDB file (files used here are found in `inputs`, where it is recommended to keep them for every project). We will use the `peptide-1` system, `open` state (i.e. PTP1B with the peptide substrate and ligand **1**) as an example here. The system folder has already been created above. To set it up, use the command:
```bash
$ ammo setup --input open-peptide-1.pdb --system peptide-1 --state open --slurm
```

The above writes a `submit.sh` file and submits it. In case of errors, the `submit.sh` and `submit.out` files are good starting points for troubleshooting. The setup produces `system.prm7` and `system_equilibrated.rst7` files, which are the system topology and equilibrated coordinates respectively. If the system was setup outside of using AMMo, the topology and coordinates can simply be placed in the appropriate `system-setup` folder and named appropriately. The `system_dry.prm7` file is the system topology excluding water and ions, for later processing.

If a topology file is provided, the system is not parameterised (as was in the case of the system including covalently linked ligand 3, which required custom parameters). Additionally, if `solvation` is `None` in the settings, it is assumed that the provided system has already been solvated.

## Steered MD
[top](#Example---PTP1B)

Steered MD is run using the `steering` command, and the main options are provided in the settings file. In this case, the steering parameters were:
* WPD loop RMSD
* Y152 χ1 angle
* P185 stacking to W179
* F196 stacking to F280 (CG distances)

AMMo uses PLUMED[[4]](#4) to apply a harmonic restraint to the system and so bias it towards a certain value of a collective variable (CV) (or multiple CVs). While BioSimSpace has support for steering with distance, torsion and RMSD based CVs[[5]](#5), in order to support larger CV flexibility AMMo uses a pseudo PLUMED input file. An example is shown below:

```bash
rmsd: RMSD REFERENCE=:179-185&(!@/H) TYPE=OPTIMAL FILE=closed_ref.pdb
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
PRINT STRIDE=2500 ARG=* FILE=COLVAR
```

In a usual PLUMED file, the `ATOMS=` arguments would point to atom indices, while here they are replaced with AMBER atom masks. This makes the pseudo-PLUMED files easier to write, and makes them more transferable. When sMD is being set up, the masks are replaced with the corresponding atom indices. In case of RMSD, the mask is used to prepare the reference file, specified by the `FILE=` argument, which will be removed in the final PLUMED input file. Additionally, the `initial` keyword is allowed when defining the steering steps. During set up it is replaced with the starting value of the CV. Simple arithmetic operations, such as `initial/2` or `initial+5.2` are also allowed.

The pseudo-PLUMED files used here are available in `.defaults`, and the references for RMSD are in `inputs`. When setting up the reference files, first the path will be searched as give, and then `inputs` directory will be searched. The `steering` part of `.defaults/config` points to these files for AMMo to use.

Once the sMD trajectories were obtained, the analysis was run in `analysis/sMD_analysis.ipynb`.

## Seeded MD
[top](#Example---PTP1B)

Each snapshot saved from sMD was resolvated, minimized and equilibrated. Then 50 ns of MD simulation was run. This was done using the `seeded` command as follows:

```bash
$ ammo seeded --system peptide-1 --state open
```

This takes the sample submission script (`.defaults/seeded-md.sh`), copies it to a remote HPC cluster along with any required files, and submits it after making some changes (such as adding the appropriate working directory and `rsync` commands for backup and copying data back to the local workstation). In this case, `AMMo` was setup on the remote cluster as well, but this is not necessary if appropriate input files for MD simulations are prepared and used in the submission script. This feature does require the setup of ssh keys between the remote cluster and the local workstation to avoid password input.

## Featurization
[top](#Example---PTP1B)

Each seeded MD trajectory was reduced to 2 features: the RMSD of the WPD loop backbone atoms, and the RMSD of the P loop backbone atoms, both using the closed WPD loop conformation structure as reference (PDB ID: 1SUG, `inputs/closed.pdb`).

```bash
$ ammo featurize --system peptide-1 --state open --slurm
```

The parameters for featurization can be found in the `features` part of the settings file.

## Markov State Modelling
[top](#Example---PTP1B)

Markov State Models were built in a jupyter notebook, which can be found in `analysis/msm.ipynb`. It makes use of the `MSMCollection` part of the `ammo` python library, which allows to build mutiple MSMs at once, so ensuring consistent clustering and metastable state assignments. More details can be found in the notebook itself.

Additionally, data analysis found in the SI is in `analysis/supplementary_information.ipynb`.

## References

<a id="1">[1]</a> Wiesmann, C.; Barr1, K. J.; Kung, J.; Zhu, J.; Erlanson, D. A.; Shen, W.; Fahr, B. J.; Zhong, M.; Taylor, L.; Randal1, M.; McDowell1, R. S.; Hansen, S. K. Nat. Struct. Mol. Biol. 2004, 11, 730–737.

<a id="2">[2]</a> Brandão, T. A. S.; Hengge, A. C.; Johnson, S. J. J. Biol. Chem. 2010, 285, 15874–15883.

<a id="3">[3]</a> Choy, M. S.; Li, Y.; Machado, L. E.; Kunze, M. B.; Connors, C. R.; Wei, X.; Lindorff-Larsen, K.; Page, R.; Peti, W. Mol. Cell 2017, 65, 644–658.

<a id="4">[4]</a>https://www.plumed.org/doc-v2.7/user-doc/html/_m_o_v_i_n_g_r_e_s_t_r_a_i_n_t.html

<a id="5">[5]</a>https://github.com/michellab/BioSimSpaceTutorials/tree/main/03_steered_md
