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
$ ammo setup --input inputs/open-peptide-1.pdb --system peptide-1 --state open --slurm
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

BioSimSpace supports 3 types of CVs: distance, dihedral angle and RMSD. This covers all of the above, except of P185 stacking to W179. The stacking is defined as the absolute difference between the following two distances: P185(CG)-W179(CE) and P185(CA)-W179(CD1). This requires a custom expression. As such, the sMD simulations here were run manually. The input files can be found in `systems/peptide-1/open/seeded-md/steering` and `systems/peptide-1/closed/seeded-md/steering` to represent the steering protocol from open to closed conformation and from closed to open conformation respectively.

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

Additionally, data analysis found in the SI is also in the `analysis` folder.

## References

<a id="1">[1]</a> Wiesmann, C.; Barr1, K. J.; Kung, J.; Zhu, J.; Erlanson, D. A.; Shen, W.; Fahr, B. J.; Zhong, M.; Taylor, L.; Randal1, M.; McDowell1, R. S.; Hansen, S. K. Nat. Struct. Mol. Biol. 2004, 11, 730–737.

<a id="1">[2]</a> Brandão, T. A. S.; Hengge, A. C.; Johnson, S. J. J. Biol. Chem. 2010, 285, 15874–15883.

<a id="1">[3]</a> Choy, M. S.; Li, Y.; Machado, L. E.; Kunze, M. B.; Connors, C. R.; Wei, X.; Lindorff-Larsen, K.; Page, R.; Peti, W. Mol. Cell 2017, 65, 644–658.