# AMMo

AMMo (Allostery in Markov Models) is a collection of a python library, command line scripts and tools, for assessing allosteric modulators. The main use for easier Markov State Model building and comparison, using steered MD simulations to explore a larger conformational space. Depending on how much customization is desired, there are 3 levels of tools available:

<img src="data/structure.png" width=500>

The highest level is the set of command line tools that create and run an [AMMo project](examples/ammo_basic/run_from_command_line_ammo.md). The lowest level is the [ammo python library](examples/ammo_basic/run_from_notebook.ipynb), with some default functions to run and analyse MD simulations (if even more customization is needed, [BioSimSpace](https://biosimspace.openbiosim.org/) is recommended to set up own simulations). The [command line scrips](examples/ammo_basic/run_from_command_line_python.md) simply run the functions available in the python library from the command line, allowing for easier execution and job scheduling.

### Examples
In-depth explanations on how to use the [python library](examples/ammo_basic/run_from_notebook.ipynb), the command line [scripts](examples/ammo_basic/run_from_command_line_python.md), and the AMMo [project](examples/ammo_basic/run_from_command_line_ammo.md) are available in the [ammo_basic](examples/ammo_basic/) of examples folder. Example folder also includes two case studies on applying AMMo to Protein Tyrosine Phosphatase 1B ([PTP1B](examples/PTP1B)) and Exchange Proteins directly Activated by cAMP ([EPACs](examples/EPACs)).

### Installation

To use `ammo`, clone this repository and run `python setup.py`.

The Python interpreter used to run AMMo commands (default: /usr/bin/python3) and the location of allostery projects (default: ~/Documents) can be set by adding the arguments --interpreter and --location when running python `python setup.py`.

These settings can also be modified in the configuration file after setup.

### Requirements

* a python evnironment containing:
    * BioSimSpace
    * pytraj[^1]
* AMBER[^1]
* GROMACS

[^1]: If AMBER is compiled with pytraj compatible with the BioSimSpace environment, a separate installation of pytraj is not needed. If AMBER is compiled with pytraj that is not compatible with the BioSimSpace environment, it needs to be removed from `PYTHONPATH` and a separate pytraj installed in the environment.
