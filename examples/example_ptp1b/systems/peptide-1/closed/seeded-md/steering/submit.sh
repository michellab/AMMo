#!/bin/bash

#SBATCH -J sMD
#SBATCH --output=submit.out
#SBATCH --gres=gpu:1

pmemd.cuda -O -i amber.cfg -o steering.out -p ../../system-setup/system.prm7 -c ../../system-setup/system_equilibrated.rst7 -r steering.rst7 -x steering.nc

