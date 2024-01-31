#!/bin/sh
#$ -N seededMD
#$ -pe gpu-titanx 1
#$ -l h_rt=23:00:00 
#$ -l h_vmem=22G

source software/amber22/amber.sh
source software/AMMo/ammo.sh
module load cuda/10.1.105
module load phys/compilers/gcc/8.4.0
module load anaconda
source activate BSS-env

python software/AMMo/ammo/scripts/seeded_md.py --folder "." --snapshot $SGE_TASK_ID --duration 50 --report 5000 --clean
