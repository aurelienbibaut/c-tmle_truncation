#!/bin/bash
# Job name:
#SBATCH --job-name=consistency_check.R
#
# Partition:
#SBATCH --partition=savio2
#
# Account:
#SBATCH --account=co_biostat
#
# Processors:
#SBATCH --ntasks=48
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=aurelien.bibaut@berkeley.edu

## Command(s) to run:
module load r
module load openmpi

#mpirun ./consistency_check.R
mpirun -n 1 R --vanilla < consistency_check.R > consistency_check.out
