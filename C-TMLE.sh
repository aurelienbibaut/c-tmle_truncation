#!/bin/bash
# Job name:
#SBATCH --job-name=C-TMLE-multi_orders
#
# Partition:
#SBATCH --partition=savio2
#
# Account:
#SBATCH --account=co_biostat
# Processors:
# SBATCH --ntasks=48
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Mail type:
# SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=aurelien.bibaut@berkeley.edu
#

module load r openmpi Rmpi
mpirun R CMD BATCH --no-save C-TMLE_truncation-multi_orders-mpi-loglik_targeting.R C-TMLE_truncation.out
