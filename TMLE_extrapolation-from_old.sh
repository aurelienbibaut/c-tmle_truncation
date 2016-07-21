#!/bin/bash
# Job name:
#SBATCH --job-name=Taylor_expansion_TMLE_raw_estimates
#
# Partition:
#SBATCH --partition=savio2
#
# Account:
#SBATCH --account=co_biostat
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
./Taylor_extrapolation-inference_mean_Yd-from_old.R