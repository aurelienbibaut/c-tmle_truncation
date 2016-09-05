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
# SBATCH --ntasks=192
#
# Wall clock limit:
#SBATCH --time=48:00:00
#SBATCH -N 2
#SBATCH --array=1-4
# Mail type:
# SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=aurelien.bibaut@berkeley.edu
#

module load gcc/4.4.7
module load r openmpi Rmpi
module load foreach
module load iterators
#export R_LIBS_USER=/global/home/users/afbibaut/Rlibs

# Rscript -e 'install.packages("speedglm_0.3-1.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

# Rscript -e 'install.packages("segmented_0.5-1.4.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

# Rscript -e 'install.packages("boot_1.3-18.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

# Rscript -e 'install.packages("doMPI_0.2.1.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

# Rscript -e 'install.packages("R.methodsS3_1.7.1.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

# Rscript -e 'install.packages("R.oo_1.20.0.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

# Rscript -e 'install.packages("R.utils_2.3.0.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

# Rscript -e 'install.packages("DEoptimR_1.0-6.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

# Rscript -e 'install.packages("robustbase_0.92-6.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

mpirun R CMD BATCH --no-save "--args task_id=${SLURM_ARRAY_TASK_ID}" find_rate-generate_datasets.R find_rate-generate_datasets-${SLURM_ARRAY_TASK_ID}.out
