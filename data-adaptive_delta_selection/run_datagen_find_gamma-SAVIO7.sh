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
#SBATCH --ntasks=48
#
# Nodes:
#SBATCH -N 2
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Mail type:
# SBATCH --mail-type=all

#SBATCH --mail-user=aurelien.bibaut@berkeley.edu
#

module swap intel/2013_sp1.4.211 gcc/4.4.7
module load r openmpi Rmpi
Rscript -e 'install.packages("speedglm_0.3-1.tar.gz", repos = NULL, type = "source", lib = "~/Rlibs", INSTALL_opts = c("--no-lock"))'

Rscript -e 'install.packages("segmented_0.5-1.4.tar.gz", repos = NULL, type = "source", lib = "~/Rlibs", INSTALL_opts = c("--no-lock"))'

Rscript -e 'install.packages("boot_1.3-18.tar.gz", repos = NULL, type = "source", lib = "~/Rlibs", INSTALL_opts = c("--no-lock"))'

Rscript -e 'install.packages("doMPI_0.2.1.tar.gz", repos = NULL, type = "source", lib = "~/Rlibs", INSTALL_opts = c("--no-lock"))'

Rscript -e 'install.packages("R.oo_1.20.0.tar.gz", repos = NULL, type = "source", lib = "~/Rlibs", INSTALL_opts = c("--no-lock"))'

Rscript -e 'install.packages("R.methodsS3_1.7.1.tar.gz", repos = NULL, type = "source", lib = "~/Rlibs", INSTALL_opts = c("--no-lock"))'

Rscript -e 'install.packages("R.utils_2.3.0.tar.gz", repos = NULL, type = "source", lib = "~/Rlibs", INSTALL_opts = c("--no-lock"))'

mpirun R CMD BATCH --no-save broken_line_find_gamma-generate_datasets-mpi.R broken_line_find_gamma-generate_datasets-mpi7.out
