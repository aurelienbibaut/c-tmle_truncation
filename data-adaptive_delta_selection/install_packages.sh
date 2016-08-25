#!/bin/bash
module load gcc/4.4.7
module load r openmpi Rmpi
module load foreach
module load iterators


Rscript -e 'install.packages("speedglm_0.3-1.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

Rscript -e 'install.packages("segmented_0.5-1.4.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

Rscript -e 'install.packages("boot_1.3-18.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

Rscript -e 'install.packages("doMPI_0.2.1.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

Rscript -e 'install.packages("R.methodsS3_1.7.1.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

Rscript -e 'install.packages("R.oo_1.20.0.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

Rscript -e 'install.packages("R.utils_2.3.0.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

Rscript -e 'install.packages("DEoptimR_1.0-6.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

Rscript -e 'install.packages("robustbase_0.92-6.tar.gz", repos = NULL, type = "source", lib = "/global/home/users/afbibaut/Rlibs", INSTALL_opts = c("--no-lock"))' --verbose

