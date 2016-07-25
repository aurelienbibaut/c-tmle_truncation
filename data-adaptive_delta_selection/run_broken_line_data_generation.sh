#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M aurelien.bibaut@berkeley.edu 

mpirun -n 1 R --vanilla < broken_line_find_gamma-generate_datasets-mpi.R  > broken_line_find_gamma-generate_datasets-mpi.out

