#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M aurelien.bibaut@berkeley.edu 

mpirun -n 1 R --vanilla < test_mpi.R > test_mpi.out

