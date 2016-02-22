#!/usr/bin/env Rscript

dumbFunction <- function(mean){
  x <- rnorm(1e7, mean = mean, sd = 1)
  mean(x)
}


library(Rmpi); library(doMPI)

cl <- startMPIcluster(32*3)
registerDoMPI(cl)

results <- foreach(i=1:3000) %dopar% {
  dumbFunction(i) 
}

print(results)
save(result, file = "test_results")

closeCluster(cl)
mpi.quit()