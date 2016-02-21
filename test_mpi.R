dumbFunction <- function(mean){
  x <- rnorm(1e7, mean = mean, sd = 1)
  mean(x)
}


library(Rmpi); library(doMPI)

cl <- startMPIcluster(32)
registerDoMPI(cl)

results <- foreach(i=1:34){
  dumbFunction(i) 
}

print(results)

closeCluster(cl)
mpi.quit()