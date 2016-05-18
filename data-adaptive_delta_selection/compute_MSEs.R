#  ------------------------------------------------------------------------
# 
# Compute the MSE as a function of n and delta in order to eventually
# compute  delta^0_n = argmin E[(\hat{\Psi}_n(\delta) - EYd)^2]
# See write up for more details and notations.
#
#  ------------------------------------------------------------------------

source('../utilities.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')
source('../true_target_parameters_derivatives_and_ICs.R')

# Target parameter definition
type <- "L0_exp"; positivity_parameter <- 2; alpha0 <- 2
beta0 <- -3; beta1 <- 1.5; beta2 <-  1
d0 <- alwaysTreated0

# Simulation parameters
M <- 1e4 # Number of data generations per (n, delta)
ns <- c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6)
deltas <- c(seq(from = 0, to = 1e-2, length = 100),
            seq(from = 1e-2, to = 2e-1, length = 100))
estimation_tasks <- expand.grid(n = ns, delta = deltas)
jobs <- sample(1:nrow(estimation_tasks)) # A job is an estimation task id

# Compute the true EYd
EYd <- compute_true_Psi0_delta(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, delta = 0)

# Perform the simulations
library(Rmpi); library(doMPI)

cl <- startMPIcluster(72)
registerDoMPI(cl)

foreach(i = 1:length(jobs)) %dopar% {
# i <- 1
  job <- jobs[i]
  estimates <- vector()
  for(m in 1:M){
    observed_data <- generate_data(type, positivity_parameter, alpha0, beta0, beta1, beta2, estimation_tasks[job, ]$n)
    estimates <- c(estimates, TMLE_truncated_target(observed_data, d0, estimation_tasks[job, ]$delta, Q_misspecified = F))
  }
  MSE <- mean((estimates - EYd)^2)
  iteration_result <- t(c(estimation_tasks[job, ]$n, estimation_tasks[job, ]$delta, MSE))
  colnames(iteration_result) <- c("n", "delta", "MSE")
  
  if(!file.exists("MSE_n,delta.csv")){
    write.table(iteration_result, file = "MSE_n,delta.csv", append = T, row.names = F, col.names = T,  sep = ",")
  }else{
    write.table(iteration_result, file = "MSE_n,delta.csv", append = T, row.names = F, col.names = F,  sep = ",")
  }
}