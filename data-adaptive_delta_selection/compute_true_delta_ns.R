#  ------------------------------------------------------------------------
# 
# Compute the MSE as a function of n and delta in order to eventually
# compute  delta^0_n = argmin E[(\hat{\Psi}_n(\delta) - EYd)^2]
# See write up for more details and notations.
#
#  ------------------------------------------------------------------------

source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')

# Target parameter definition
type <- "L0_exp"; positivity_parameter <- 2; alpha0 <- 2
beta0 <- -3; beta1 <- 1.5; beta2 <-  1
d0 <- alwaysTreated0

# Co

# Simulation parameters
ns <- c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6)
deltas <- c(seq(from = 0, to = 1e-2, length = 100),
            seq(from = 1e-2, to = 2e-1, length = 100))
estimation_tasks <- expand.grid(n = ns, delta = deltas)
jobs <- sample(1:nrow(estimation_tasks)) # A job is an estimation task id

delta_ns <- matrix(nrow = length(ns), ncol = 4)

for(n in ns){
  
}