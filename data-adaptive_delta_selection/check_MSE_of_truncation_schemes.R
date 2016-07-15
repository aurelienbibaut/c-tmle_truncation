source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

# library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)
# library(boot); library(segmented)

# Specify data-generating distribution ------------------------------------

# Set of parameters 1
# lambda <- 2; alpha0 <- 2; beta0 <- -3; beta1 <- 1.5; beta2 <- 1
# beta <- 0.25
# gamma <- 0.125
# kappa <- 1 / (2 * (gamma + 1 - beta))
# Set of parameters 2
# lambda <- 2; alpha0 <- 4; beta0 <- -3; beta1 <- 1.5; beta2 <- 0.5
# kappa <- 5 / 4
# beta <- 2 - kappa
# gamma <- 1 - kappa / 2
# Set of parameters 3
lambda <- 2; alpha0 <- 4; beta2 <- -3; beta0 <- -1; beta1 <- 1
beta <- 7/8; gamma <- 5/16

# Compute true target parameter
EY1 <- compute_true_Psi0_delta('L0_exp', lambda, alpha0, beta0, beta1, beta2, alwaysTreated0, 0)

# Simulate data
n <- 10^4.5

# Simulations parameters
nb_samples <- 10^5

# Set up cluster
cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()), outfile = '')
registerDoParallel(cl)

# Set up tasks
results <- foreach(i=1:nb_samples, .combine = rbind,
                   .packages = c('speedglm', 'boot'), .verbose = T, .inorder = F) %dopar% {
                     
                     observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, n)
                     Psi_n_0 <- TMLE_EY1_speedglm(observed_data, 0, F)
                     Psi_n_0.05 <- TMLE_EY1_speedglm(observed_data, 0.05, F)
                     delta_n <- 0.05 * (n / n0)^-optimal_rate
}