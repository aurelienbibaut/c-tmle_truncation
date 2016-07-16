source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
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
lambda <- 2; alpha0 <- 4; beta0 <- -3; beta1 <- 1.5; beta2 <- 0.5
kappa <- 5 / 4
beta <- 2 - kappa
gamma <- 1 - kappa / 2
# Set of parameters 3
# lambda <- 2; alpha0 <- 4; beta2 <- -3; beta0 <- -1; beta1 <- 1
# beta <- 7/8; gamma <- 5/16

# Define optimal rate based truncation parameters
optimal_rate <- 1 / (2 * (gamma + 1 - beta))
n0 <- 1000

# Compute true target parameter
EY1 <- compute_true_Psi0_delta('L0_exp', lambda, alpha0, beta0, beta1, beta2, alwaysTreated0, 0)

# Simulate data
ns <- c(1e3, 5e3, 10^4, 2e4, 3e4, 4e4)

# Simulations parameters
nb_samples <- 10^3

# Set up cluster
cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()), outfile = '')
registerDoParallel(cl)

# Set up tasks
jobs <- kronecker(ns, rep(1, nb_samples))

results <- foreach(n=jobs, .combine = rbind,
                   .packages = c('speedglm', 'boot'), .verbose = T, .inorder = F) %dopar% {
                     
                     observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, n)
                     Psi_n_0 <- TMLE_EY1_speedglm(observed_data, 0, F)$Psi_n
                     Psi_n_0.05 <- TMLE_EY1_speedglm(observed_data, 0.05, F)$Psi_n
                     delta_n_1 <- 0.05 * (n / n0)^(-optimal_rate)
                     delta_n_2 <- 0.02 * (n / n0)^(-optimal_rate)
                     delta_n_3 <- 0.05 * (n / n0)^(-1)
                     delta_n_4 <- 0.05 * (n / n0)^(-1)
                     
                     Psi_n_delta_n_1 <- TMLE_EY1_speedglm(observed_data, delta_n_1, F)$Psi_n
                     Psi_n_delta_n_2 <- TMLE_EY1_speedglm(observed_data, delta_n_2, F)$Psi_n
                     Psi_n_delta_n_3 <- TMLE_EY1_speedglm(observed_data, delta_n_3, F)$Psi_n
                     Psi_n_delta_n_4 <- TMLE_EY1_speedglm(observed_data, delta_n_4, F)$Psi_n
                     
                     rbind(c(n, '0', Psi_n_0),
                           c(n, '0.05', Psi_n_0.05),
                           c(n, 'delta_n_1', Psi_n_delta_n_1),
                           c(n, 'delta_n_2', Psi_n_delta_n_2),
                           c(n, 'delta_n_3', Psi_n_delta_n_3),
                           c(n, 'delta_n_4', Psi_n_delta_n_4))
                     
                   }
colnames(results) <- c("n", "truncation_scheme", "estimate")
results_df <- transform(results, as.data.frame(results), n = as.numeric(as.character(n)),
                        truncation_scheme = as.character(truncation_scheme),
                        estimate = as.numeric(as.character(estimate)))
estimation_tasks <- unique(results_df[, 1:2])
row.names(estimation_tasks) <- NULL
MSEs <- vector()
for(i in 1:nrow(estimation_tasks)){
  MSEs <- rbind(MSEs,
                c(estimation_tasks[i, 1], estimation_tasks[i, 2],
                  mean((subset(results_df, n == estimation_tasks[i, 1] & 
                                 truncation_scheme == estimation_tasks[i, 2],
                               select = estimate)[[1]] - EY1)^2)))
}
MSEs_df <- data.frame(n = as.numeric(MSEs[, 1]), truncation_scheme = MSEs[, 2],
                      MSE = as.numeric(MSEs[, 3]))
MSE.plot <- ggplot(subset(MSEs_df, truncation_scheme != 0), 
                   aes(x = n, y = MSE, colour = truncation_scheme)) + geom_point() + geom_line()
print(MSE.plot)
