source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra)
library(speedglm)
library(foreach); library(doParallel)
library(robustbase)

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
# beta <- 7/8; gamma <- 1/16
true_rate <- 1 / (2 * (gamma + 1 - beta))

# Define Estimation tasks
eta <- 2
ns <- 10^seq(from = 3, to = 7.5, by = 0.5)
candidate_rates <- c(0.8 * true_rate, true_rate, 1.2 * true_rate)
estimation_tasks <- expand.grid(n = ns, candidate_rate = candidate_rates)

# Generate one sample
observed_data.list <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, max(ns))
observed_data <- data.frame(L0 = observed_data.list$L0,
                            A0 = observed_data.list$A0,
                            L1 = observed_data.list$L1)

cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
registerDoParallel(cl)

results <- foreach(job = 1:nrow(estimation_tasks), .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     n <- estimation_tasks[job, ]$n
                     candidate_rate <- estimation_tasks[job, ]$candidate_rate
                     cat("Beginning of n = ", n, "\n")
                     delta <- n^(-candidate_rate / eta)
                     Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
                     TMLE_delta <- TMLE_EY1_speedglm(observed_data[1:n, ], delta)
                     TMLE_delta_plus_Delta <- TMLE_EY1_speedglm(observed_data[1:n, ], delta + Delta)
                     sigma_n <- sqrt(TMLE_delta$var_IC)
                     fin_diff <- (TMLE_delta_plus_Delta$Psi_n - TMLE_delta$Psi_n) / Delta
                     cat("n = ", n, ", delta = ", delta, ", sigma_n = ", sigma_n, "\n")
                     c(n = n, candidate_rate = candidate_rate, delta = delta, sigma_n = sigma_n, fin_diff = fin_diff)
                   }

log_sigma_n_log_delta <- ggplot(as.data.frame(results), aes(x = log(delta) / log(10),
                                                   y = log(sigma_n) / log(10),
                                                   coulour = factor(candidate_rate))) +
  geom_point() + geom_line() +
  geom_abline(intercept = -0.7, slope = -gamma) +
  geom_abline(intercept = -0.8, slope = -gamma)

log_sigma_n_log_n <- ggplot(as.data.frame(results), aes(x = log(n) / log(10),
                                                        y = log(sigma_n) / log(10),
                                                        coulour = factor(candidate_rate))) +
  geom_point() + geom_line()

log_fin_diff_log_delta <- ggplot(as.data.frame(results), aes(x = log(delta) / log(10),
                                                             y = log(abs(fin_diff)) / log(10),
                                                             coulour = candidate_rate)) +
  geom_point() + geom_line()

log_fin_diff_log_n <- ggplot(as.data.frame(results), aes(x = log(n) / log(10),
                                                             y = log(abs(fin_diff)) / log(10),
                                                             coulour = factor(candidate_rate))) +
  geom_point() + geom_line()

grid.arrange(nrow = 2, ncol = 2, log_sigma_n_log_delta, log_sigma_n_log_n,
             log_fin_diff_log_delta, log_fin_diff_log_n)