source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(speedglm)
library(foreach); library(doParallel)

# debug(TMLE_EY1)
# deltas <- seq(from = 0, to = 10^-2, length = 100)
# 
# Psi0s <- sapply(deltas, function(delta) compute_true_Psi0_delta("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta))
# 
# Psi0_0 <- Psi0s[1]
# b0 <- abs(Psi0s - Psi0_0)
# 
# rate_b0 <- (log(b0[100]) - log(b0[2])) / (log(deltas[100]) - log(deltas[2]))
# 
# sigma0s <- sapply(deltas, function(delta) sqrt(compute_true_var_IC_Psi0_delta("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta)))
# 
# rate_sigma0 <- (log(sigma0s[100]) - log(sigma0s[2])) / (log(deltas[100]) - log(deltas[2]))
# 
# par(mfrow = c(1, 2))
# plot(log(deltas), log(b0))
# abline(-4, 0.75)
# 
# plot(log(deltas), log(sigma0s))
# print(rate_b0)
# print(rate_sigma0)
# 
# cat("beta = ", 1 - rate_b0, "\ngamma = ", abs(rate_sigma0))

# Now let's check if (Psi_n(delta + Delta) - Psi(delta) + n^(-1/2) * ((delta+Delta)^-gamma - delta^-gamma)) / Delta
beta <- 0.25
gamma <- 0.125
etas <- c(1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
finite_diffs <- vector(); finite_diffs_bias <- vector(); true_finite_diffs_bias <- vector()
ns <- floor(10^seq(from = 2, to = 8, by = 0.5))
deltas <- vector()

observed_data <- generate_data("L0_exp", 2, 2, -3, 1.5, 1, max(ns))
# data <- data.frame(L0 = observed_data$L0, A0 = observed_data$A0, L1 = observed_data$L1)



jobs <- expand.grid(n = ns, eta = etas)

# Set up cluster
cl <- makeCluster(getOption("cl.cores", 4), outfile = '')
registerDoParallel(cl)

results <- foreach(i=1:nrow(jobs), .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     n <- jobs[i, ]$n; eta <- jobs[i, ]$eta
                     
                     delta_n_plus <- n^(-1 / (2 * eta * (gamma + 1 - beta)))
                     
                     deltas <- c(deltas, delta_n_plus)
                     
                     Delta <- n^(-0.25) * delta_n_plus^((beta + 1 - gamma) / 2)
                     
                     cat('Job ', i, ', n = ', n, ', eta = ', eta, ', delta =', delta_n_plus, ', Delta = ', Delta, '\n')
                     Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_plus + Delta)
                     cat('Psi_n(delta + Delta) = ', Psi_n_delta_plus_Delta, '\n')
                     Psi_n_delta <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_plus)
                     cat('Psi_n(delta) = ', Psi_n_delta, '\n')
                     
                     finite_diffs <- c(finite_diffs,
                                       (Psi_n_delta_plus_Delta - Psi_n_delta + 
                                          n^(-0.5) * ((delta_n_plus + Delta)^(-gamma) - delta_n_plus^(-gamma))) / Delta)
                     
                     finite_diffs_bias <- c(finite_diffs_bias, (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta)
                     Psi_0_delta_plus_Delta <- compute_true_Psi0_delta("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta_n_plus + Delta)
                     Psi_0_delta <- compute_true_Psi0_delta("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta_n_plus)
                     
                     true_finite_diffs_bias <- c(true_finite_diffs_bias, (Psi_0_delta_plus_Delta - Psi_0_delta) / Delta)
                     
                     print(finite_diffs)
                     print(finite_diffs_bias)
                     
                     #   par(mfrow = c(1, 2))
                     #   # plot(log(deltas), log(abs(finite_diffs)))
                     #   plot(log(deltas), log(abs(finite_diffs_bias)))
                     #   plot(log(deltas), log(abs(true_finite_diffs_bias)))
                     c(eta, n, finite_diffs, finite_diffs_bias, true_finite_diffs_bias)
                   }

stopCluster(cl)

plot(log(deltas), log(abs(finite_diffs_bias)), ylim = c(min(c(log(abs(finite_diffs_bias)), log(abs(true_finite_diffs_bias)))),
                                                        max(c(log(abs(finite_diffs_bias)), log(abs(true_finite_diffs_bias))))),
     main= substitute(list(eta) == list(x),
                      list(x = eta)))
lines(log(deltas), log(abs(true_finite_diffs_bias)))
