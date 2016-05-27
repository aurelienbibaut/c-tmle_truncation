source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

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
finite_diffs <- vector(); finite_diffs_bias <- vector()
ns <- 10^(2:7)
deltas <- vector()

for(n in ns){
  cat("n = ", n, "\n")
  observed_data <- generate_data("L0_exp", 2, 2, -3, 1.5, 1, n)
  delta_n_plus <- n^(-1 / (2 * (gamma + 1 - beta))) * n^0.1
  
  deltas <- c(deltas, delta_n_plus)
  
  Delta <- n^(-0.25) * delta_n_plus^((beta + 1 - gamma) / 2)
  cat("Delta = ", Delta, "\n")
  
  Psi_n_delta_plus_Delta <- TMLE_EY1(observed_data, delta_n_plus + Delta)
  Psi_n_delta <- TMLE_EY1(observed_data, delta_n_plus)
  
  finite_diffs <- c(finite_diffs,
                    (Psi_n_delta_plus_Delta - Psi_n_delta + 
                       n^(-0.5) * ((delta_n_plus + Delta)^(-gamma) - delta_n_plus^(-gamma))) / Delta)
  
  finite_diffs_bias <- c(finite_diffs_bias, (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta)
  
  
  print(finite_diffs)
  print(finite_diffs_bias)
  
  par(mfrow = c(1, 2))
  plot(log(deltas), log(abs(finite_diffs)))
  plot(log(deltas), log(abs(finite_diffs_bias)))
}

