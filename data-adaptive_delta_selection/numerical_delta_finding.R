source('../utilities.R')
# source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(speedglm)

Rn <- function(observed_data, delta){
  n <- length(observed_data$L0)
  Delta <- n^-0.25 * delta^0.8
  
  TMLE_delta.result <- TMLE_EY1_speedglm(observed_data, delta, F)
  Psi_n_delta <- TMLE_delta.result$Psi_n
  sigma_n <- sqrt(TMLE_delta.result$var_IC)
  
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(observed_data, delta + Delta, F)$Psi_n
  
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta + 1 / sqrt(n) * sigma_n
}


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
n <- 5000

observed_data <- generate_data('L0_exp', lambda, alpha0, beta0, beta1, beta2, n)
deltas <- 10^seq(from = -6, to = -0.5, length = 100)
Rns <- sapply(deltas, Vectorize(function(delta) Rn(observed_data, delta)))

plot(deltas, Rns)
