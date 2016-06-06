source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')


beta <- 0.25
gamma <- 0.125
finite_diffs <- vector(); finite_diffs_bias <- vector()
ns <- 10^(2:7)
deltas <- vector()

finite_diff_estimator <- function(observed_data, delta, n){
  Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
  
  Psi_n_delta_plus_Delta <- TMLE_EY1(observed_data, delta + Delta)
  Psi_n_delta <- TMLE_EY1(observed_data, delta)
  
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}

current_solver_result <- NA
delta_n_star <- NA
etas_n <- NA

for(n in ns){
  cat("n = ", n, "\n")
  observed_data <- generate_data("L0_exp", 2, 2, -3, 1.5, 1, n)
  
  eta_n <- 1 - 1e-2
  class(current_solver_result) <- "try-error"
  while(class(current_solver_result) == "try-error"){
    eta_n <- eta_n + 1e-2
    current_solver_result <- try(uniroot(Vectorize(function(delta) (finite_diff_estimator(observed_data, delta, n) * delta^2)^eta - 1 / n,
                                         interval = c(delta_n_a, delta_n_b),
                                         extendInt = "yes")$root))
  }
  delta_n_star <- c(delta_n_star, current_solver_result)
  etas_n <- c(etas_n, eta_n)
}