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
delta_n_star <- vector()
etas_n <- vector()

g_a <- vector(); g_b <- vector()

for(n in ns){
  cat("n = ", n, "\n")
  observed_data <- generate_data("L0_exp", 2, 2, -3, 1.5, 1, n)
  
  
  
  # eta_n <- 1 - 1e-2
  eta_n <- 2
  #   class(current_solver_result) <- "try-error"
  #   while(class(current_solver_result) == "try-error"){
  #     
  #     eta_n <- eta_n + 1e-2
  #     
  #     delta_n_a <- n^(-1 / (2 * 0.5 * eta_n * (gamma + 1 - beta)))
  #     delta_n_b <- n^(-1 / (2 * 2 * eta_n * (gamma + 1 - beta)))
  #     
  # #     eq_func <- Vectorize(function(delta) (abs(finite_diff_estimator(observed_data, delta, n)) * delta^2)^eta_n - 1 / n)
  # #     LHS <- Vectorize(function(delta) abs(finite_diff_estimator(observed_data, delta, n)) * delta^2)
  # #     deltas <- seq(from = delta_n_a, to = delta_n_b, length = 100)
  # #     plot(deltas, eq_func(deltas))
  # #     abline(0, 0)
  #     current_solver_result <- try(uniroot(Vectorize(function(delta) abs(finite_diff_estimator(observed_data, delta, n) * delta^2)^eta_n - 1 / n),
  #                                                    interval = c(delta_n_a, delta_n_b),
  #                                                    extendInt = "yes")$root)
  #   }
  
  delta_n_a <- n^(-1 / (2 * 0.5 * eta_n * (gamma + 1 - beta)))
  delta_n_star <- n^(-1 / (2 *  eta_n * (gamma + 1 - beta)))
  delta_n_b <- n^(-1 / (2 * 4 * eta_n * (gamma + 1 - beta)))
  
  g_a <- c(g_a, abs(finite_diff_estimator(observed_data, delta_n_a, n) * delta_n_a^2)^eta_n * n)
  g_b <- c(g_b, abs(finite_diff_estimator(observed_data, delta_n_b, n) * delta_n_b^2)^eta_n * n)
  g_star <- c(g_star, abs(finite_diff_estimator(observed_data, delta_n_star, n) * delta_n_star^2)^eta_n * n)
  
  plot(log(ns[which(ns <= n)]), log(g_a), type = 'l', col = 'green', ylim = c(min(c(log(g_a[g_a != 0]), log(g_b[g_b != 0]))),
                                                                              max(c(log(g_a[g_a != 0]), log(g_b[g_b != 0])))))
  lines(log(ns[which(ns <= n)]), log(g_b), col = 'blue')
  lines(log(ns[which(ns <= n)]), log(g_star), col = 'red')
  
  #   delta_n_star <- c(delta_n_star, current_solver_result)
  #   etas_n <- c(etas_n, eta_n)
}