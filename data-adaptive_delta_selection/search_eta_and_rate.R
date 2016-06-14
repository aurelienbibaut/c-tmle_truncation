library(doParallel); library(foreach)
library(ggplot2); library(gridExtra)
library(boot)
library(stats)

source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')


beta <- 0.25
gamma <- 0.125
finite_diffs <- vector(); finite_diffs_bias <- vector()

deltas <- vector()

finite_diff_estimator.bootstrap <- function(observed_data, delta, n){
  Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
  
  Psi_n_delta_plus_Delta <- TMLE_EY1(observed_data, delta + Delta)
  result_TMLE_delta <- TMLE_EY1_bootstrap(observed_data, delta, nb_boostrap_samples = 4)
  Psi_n_delta <- result_TMLE_delta$Psi_n
  
  return(list(fin_diff = (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta,
              p_val = result_TMLE_delta$shapiro.p_value))
  # TMLE_fin_diff.boot <- TMLE_fin_diff_EY1_bootstrap(observed_data, delta, 
  #                                                   Delta, nb_boostrap_samples = 4)
  # return(list(fin_diff = TMLE_fin_diff.boot$Psi_n / Delta, p_val = TMLE_fin_diff.boot$shapiro.p_value))
}

finite_diff_estimator <- function(observed_data, delta, n){
  Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
  Psi_n_delta_plus_Delta <- TMLE_EY1(observed_data, delta + Delta)
  Psi_n_delta <- TMLE_EY1(observed_data, delta)
  
  (Psi_n_delta_plus_Delta - Psi_n_delta)/Delta
}

finite_diff_estimator.subsampling <- function(observed_data, base_rate, eta, n, nb_subsamples = 5, C = 1){
  if(rate > 0) cat('Warning: positive rate passed\n')
  # Full sample estimation
  delta <- C * n^(base_rate / eta)
  full_sample_fin_diff <- finite_diff_estimator(observed_data, delta, n)
  
  # Subsampling
  n_sub <- floor(n / sqrt(10))
  delta_sub <- C * n_sub^(base_rate / eta)
  indices <- t(replicate(nb_subsamples, sample(1:n, n_sub, F)))
  subsampled_fin_diffs <- apply(indices, 1, function(y) finite_diff_estimator(lapply(observed_data, function(x) x[y]),
                                                                              delta_sub, n_sub))
  # Compute median difference in LHS between full sample and subsamples
  full_sample.LHS <- (-full_sample_fin_diff * delta^2)^eta * n
  subsamples.LHS <- (-subsampled_fin_diffs * delta^2)^eta * n
  
  if(full_sample_fin_diff == 0 | any(subsamples.LHS) == 0){
    return(NA)
  }else{
    return(median(full_sample.LHS - subsamples.LHS))
  }
}

# debug(finite_diff_estimator.subsampling)


n <- 2e5
# etas <- c(1, 1.1, 1.5, 2, 3, 4)
etas <- c(1, 2, 3, 4)
optimal_rate <- -1 / (2 * (gamma + 1 - beta))
# rates <- c(0.95 * optimal_rate, 0.99 * optimal_rate,
#            optimal_rate, 
#            1.01 * optimal_rate, 1.05 * optimal_rate)
rates <- seq(from = -0.7, to = -0.4, length = 20)
etas <- seq(from = 3.2, to = 5, length = 10)
# rates <- sort(c(seq(from = -1 / (2 * 0.8 * (gamma + 1 - beta)), 
#                     to =  -1 / (2 * 4 * (gamma + 1 - beta)), 
#                     length = 5), 
#                 0.95 * optimal_rate, 0.99 * optimal_rate,
#                 optimal_rate, 
#                 1.01 * optimal_rate, 1.05 * optimal_rate))

jobs <- expand.grid(eta = etas, rate = rates)

observed_data <- generate_data("L0_exp", 2, 2, -3, 1.5, 1, n)

# Set up cluster
cl <- makeCluster(getOption("cl.cores", 40), outfile = '')
registerDoParallel(cl)


results <- foreach(i=1:nrow(jobs), .combine = rbind, 
                   .packages = c("boot", "stats"), .verbose = T) %dopar% {
                     rate <- jobs[i, ]$rate;  eta <- jobs[i, ]$eta
                     
                     cat('About to compute the finite difference in LHS for n=', n, ', rate = ', rate, 'eta =', eta, '\n')
                     LHS.fin_diff <- finite_diff_estimator.subsampling(observed_data, base_rate = rate, eta = eta, n)
                     
                     c(n, rate, eta, LHS.fin_diff, 0)
                   }

colnames(results) <- c("n", "base_rate", "eta", "LHS", "p_value")
# Plot the results
results_df <- as.data.frame(results)
results_df <- transform(results_df, base_rate = factor(round(base_rate, 3)), eta = factor(eta))

# LHS.fin_diff.plot <- ggplot(subset(results_df, as.numeric(as.character(eta)) > 3.8), aes(base_rate, eta)) + geom_tile(aes(fill = tanh(1e6 * LHS)))
# print(LHS.fin_diff.plot)

filled.contour(etas, rates, matrix(results_df$LHS, ncol = length(rates)))
# grid.arrange(LHS_plots[[1]], LHS_plots[[2]],
#              LHS_plots[[3]], LHS_plots[[4]],
#              LHS_plots[[5]], LHS_plots[[6]], nrow = 2, ncol = 3)

# grid.arrange(LHS_plots[[1]], LHS_plots[[2]], LHS_plots[[3]], LHS_plots[[4]],
#              nrow = 2, ncol = 2)