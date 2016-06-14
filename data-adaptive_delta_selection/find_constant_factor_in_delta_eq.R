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
  
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}

finite_diff_estimator.subsampling <- function(observed_data, base_rate, eta, n, nb_subsamples = 5, C = 1 / sqrt(10)){
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
  full_sample.LHS <- (full_sample_fin_diff * delta^2)^eta * n
  subsamples.LHS <- (subsampled_fin_diffs * delta^2)^eta * n
  median(full_sample.LHS - subsamples.LHS)
}

# debug(finite_diff_estimator.subsampling)

# ns <- c(10^(2:6) , 2e6)
# ns <- 10^(2:6)
# ns <- 10^(5:6)
# ns <- c(10^(2:5), 2e5)
ns <- floor(c(10^seq(from = 2, to = log(2e6) / log(10), length = 20)))
# etas <- c(1, 1.1, 1.5, 2, 3, 4)
# etas <- c(1, 2, 3, 4)
etas <- c(1, 1.1, 1.5, 2, 3, 4)
optimal_rate <- -1 / (2 * (gamma + 1 - beta))
# rates <- c(0.8 * optimal_rate, 0.95 * optimal_rate, 0.99 * optimal_rate,
#            optimal_rate, 
#            1.01 * optimal_rate, 1.05 * optimal_rate, 1.2 * optimal_rate)
rates <- seq(from = -0.7, to = -1/3, length = 9)
# rates <- sort(c(seq(from = -1 / (2 * 0.8 * (gamma + 1 - beta)), 
#                     to =  -1 / (2 * 4 * (gamma + 1 - beta)), 
#                     length = 5), 
#                 0.95 * optimal_rate, 0.99 * optimal_rate,
#                 optimal_rate, 
#                 1.01 * optimal_rate, 1.05 * optimal_rate))
C <- 1 / sqrt(10)
jobs <- expand.grid(n = ns, eta = etas, rate = rates, replicate = 1)

observed_data <- generate_data("L0_exp", 2, 2, -3, 1.5, 1, max(ns))


# Set up cluster
cl <- makeCluster(getOption("cl.cores", 40), outfile = '')
registerDoParallel(cl)


results <- foreach(i=1:nrow(jobs), .combine = rbind, 
                   .packages = c("boot", "stats"), .verbose = T) %dopar% {
                     n <- jobs[i, ]$n; rate <- jobs[i, ]$rate;  eta <- jobs[i, ]$eta
                     replicate <- jobs[i, ]$replicate
                     
                     # delta <- 1 / sqrt(10) * n^(rate * 1 / eta)
                     delta <- C * n^(rate * 1 / eta)
                     
                     cat('About to compute the LHS for n=', n, ', rate = ', rate, 'eta =', eta, '\n')
                     # indices <- sample(1:length(observed_data$L0), n, F)
                     indices <- 1:n
                     # fin_diff_result <- finite_diff_estimator.bootstrap(lapply(observed_data, function(x) x[indices]), delta, n)
                     fin_diff <- finite_diff_estimator(lapply(observed_data, function(x) x[indices]), delta, n)
                     LHS <- abs(fin_diff * delta^2)^eta * n
                     
                     c(n, rate, eta, replicate, LHS)
                   }

colnames(results) <- c("n", "base_rate", "eta", "replicate", "LHS")
# Plot the results
results_df <- as.data.frame(results)

results_df_bis <- data.frame(n = results_df$n, 
                             track = apply(cbind(results_df$base_rate, results_df$replicate), 
                                           1, function(x) paste(c(x[1], x[2]), collapse = "_")),
                             eta = results_df$eta,
                             base_rate = results_df$base_rate,
                             LHS = results_df$LHS)

results_df_bis <- transform(results_df_bis, base_rate = factor(round(base_rate, 3)))

LHS_plots <- list()
for(i in 1:length(etas)){
  LHS_plots[[i]] <- ggplot(data = subset(results_df_bis, eta == etas[i] &
                                           n >= 1e2), aes(x = log(n)/log(10), y = log(LHS),
                                                          group = track,
                                                          colour = base_rate)) +
    geom_line() + 
    # geom_point(aes(size = p_value)) + 
    ggtitle(substitute(list(eta) == list(x),
                       list(x = etas[i])))
}
# 
grid.arrange(LHS_plots[[1]], LHS_plots[[2]],
             LHS_plots[[3]], LHS_plots[[4]],
             LHS_plots[[5]], LHS_plots[[6]], nrow = 2, ncol = 3, 
             top = paste(c("C = ", C), collapse = ''))
# 
# grid.arrange(LHS_plots[[1]], LHS_plots[[2]], LHS_plots[[3]], LHS_plots[[4]],
#              nrow = 2, ncol = 2)
# grid.arrange(LHS_plots[[1]], LHS_plots[[2]], ncol = 2, nrow = 1,
#              top = paste(c("C = ", C), collapse = ''))

