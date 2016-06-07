library(doParallel); library(foreach)
library(ggplot2); library(gridExtra)

source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')


beta <- 0.25
gamma <- 0.125
finite_diffs <- vector(); finite_diffs_bias <- vector()

deltas <- vector()

finite_diff_estimator <- function(observed_data, delta, n){
  Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
  
  Psi_n_delta_plus_Delta <- TMLE_EY1(observed_data, delta + Delta)
  Psi_n_delta <- TMLE_EY1(observed_data, delta)
  
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}


ns <- c(10^(2:6) , 2e6)
etas <- c(1, 1.1, 1.5, 2, 3, 4)
rates <- seq(from = -1 / (2 * 0.8 * eta_n * (gamma + 1 - beta)), 
             to =  -1 / (2 * 4 * eta_n * (gamma + 1 - beta)), 
             length = 5)
jobs <- expand.grid(n = ns, eta = etas, rate = rates)


# Set up cluster
cl <- makeCluster(getOption("cl.cores", 12), outfile = '')
registerDoParallel(cl)


results <- foreach(i=1:nrow(jobs), .combine = rbind) %dopar% {
  n <- jobs[i, ]$n; rate <- jobs[i, ]$rate;  eta <- jobs[i, ]$eta
  observed_data <- generate_data("L0_exp", 2, 2, -3, 1.5, 1, n)
  delta <- n^rate
  
  cat('About to compute the LHS for n=', n, ', rate = ', rate, 'eta =', eta, '\n')
  LHS <- abs(finite_diff_estimator(observed_data, delta, n) * delta^2)^eta * n
  
  c(n, rate, eta, LHS)
}

# Plot the results
results_df <- as.data.frame(results)
LHS_plots <- list()
for(i in 1:length(etas)){
  LHS_plots[[i]] <- ggplot(data = subset(results_df, eta == etas[i]), aes(x = log(n), y = log(LHS), colour = rate)) +
    geom_line()
}