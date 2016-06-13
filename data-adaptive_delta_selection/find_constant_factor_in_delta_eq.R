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

finite_diff_estimator <- function(observed_data, delta, n){
  Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
  
  Psi_n_delta_plus_Delta <- TMLE_EY1(observed_data, delta + Delta)
  result_TMLE_delta <- TMLE_EY1_bootstrap(observed_data, delta)
  Psi_n_delta <- result_TMLE_delta$Psi_n
  
  return(list(fin_diff = (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta,
              p_val = result_TMLE_delta$shapiro.p_value))
}


# ns <- c(10^(2:6) , 2e6)
ns <- c(10^(2:5), 2e5)
# ns <- floor(c(10^seq(from = 2, to = log(2e6) / log(10), length = 20)))
# etas <- c(1, 1.1, 1.5, 2, 3, 4)
etas <- c(1, 2, 4)
optimal_rate <- -1 / (2 * (gamma + 1 - beta))
rates <- c(0.95 * optimal_rate, 0.99 * optimal_rate,
           optimal_rate, 
           1.01 * optimal_rate, 1.05 * optimal_rate)
# rates <- sort(c(seq(from = -1 / (2 * 0.8 * (gamma + 1 - beta)), 
#                     to =  -1 / (2 * 4 * (gamma + 1 - beta)), 
#                     length = 5), 
#                 0.95 * optimal_rate, 0.99 * optimal_rate,
#                 optimal_rate, 
#                 1.01 * optimal_rate, 1.05 * optimal_rate))

jobs <- expand.grid(n = ns, eta = etas, rate = rates)

observed_data <- generate_data("L0_exp", 2, 2, -3, 1.5, 1, max(ns))

# Set up cluster
cl <- makeCluster(getOption("cl.cores", 40), outfile = '')
registerDoParallel(cl)


results <- foreach(i=1:nrow(jobs), .combine = rbind, 
                   .packages = c("boot", "stats"), .verbose = T) %dopar% {
                     n <- jobs[i, ]$n; rate <- jobs[i, ]$rate;  eta <- jobs[i, ]$eta
                     
                     delta <- n^(rate * 1 / eta)
                     
                     cat('About to compute the LHS for n=', n, ', rate = ', rate, 'eta =', eta, '\n')
                     fin_diff_result <- finite_diff_estimator(lapply(observed_data, function(x) x[1:n]), delta, n)
                     LHS <- abs(fin_diff_result$fin_diff * delta^2)^eta * n
                     
                     c(n, rate, eta, LHS, fin_diff_result$p_val)
                   }

colnames(results) <- c("n", "base_rate", "eta", "LHS", "p_value")
# Plot the results
results_df <- as.data.frame(results)
results_df <- transform(results_df, base_rate = factor(round(base_rate, 3)))
LHS_plots <- list()
for(i in 1:length(etas)){
  LHS_plots[[i]] <- ggplot(data = subset(results_df, eta == etas[i]), aes(x = log(n)/log(10), y = log(LHS), colour = base_rate)) +
    geom_line() + geom_point(aes(size = p_value)) + 
    ggtitle(substitute(list(eta) == list(x),
                                     list(x = etas[i])))
}

# grid.arrange(LHS_plots[[1]], LHS_plots[[2]],
#              LHS_plots[[3]], LHS_plots[[4]],
#              LHS_plots[[5]], LHS_plots[[6]], nrow = 2, ncol = 3)

grid.arrange(LHS_plots[[1]], LHS_plots[[2]], LHS_plots[[3]], nrow = 2, ncol = 2)
