source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(speedglm); library(ggplot2)
library(foreach); library(doParallel)
library(robustbase)

# Now let's check if (Psi_n(delta + Delta) - Psi(delta) + n^(-1/2) * ((delta+Delta)^-gamma - delta^-gamma)) / Delta
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

# etas <- c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
etas <- c(0.95, 1)
finite_diffs <- vector(); finite_diffs_bias <- vector(); true_finite_diffs_bias <- vector()
# ns <- floor(10^seq(from = 2, to = 4, by = 0.5))
ns <- c(5e3, 1e4)
deltas <- 10^seq(from = -1, to = -7.5, by = -0.05)

observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, max(ns))
# data <- data.frame(L0 = observed_data$L0, A0 = observed_data$A0, L1 = observed_data$L1)


jobs <- expand.grid(n = ns, delta = deltas)

# Set up cluster
cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
registerDoParallel(cl)

results <- foreach(i=1:nrow(jobs), .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     n <- jobs[i, ]$n; delta <- jobs[i, ]$delta
                     
                     cat('Job ', i, ', n = ', n, ', delta =', delta, '\n')
                     
                     indices <- 1:n
                     TMLE_delta.results <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[indices]), delta)
                     var_IC_delta <- TMLE_delta.results$sigma_n
                     gns.quantiles_0.025 <- TMLE_delta.results$gns.quantile_0.025
                     gns.quantiles_0.975 <- TMLE_delta.results$gns.quantile_0.975
                     
                     true_var_IC.result <- compute_true_var_IC_Psi0_delta("L0_exp", lambda, alpha0, beta0, beta1, beta2, 
                                                                          alwaysTreated0, delta)
                     true_var_IC_delta <- true_var_IC.result$var_IC0
                     true_var.first_comp <- true_var_IC.result$first_component
                     true_var.first_comp.a <- true_var_IC.result$first_component.a
                     
                     iteration.result <- c(n, delta, var_IC_delta, 
                                           true_var_IC_delta, true_var.first_comp, true_var.first_comp.a,
                                           gns.quantiles_0.025, gns.quantiles_0.975)
                     cat('\n\n Iteration results\nJob ', i, ', delta = ', delta, ', n = ', n, 
                         '\nvar_IC_delta = ', var_IC_delta, 
                         '\ntrue_var_IC_delta = ', true_var_IC_delta, '\n')
                     iteration.result
                   }

stopCluster(cl)

row.names(results) <- NULL

colnames(results) <- c('n', 'delta', 'var_IC_delta', 
                       'true_var_IC_delta', 'true_var.first_comp', 'true_var.first_comp.a',
                       'gns.quantiles_0.025', 'gns.quantiles_0.975')

results_with_derivatives <- vector()

shrunk_delta_range <- vector()

for(n in ns){
  indices <- results[, 'n'] == n
  results_n <- results[indices, ]
  first_derivative <- c(log(results_n[2:nrow(results_n), 'var_IC_delta']) - 
                          log(results_n[1:(nrow(results_n) - 1), 'var_IC_delta']), NA)
  second_derivative <- c(first_derivative[2:length(first_derivative)] - 
                           first_derivative[1:(length(first_derivative) - 1)], NA)
  log_var_IC_delta <- log(results_n[, 'var_IC_delta']) / log(10)
  log_var_min <- 0.9 * min(log_var_IC_delta) + 0.1 * max(log_var_IC_delta)
  log_var_max <- 0.1 * min(log_var_IC_delta) + 0.9 * max(log_var_IC_delta)
  delta_min_n <- max(results_n[log_var_IC_delta >= log_var_max, 'delta'])
  delta_max_n <- max(results_n[log_var_IC_delta >= log_var_min, 'delta'])
  results_n <- cbind(results_n, first_derivative = first_derivative,
                     second_derivative = second_derivative)
  results_with_derivatives <- rbind(results_with_derivatives,
                                    results_n)
  regression_df <- subset(data.frame(results_n), delta <= delta_max_n & delta >= delta_min_n)
  regression_df <- cbind(regression_df, log_delta = log(regression_df$delta) / log(10),
                         log_var_IC = log(regression_df$var_IC_delta) / log(10))
  line_fit <- ltsReg(regression_df$log_delta, regression_df$log_var_IC)
  shrunk_delta_range <- rbind(shrunk_delta_range,
                              c(n = n, delta_min = delta_min_n, delta_max = delta_max_n,
                                log_var_min = log_var_min,
                                log_var_max = log_var_max))
}


results_df <- as.data.frame(results_with_derivatives)


sigmas_ns.plot <- ggplot(results_df, aes(x = log(delta) / log(10), 
                                         y = log(var_IC_delta) / log(10), 
                                         colour = factor(n))) + 
  geom_line() +
  geom_abline(intercept = -1, slope = - 2 * gamma) +
  geom_hline(mapping = aes(yintercept =  log_var_max, colour = factor(n)), data = as.data.frame(shrunk_delta_range)) +
  geom_hline(mapping = aes(yintercept =  log_var_min, colour = factor(n)), data = as.data.frame(shrunk_delta_range)) +
  geom_vline(mapping = aes(xintercept = log(delta_min) / log(10), colour = factor(n)), data = as.data.frame(shrunk_delta_range)) +
  geom_vline(mapping = aes(xintercept = log(delta_max) / log(10), colour = factor(n)), data = as.data.frame(shrunk_delta_range)) +
  xlab(expression(log[10](delta))) + 
  ylab(expression(log(sigma[n](delta)))) +
  ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                       group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                     list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))
print(sigmas_ns.plot)
# first_derivative_sigmas_ns.plot <- ggplot(results_df, aes(x = log(delta) / log(10), 
#                                                          y = first_derivative,
#                                                          colour = factor(n))) +
#   geom_line()
# 
# second_derivative_sigmas_ns.plot <- ggplot(results_df, aes(x = log(delta) / log(10), 
#                                                           y = second_derivative,
#                                                           colour = factor(n))) +
#   geom_line()
# 
# grid.arrange(sigmas_ns.plot, 
#              first_derivative_sigmas_ns.plot,
#              second_derivative_sigmas_ns.plot,
#              nrow = 3, ncol = 1)

