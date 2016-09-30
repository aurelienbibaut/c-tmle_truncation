source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(speedglm); library(ggplot2)
library(foreach); library(doParallel)

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

# etas <- c(2.5, 3, 3.5, 4, 4.5, 5)
etas <- c(2.5, 3, 4, 5)
finite_diffs <- vector(); finite_diffs_bias <- vector(); true_finite_diffs_bias <- vector()
ns <- floor(10^seq(from = 2, to = 6.5, by = 0.5))
deltas <- vector()

observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, max(ns))
# data <- data.frame(L0 = observed_data$L0, A0 = observed_data$A0, L1 = observed_data$L1)



jobs <- expand.grid(n = ns, eta = etas)

# Set up cluster
cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
registerDoParallel(cl)

results <- foreach(i=1:nrow(jobs), .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     n <- jobs[i, ]$n; eta <- jobs[i, ]$eta
                     
                     delta_n_plus <- n^(-1 / (2 * eta * (gamma + 1 - beta)))
                     delta_n_a <- delta_n_plus * n^(-0.05)
                     delta_n_b <- delta_n_plus * n^0.05
                     
                     Delta <- n^(-0.25) * delta_n_plus^((beta + 1 - gamma) / 2)
                     Delta_a <- n^(-0.25) * delta_n_a^((beta + 1 - gamma) / 2)
                     Delta_b <- n^(-0.25) * delta_n_b^((beta + 1 - gamma) / 2)
                     
                     cat('Job ', i, ', n = ', n, ', eta = ', eta, ', delta =', delta_n_plus, ', Delta = ', Delta, '\n')
                     Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_plus + Delta)$Psi_n
                     Psi_n_delta_plus_Delta_a <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_a + Delta_a)$Psi_n
                     Psi_n_delta_plus_Delta_b <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_b + Delta_b)$Psi_n
                     cat('Job ', i, ', Psi_n(delta + Delta) = ', Psi_n_delta_plus_Delta, '\n')
                     
                     Psi_n_delta <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_plus)$Psi_n
                     Psi_n_delta_a <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_a)$Psi_n
                     Psi_n_delta_b <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_b)$Psi_n
                     
                     cat('Job ', i,', Psi_n(delta) = ', Psi_n_delta, '\n')
                     
                     cat('Job ', i, ',Psi_n(', delta_n_plus, ') = ', Psi_n_delta)
                     
                     
                     finite_diffs <- (Psi_n_delta_plus_Delta - Psi_n_delta + 
                                        n^(-0.5) * ((delta_n_plus + Delta)^(-gamma) - delta_n_plus^(-gamma))) / Delta
                     
                     
                     finite_diffs_bias <- (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
                     finite_diffs_bias_a <- (Psi_n_delta_plus_Delta_a - Psi_n_delta_a) / Delta_a
                     finite_diffs_bias_b <- (Psi_n_delta_plus_Delta_b - Psi_n_delta_b) / Delta_b
                     
                     # True quantities
                     
                     #                      Psi_0_delta_plus_Delta <- compute_true_Psi0_delta("L0_exp", lambda, alpha0, beta0, beta1, beta2,
                     #                                                                        alwaysTreated0, delta_n_plus + Delta)
                     #                      Psi_0_delta <- compute_true_Psi0_delta("L0_exp", lambda, alpha0, beta0, beta1, beta2, 
                     #                                                             alwaysTreated0, delta_n_plus)
                     #                      
                     #                      true_finite_diffs_bias <- (Psi_0_delta_plus_Delta - Psi_0_delta) / Delta
                     
                     #                      true_var_IC_delta <- true_var_IC.result$var_IC0
                     #                      true_var.first_comp <- true_var_IC.result$first_component
                     #                      true_var.first_comp.a <- true_var_IC.result$first_component.a
                     
                     #                      c(eta, n, finite_diffs, 
                     #                        finite_diffs_bias,
                     #                        finite_diffs_bias_a,
                     #                        finite_diffs_bias_b,
                     #                        true_finite_diffs_bias)
                     
                     iteration.result <- rbind(c(eta, n, delta_n_plus, finite_diffs_bias, 'optimal_rate'),
                                               c(eta, n, delta_n_a, finite_diffs_bias_a, 'too_slow_rate'),
                                               c(eta, n, delta_n_b, finite_diffs_bias_b, 'too_fast_rate'))
                   }

stopCluster(cl)

row.names(results) <- NULL
colnames(results) <- c('eta', 'n', 'delta', 'fin_diff', 'rate.type')
# results_df <- as.data.frame(rbind(cbind(results[, c('eta', 'n', 'delta')], fin_diff = results[, 'bias.fin_diff'], type = 'bias.fin_diff'),
#                                   cbind(results[, c('eta', 'n', 'delta')], fin_diff = results[, 'bias.fin_diff.a'], type = 'bias.fin_diff.a'),
#                                   cbind(results[, c('eta', 'n', 'delta')], fin_diff = results[, 'bias.fin_diff.b'], type = 'bias.fin_diff.b'),
#                                   cbind(results[, c('eta', 'n', 'delta')], fin_diff = results[, 'true_bias.fin_diff'], type = 'true_bias.fin_diff')))
results_df <- as.data.frame(results)
results_df <- transform(results_df, n = as.numeric(as.character(n)), 
                        eta = as.numeric(as.character(eta)),
                        delta = as.numeric(as.character(delta)),
                        fin_diff = as.numeric(as.character(fin_diff)))

results_df <- cbind(results_df, track = apply(cbind(results_df$eta, 
                                                    as.character(results_df$rate.type)), 1, function(x) paste(c(x[1], x[2]), 
                                                                                                              collapse = '')))
all_etas.plot_bis <- ggplot(results_df, 
                            aes(x = log(delta), y = eta * log(abs(fin_diff) * delta^2) + log(n), 
                                group = factor(track),
                                colour = factor(eta))) + 
  geom_line()

fin_diffs.all_etas.plot_bis <- ggplot(results_df, 
                                      aes(x = log(delta), y = log(abs(fin_diff)), 
                                          group = factor(track),
                                          colour = factor(eta), label = round(log(n)/log(10), 1))) + 
  geom_line() + geom_text() +
  geom_abline(intercept = -4, slope = -beta) +
  geom_abline(intercept = -4.1, slope = -beta) +
  geom_abline(intercept = -4.2, slope = -beta)
# all_etas.sigmas.plot <- ggplot(subset(results_df, n < 1e8 & eta > 0), aes(x = log(delta), y = log(var_IC_delta), label = round(log(n) / log(10), 1), colour = factor(eta))) +
#   geom_line() + geom_point() + geom_text() + 
#   geom_line(aes(x = log(delta), y = log(true_var_IC_delta) + 0.3, colour = factor(eta), size = 1.1)) + 
#   geom_line(aes(x = log(delta), y = log(true_var.first_comp), colour = factor(eta), linetype = '12345678')) +
#   geom_line(aes(x = log(delta), y = log(true_var.first_comp.a) + 0.8, colour = factor(eta), linetype = 'dotted')) +
#   geom_abline(intercept = -3, slope = -2 * gamma) +
#   geom_abline(intercept = -3.5, slope = -2 * gamma) +
#   geom_abline(intercept = -4, slope = -2 * gamma)
# 
# print(all_etas.sigmas.plot)
# plot(log(deltas), log(abs(finite_diffs_bias)), ylim = c(min(c(log(abs(finite_diffs_bias)), log(abs(true_finite_diffs_bias)))),
#                                                         max(c(log(abs(finite_diffs_bias)), log(abs(true_finite_diffs_bias))))),
#      main= substitute(list(eta) == list(x),
#                       list(x = eta)))
# lines(log(deltas), log(abs(true_finite_diffs_bias)))
