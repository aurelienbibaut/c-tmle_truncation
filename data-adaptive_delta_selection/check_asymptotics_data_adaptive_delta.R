source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(speedglm)
library(foreach); library(doParallel)

# debug(TMLE_EY1)
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

etas <- c(3, 3.5, 4, 4.5, 5)
finite_diffs <- vector(); finite_diffs_bias <- vector(); true_finite_diffs_bias <- vector()
ns <- floor(10^seq(from = 2, to = 7, by = 0.5))
deltas <- vector()

observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, max(ns))
# data <- data.frame(L0 = observed_data$L0, A0 = observed_data$A0, L1 = observed_data$L1)



jobs <- expand.grid(n = ns, eta = etas)

# Set up cluster
cl <- makeCluster(getOption("cl.cores", 5), outfile = '')
registerDoParallel(cl)

results <- foreach(i=1:nrow(jobs), .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     n <- jobs[i, ]$n; eta <- jobs[i, ]$eta
                     
                     delta_n_plus <- n^(-1 / (2 * eta * (gamma + 1 - beta)))
                     # delta_n_a <- delta_n_plus * n^(-0.05)
                     # delta_n_b <- delta_n_plus * n^0.05
                     
                     Delta <- n^(-0.25) * delta_n_plus^((beta + 1 - gamma) / 2)
                     # Delta_a <- n^(-0.25) * delta_n_a^((beta + 1 - gamma) / 2)
                     # Delta_b <- n^(-0.25) * delta_n_b^((beta + 1 - gamma) / 2)
                     
                     cat('Job ', i, ', n = ', n, ', eta = ', eta, ', delta =', delta_n_plus, ', Delta = ', Delta, '\n')
                     TMLE_delta_plus_Delta <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_plus + Delta)
                     Psi_n_delta_plus_Delta <- TMLE_delta_plus_Delta$Psi_n
                     # Psi_n_delta_plus_Delta_a <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_a + Delta_a)
                     # Psi_n_delta_plus_Delta_b <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_b + Delta_b)
                     cat('Job ', i, ', Psi_n(delta + Delta) = ', Psi_n_delta_plus_Delta, '\n')
                     TMLE_delta <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_plus)
                     Psi_n_delta <- TMLE_delta$Psi_n
                     sigma_n_delta <- TMLE_delta$sigma_n
                     cat('Job ', i,'Psi(delta)', Psi_n_delta, '\n')
                     # Psi_n_delta_a <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_a)
                     # Psi_n_delta_b <- TMLE_EY1_speedglm(lapply(observed_data, function(x) x[1:n]), delta_n_b)
                     cat('Job ', i, ',Psi_n(', delta_n_plus, ') = ', Psi_n_delta, 
                         'sigma_n(', delta_n_plus, ') = ', sigma_n_delta, '\n')
                     
                     
                     finite_diffs <- (Psi_n_delta_plus_Delta - Psi_n_delta + 
                                        n^(-0.5) * ((delta_n_plus + Delta)^(-gamma) - delta_n_plus^(-gamma))) / Delta
                     
                     
                     finite_diffs_bias <- (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
                     
                     # finite_diffs_bias_a <- (Psi_n_delta_plus_Delta_a - Psi_n_delta_a) / Delta_a
                     # finite_diffs_bias_b <- (Psi_n_delta_plus_Delta_b - Psi_n_delta_b) / Delta_b
                     finite_diffs_bias_a <- finite_diffs_bias
                     finite_diffs_bias_b <- finite_diffs_bias
                     
                     Psi_0_delta_plus_Delta <- compute_true_Psi0_delta("L0_exp", 2, 2, -3, 1.5, 1, 
                                                                       alwaysTreated0, delta_n_plus + Delta)
                     Psi_0_delta <- compute_true_Psi0_delta("L0_exp", 2, 2, -3, 1.5, 1, 
                                                            alwaysTreated0, delta_n_plus)
                     
                     true_finite_diffs_bias <- (Psi_0_delta_plus_Delta - Psi_0_delta) / Delta
                     
                     
                     #                      c(eta, n, finite_diffs, 
                     #                        finite_diffs_bias,
                     #                        finite_diffs_bias_a,
                     #                        finite_diffs_bias_b,
                     #                        true_finite_diffs_bias)
                     c(eta, n, finite_diffs_bias, sigma_n_delta)
                   }

stopCluster(cl)

row.names(results) <- NULL
# colnames(results) <- c('eta', 'n', 'fin_diff',
#                        'bias.fin_diff',
#                        'bias.fin_diff.a',
#                        'bias.fin_diff.b',
#                        'true_bias.fin_diff')
colnames(results) <- c('eta', 'bias.fin_diff', 'sigma_n')
results <- cbind(results, delta = results[, 'n']^(- 1 / (2 * results[, 'eta'] * (gamma + 1 - beta))))

results_df <- as.data.frame(rbind(cbind(results[, c('eta', 'n', 'delta')], fin_diff = results[, 'bias.fin_diff'], type = 'bias.fin_diff'),
                                  cbind(results[, c('eta', 'n', 'delta')], fin_diff = results[, 'bias.fin_diff.a'], type = 'bias.fin_diff.a'),
                                  cbind(results[, c('eta', 'n', 'delta')], fin_diff = results[, 'bias.fin_diff.b'], type = 'bias.fin_diff.b'),
                                  cbind(results[, c('eta', 'n', 'delta')], fin_diff = results[, 'true_bias.fin_diff'], type = 'true_bias.fin_diff')))
results_df <- transform(results_df, n = as.numeric(as.character(n)), 
                        eta = as.numeric(as.character(eta)),
                        delta = as.numeric(as.character(delta)),
                        fin_diff = as.numeric(as.character(fin_diff)))


# plots <- list()
# for(i in 1:length(etas)){
#   data <- subset(results_df, eta == etas[i] & type == 'bias.fin_diff')
#   x1 <- min(log(data$delta)); x2 <- 0
#   y1 <- log(abs(data$fin_diff[data$delta == min(data$delta)]))
#   y2 <- y1 - 0.25 * (x2 - x1)
#   df <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)
#   plots[[i]] <- ggplot(data = data, aes(x = log(delta), y = log(abs(fin_diff)), colour = type)) + 
#     geom_line() +
#     geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df) + 
#     ggtitle(substitute(list(eta) == list(x),
#                        list(x = etas[i])))
# }

all_etas.plot <- ggplot(subset(results_df, type %in% c('bias.fin_diff', 'bias.fin_diff.a', 'bias.fin_diff.b') & eta > 2.5), 
                        aes(x = log(delta), y = log(abs(fin_diff)), colour =factor(eta))) + 
  geom_line() +
  geom_abline(intercept = -4, slope = -beta) +
  geom_abline(intercept = -4.1, slope = -beta) +
  geom_abline(intercept = -4.2, slope = -beta)

results_df <- cbind(results_df, track = apply(cbind(results_df$eta, 
                                                    as.character(results_df$type)), 1, function(x) paste(c(x[1], x[2]), 
                                                                                                         collapse = '')))
all_etas.plot_bis <- ggplot(subset(results_df, type %in% c('bias.fin_diff', 'bias.fin_diff.a', 'bias.fin_diff.b') & eta > 2.5), 
                            aes(x = log(delta), y = eta * log(abs(fin_diff) * delta^2) + log(n), 
                                group = factor(track),
                                colour = factor(eta))) + 
  geom_line() +
  geom_abline(intercept = -4, slope = -beta) +
  geom_abline(intercept = -4.1, slope = -beta) +
  geom_abline(intercept = -4.2, slope = -beta)
# plot(log(deltas), log(abs(finite_diffs_bias)), ylim = c(min(c(log(abs(finite_diffs_bias)), log(abs(true_finite_diffs_bias)))),
#                                                         max(c(log(abs(finite_diffs_bias)), log(abs(true_finite_diffs_bias))))),
#      main= substitute(list(eta) == list(x),
#                       list(x = eta)))
# lines(log(deltas), log(abs(true_finite_diffs_bias)))
