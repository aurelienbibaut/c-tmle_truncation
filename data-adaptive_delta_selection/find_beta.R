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

# Define finite difference function
finite_difference <- function(data, delta, n){
  Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(data, delta + Delta)$Psi_n
  Psi_n_delta <- TMLE_EY1_speedglm(data, delta)$Psi_n
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}


# etas <- c(2.5, 3, 3.5, 4, 4.5, 5)
etas <- seq(from = 2, to = 5, by = 0.5)
finite_diffs <- vector(); finite_diffs_bias <- vector(); true_finite_diffs_bias <- vector()
ns <- floor(10^seq(from = 2, to = 4, by = 0.2))
nb_observations <- 1e4

delta_min <- 1e-5

observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, nb_observations)

jobs <- expand.grid(n = ns, eta = etas)

# Set up cluster
cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
registerDoParallel(cl)

results <- foreach(i=1:nrow(jobs), .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     n <- jobs[i, ]$n; eta <- jobs[i, ]$eta
                     delta <- n^(-1 / (2 * eta * (gamma + 1 - beta)))
                     
                     if(n == nb_observations){
                       indices <- matrix(1:nb_observations, nrow = 1)
                     }else if(n < nb_observations){
                       indices <- t(replicate(1000, sample(1:nb_observations, n, replace = F)))
                     }else{
                       indices <- t(replicate(1000, sample(1:nb_observations, n, replace = T)))
                     }
                     # indices <- matrix(1:n, nrow = 1)
                     
                     fin_diffs <- apply(indices, 1, function(y) finite_difference(lapply(observed_data, function(x) x[y]), delta, n))
                     fin_diff <- median(fin_diffs)
                     
                     cat('Job ', i, ', n = ', n, ', eta = ', eta, 
                         ', delta = ', delta, ', fin_diff = ', fin_diff, '\n')
                     c(n, eta, delta, fin_diff)
                   }


stopCluster(cl)

row.names(results) <- NULL
colnames(results) <- c('n', 'eta', 'delta', 'fin_diff')

results_df <- as.data.frame(results)
results_df <- transform(results_df, n = as.numeric(as.character(n)), 
                        eta = as.numeric(as.character(eta)),
                        delta = as.numeric(as.character(delta)),
                        fin_diff = as.numeric(as.character(fin_diff)))

# results_df <- cbind(results_df, track = apply(cbind(results_df$eta, 
#                                                     as.character(results_df$rate.type)), 1, function(x) paste(c(x[1], x[2]), 
#                                                                                                               collapse = '')))
# all_etas.plot <- ggplot(results_df, 
#                         aes(x = log(delta) / log(10), y = eta * log(abs(fin_diff) * delta^2) + log(n), 
#                             group = factor(track),
#                             colour = factor(eta))) + 
#   geom_line()

# print(all_etas.plot)

regression_df <- as.data.frame(cbind(log_fin_diff = log(abs(results_df$fin_diff))/log(10), 
                                     log_delta = log(results_df$delta) / log(10)))
line_fit <- lm(formula = log_fin_diff ~ log_delta, regression_df)
lts.line_fit <- ltsReg(regression_df$log_delta, regression_df$log_fin_diff)

fin_diffs.all_etas.plot <- ggplot(results_df, 
                                  aes(x = log(delta) / log(10), 
                                      y = log(abs(fin_diff)) / log(10),
                                      colour = factor(eta),
                                      label = round(log(n) / log(10), 1))) + 
  geom_line() + geom_text() +
  geom_abline(intercept = -2, slope = -beta, linetype = 'dashed') +
  geom_abline(intercept = -2.1, slope = -beta, linetype = 'dashed') +
  geom_abline(intercept = -2.2, slope = -beta, linetype = 'dashed') +
  geom_abline(intercept = line_fit$coefficients["(Intercept)"],
              slope = line_fit$coefficients["log_delta"], colour = 'blue', size = 1) +
  geom_abline(intercept = lts.line_fit$raw.coefficients[1],
              slope = lts.line_fit$raw.coefficients[2], colour = 'green', size = 1) +
  ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                       group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                     list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))
print(fin_diffs.all_etas.plot)