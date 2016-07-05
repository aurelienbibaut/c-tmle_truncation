source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)
library(boot)


# Specify data-generating distribution ------------------------------------

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


# Define inference functions ----------------------------------------------

# Define finite difference functions
finite_difference_with_bootstrap_of_targeting_step <- function(data, delta, Delta){
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(data, delta + Delta)$Psi_n
  TMLE_Psi_0_delta.bootstrap <- TMLE_EY1_bootstrap_speedglm(data, delta, nb_boostrap_samples = 1000)
  Psi_n_delta <- TMLE_Psi_0_delta.bootstrap$Psi_n
  list(fin_diff = (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta,
       shapiro.p_value = TMLE_Psi_0_delta.bootstrap$shapiro.p_value)
}

finite_difference_for_boot <- function(data, indices, delta, Delta){
  replicate.data <- data[indices, ]
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(replicate.data, delta + Delta)$Psi_n
  Psi_n_delta <- TMLE_EY1_speedglm(replicate.data, delta)$Psi_n
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}

# Set up cluster
cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
registerDoParallel(cl)

# Simulate data
n <- 1e4
observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, n)

# Set up tasks
deltas <- 10^seq(from = -5, to = -0.8, by = 0.1)
Delta.delta_rates <- c(0.8, 1, 1.1, 1.375, 1.5)
jobs <- expand.grid(delta = deltas, Delta.delta_rate = Delta.delta_rates)

results <- foreach(i=1:nrow(jobs), .combine = rbind,
                   .packages = c('speedglm', 'boot'), .verbose = T, .inorder = T) %dopar% {
                     
                     delta <- jobs[i, ]$delta; Delta <- n^(-0.25) * delta^jobs[i, ]$Delta.delta_rate
                     
                     #                      fin_diff.result <- finite_difference(observed_data, delta, Delta)
                     observed_data.df <- data.frame(L0 = observed_data$L0,
                                                    A0 = observed_data$A0,
                                                    L1 = observed_data$L1)
                     fin_diff <- finite_difference_for_boot(observed_data.df, 1:n, delta, Delta)
                     fin_diff.bootstrap <- boot(data = observed_data.df,
                                                statistic = finite_difference_for_boot,
                                                R = 1000, sim = 'ordinary',
                                                delta = delta, Delta = Delta)$t
                     shapiro.p_value <- 0
                     try(shapiro.p_value <- shapiro.test(fin_diff.bootstrap)$p.value)
                     if(is.null(shapiro.p_value)) shapiro.p_value <- 0
                     
                     c(delta, jobs[i, ]$Delta.delta_rate, fin_diff, shapiro.p_value)
                   }
stopCluster(cl)

row.names(results) <- NULL
colnames(results) <- c('delta', 'Delta.delta_rate', 'fin_diff', 'p_value')

results_df <- as.data.frame(results)

# Make plot of finite differences
fin_diffs.plot <- ggplot(results_df, 
                         aes(x = log(delta) / log(10), 
                             y = log(abs(fin_diff)) / log(10),
                             colour = factor(Delta.delta_rate))) + 
  geom_line() + geom_point(aes(size = log(p_value) / log(10))) + 
  geom_abline(intercept = -0.5, slope = -beta, linetype = 'dotted') +
  geom_abline(intercept = -1, slope = -beta, linetype = 'dotted') +
  geom_abline(intercept = -1.5, slope = -beta, linetype = 'dotted') +
  xlab(expression(log[10](delta))) + 
  ylab(expression(log[10](widehat(Delta * b[n])(delta)))) +
  ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                       group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                     list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))

print(fin_diffs.plot)