source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)

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
# Set of parameters 3
# lambda <- 2; alpha0 <- 4; beta2 <- -3; beta0 <- -1; beta1 <- 1
# beta <- 7/8; gamma <- 1/16
true_rate <- 1 / (2 * (gamma + 1 - beta))

# Define candidate rates, etas and base constant C0 in delta_n
candidate_rates <- seq(from = 0.1, to = 2.5, length = 10)
etas <- c(1.1, 1.5, 2, 3)
C0 <- 0.05; n0 <- 100

# Define estimation tasks
n <- floor(10^7)
deltas <- 10^seq(from = -7, to = -1, length = 10)

# Generate one sample
observed_data.list <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, n)
observed_data <- data.frame(L0 = observed_data.list$L0,
                            A0 = observed_data.list$A0,
                            L1 = observed_data.list$L1)
# Perform tasks
# Set up cluster
cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
registerDoParallel(cl)

results <- foreach(delta = deltas, .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
                     
                     Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(observed_data, delta + Delta)$Psi_n
                     TMLE_delta <- TMLE_EY1_speedglm(observed_data, delta)
                     
                     fin_diff <- (Psi_n_delta_plus_Delta - TMLE_delta$Psi_n) / Delta
                     
                     c(delta = delta, fin_diff = fin_diff, sigma_n = sqrt(TMLE_delta$var_IC.first_part))
                   }

colnames(results)[1] <- "delta"
fin_diff_plot <- ggplot(as.data.frame(results), aes(x = log(delta) / log(10),
                                                    y = log(abs(fin_diff)) / log(10))) +
  geom_point() + geom_line() + geom_abline(intercept = 0, slope = -beta)

sigma_n_plot <- ggplot(as.data.frame(results), aes(x = log(delta) / log(10),
                                                    y = log(sigma_n) / log(10))) +
  geom_point() + geom_line() + geom_abline(intercept = -0.25, slope = -gamma)

grid.arrange(nrow = 2, ncol = 1, fin_diff_plot, sigma_n_plot)

