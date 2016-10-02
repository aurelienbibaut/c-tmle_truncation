source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)
library(data.tree)

# Specify data-generating distribution ------------------------------------

# Set of parameters 1
# lambda <- 2; alpha0 <- 2; beta0 <- -3; beta1 <- 1.5; beta2 <- 1
# beta <- 0.25
# gamma <- 0.125
# kappa <- 1 / (2 * (gamma + 1 - beta))
# Set of parameters 2
# lambda <- 2; alpha0 <- 4; beta0 <- -3; beta1 <- 1.5; beta2 <- 0.5
# kappa <- 5 / 4
# beta <- 2 - kappa
# gamma <- 1 - kappa / 2
# Set of parameters 3
lambda <- 2; alpha0 <- 4; beta2 <- -3; beta0 <- -1; beta1 <- 1
beta <- 7/8; gamma <- 1/16
true_rate <- 1 / (2 * (gamma + 1 - beta))

# Define candidate rates, etas and base constant C0 in delta_n
candidate_rates <- seq(from = 0.1, to = 2.5, length = 10)
etas <- c(1.1, 1.5, 2, 3)
C0 <- 0.05; n0 <- 100

# Define estimation tasks
n_full <- 1e7
deltas <- 10^seq(from = -10, to = -1, length = 10)
estimation_tasks <- expand.grid(eta = etas, candidate_rate = candidate_rates)
jobs <- 1:nrow(estimation_tasks)

# Define finite difference function
compute_fin_diff <- function(data, indices, delta, Delta){
  replicate.data <- data[indices, ]
  
  # Compute finite differnece
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(replicate.data, delta + Delta)$Psi_n
  Psi_n_delta <- TMLE_EY1_speedglm(replicate.data, delta)$Psi_n
  
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Deltas
}

# Generate one sample
observed_data.list <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, n_full)
observed_data <- data.frame(L0 = observed_data.list$L0,
                            A0 = observed_data.list$A0,
                            L1 = observed_data.list$L1)
# Perform tasks
# Set up cluster
cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
registerDoParallel(cl)

results <- foreach(i=jobs, .combine = rbind, 
                   .packages = c("speedglm"), .verbose = T, .inorder = T) %dopar% {
                     
                     
                     eta <- estimation_tasks[i, ]$eta
                     delta_n_full <- (C0 * (n_full / n0)^(-estimation_tasks[i, ]$candidate_rate))^(1 / eta)
                     delta_n_subsample <- (C0 * (n_subsample / n0)^(-estimation_tasks[i, ]$candidate_rate))^(1 / eta)
                     Delta_n_full <- n_full^(-0.25) * delta_n_full^((beta + 1 - gamma) / 2)
                     Delta_n_subsample <- n_subsample^(-0.25) * delta_n_subsample^((beta + 1 - gamma) / 2)
                     
                     full_sample.pivot <- compute_pivot(observed_data, 1:n_full, delta_n_full, Delta_n_full, eta, n_full)
                     subsambles.pivots <- apply(subsamples.indices, 1, 
                                                function(x) compute_pivot(observed_data, x, delta_n_subsample, 
                                                                          Delta_n_subsample, eta, n_subsample))
                     
                     c(eta = estimation_tasks[i, ]$eta, 
                       candidate_rate = estimation_tasks[i, ]$candidate_rate,
                       full_sample.pivot = full_sample.pivot,
                       subsample.pivot = median(subsambles.pivots))
                   }


# Make plots
pivots.plots <- list()
for(i in 1:length(etas)){
  pivots.plots[[i]] <- ggplot(subset(as.data.frame(results), eta == etas[i]),
                              aes(x = 0, xend = 1, y = subsample.pivot, yend = full_sample.pivot,
                                  colour = factor(candidate_rate))) + geom_segment() +
    ggtitle(paste("eta =", etas[i], sep = ''))
}
