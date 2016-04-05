#!/usr/bin/env Rscript

# Treatment rule(s)
alwaysTreated0 <- function(L0){
  1
}

neverTreated0 <- function(L0){
  0
}


# Functions to be used later
logit <- function(x){
  log(x/(1-x))
}

expit<-function(x){
  result<-exp(x)/(1+exp(x))
  result[is.nan(result)] <- 1
  result
}

g_to_g_delta<-function(delta, g){
  (g < delta) * delta + (g >= delta) * g
}

loss <- function(coefficients, outcome, covariates, boolean_subset, offset = 0){
  Q_k <- as.vector(expit(covariates %*% coefficients + offset))
  sum(-boolean_subset * (outcome * log(Q_k) + (1 - outcome) * log(1 - Q_k)))
}

build_univariate_loss <- function(outcome, covariate, boolean_subset, offset = 0){
  univariate_loss <- function(coefficient){
    Q_k <- as.vector(expit(covariate * coefficient + offset))
    sum(-boolean_subset * (outcome * log(Q_k) + (1-outcome) * log(1 - Q_k)))
  }
}

# Data generation
generate_data <- function(type = "L0_unif", positivity_parameter, alpha0, beta0, beta1, beta2, n){
  if(type == "L0_unif") 
    L0 <- runif(n, min = -positivity_parameter, max = positivity_parameter)
  else 
    L0 <- rexp(n, rate = 1 / positivity_parameter) * (1 - 2 * rbinom(n, 1, prob = 0.5))
  L0 <- runif(n, min = -positivity_parameter, max = positivity_parameter)
  g00 <- expit(alpha0*L0)
  A0 <- rbinom(n, 1, g00)
  PL1givenA0L0 <- expit(beta0 + beta1 * A0 + beta2 * L0)
  L1 <- rbinom(n, 1, PL1givenA0L0)
  list(L0 = L0, A0 = A0, L1 = L1)
}

# TMLE of truncated target parameter
TMLE_truncated_target <- function(observed_data, d0, delta, Q_misspecified = F){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  
  # 0. Fit models for g_{n,k=0}
  g0_fit <- glm(A0 ~ 1 + L0, family=binomial)
  g0 <- as.vector(predict(g0_fit, type="response"))
  g1 <- 1 - g0
  PA0givenL0 <- g0 * A0 + g1 * (1 - A0)
  PdL0_givenL0 <- g0 * d0(L0) + g1 * (1 - d0(L0))
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  if(Q_misspecified == FALSE){
#     coeffs_Q1d_bar_0n <- optim(par=c(0,0,0), fn=loss, outcome=L1, 
#                                covariates=cbind(1,L0,A0), 
#                                boolean_subset=(A0==d0(L0)))$par
#     offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0, d0(L0)) %*% coeffs_Q1d_bar_0n)
    coeffs_Q1d_bar_0n <- glm(L1 ~ A0 + L0, family = binomial())$coefficients
#     offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0, d0(L0)) %*% coeffs_Q1d_bar_0n)
    coeffs_Q1d_bar_0n[is.na(coeffs_Q1d_bar_0n)] <- 0
    offset_vals_Q1d_bar_0n <- as.vector(cbind(1, d0(L0), L0) %*% coeffs_Q1d_bar_0n)
  }else{
    offset_vals_Q1d_bar_0n <- rep(logit(mean(L1[A0==d0(L0)])), n)
  }
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Compute clever covariate
  PA0givenL0_delta <- g_to_g_delta(delta, PA0givenL0)
  PdL0_givenL0_delta <- g_to_g_delta(delta, PdL0_givenL0)
  
  H_delta <- (A0 == d0(L0)) / PA0givenL0_delta
  H_delta_setting_A_to_d <- 1 / PdL0_givenL0
  
  # Fit parametric submodel to training set
  epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = logit(Q1d_bar_0n),
                 subset = which(A0 == d0(L0)))$coefficients[1]
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
#   # 1.c Refit model for L1 given A1, A0, L1, L0 with clever cov. and offset:
#   interval_epsilon_1_n <- c(-10,10)
#   univariate_loss_truncated_g <- build_univariate_loss(outcome = L1, covariate = H_delta, 
#                                                        boolean_subset = (A0 == d0(L0)), 
#                                                        offset=offset_vals_Q1_d)
#   epsilon_1_n_truncated_g <- optim(fn = univariate_loss_truncated_g, par=0, method="Brent", lower=interval_epsilon_1_n[1], upper=interval_epsilon_1_n[2])$par
# 
#   Q_1_dn_star_truncated_g <- expit(offset_vals_Q1_d+epsilon_1_n_truncated_g*C1_truncated_g_delta_setting_A_to_d)


  # Return estimator
  mean(PdL0_givenL0 / PdL0_givenL0_delta * Q1d_bar_star_n)
#   mean(Q1d_bar_0n)
}

# debug(TMLE_truncated_target)

IPTW_ATE <- function(observed_data){
    L1 <- observed_data$L1; A0 <- observed_data$A0; L0 <- observed_data$L0
    
    g0_fit <- glm(A0 ~ 1 + L0, family=binomial)
    g0 <- as.vector(predict(g0_fit, type="response"))
    
    g1 <- 1 - g0
    
    mean(A0 * L1/ g0 - (1 - A0) * L1 / g1)
}

# Simulations -------------------------------------------------------------
# set.seed(0)
delta0s <- c(1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 2e-1)
# delta0s <- c(1e-4, 5e-4, (1:9)*1e-3, (1:9)*1e-2, (1:4)*1e-1)
delta0s <- 1e-4
orders <- 0:8
# orders <- 9

# Compute true value of EY^d
compute_Psi_d_MC <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, M){
  # Monte-Carlo estimation of the true value of mean of Yd
  if(type == "L0_unif") 
    L0_MC <- runif(M, min = -positivity_parameter, max = positivity_parameter)
  else 
    L0_MC <- rexp(M, rate = 1 / positivity_parameter) * (1 - 2 * rbinom(M, 1, prob = 0.5))
  A0_MC <- d0(L0_MC)
  g0_MC <- expit(alpha0 * L0_MC)
  A0_MC_bis <- 
  PL1givenA0L0_MC <- expit(beta0 + beta1 * A0_MC + beta2 * L0_MC)
  mean(PL1givenA0L0_MC)
}

# Specify the jobs. A job is the computation of a batch. 
# It is fully characterized by the parameters_tuple_id that the batch corresponds to.
ns <- c((1:9)*100, c(1:9)*1000, 2*c(1:5)*1e4)
# parameters_grid <- rbind(expand.grid(type = "L0_unif", positivity_parameter = c(2, 4), 
#                                      alpha0 = 2, beta0 = -3, beta1 = +1.5, beta2 = 1, n = ns),
#                          expand.grid(type = "L0_exp", positivity_parameter = c(2, 4), 
#                                      alpha0 = 2, beta0 = -3, beta1 = +1.5, beta2 = 1, n = ns))
parameters_grid <- expand.grid(type = "L0_unif", positivity_parameter = 4, 
                                     alpha0 = 2, beta0 = -3, beta1 = +1.5, beta2 = -1, n = ns)

batch_size <- 2; nb_batchs <- 500
jobs <- kronecker(1:nrow(parameters_grid), rep(1, nb_batchs))
first_seed_batch <- 1:length(jobs) * batch_size
jobs_permutation <- sample(1:length(jobs))
jobs <- jobs[jobs_permutation]
first_seed_batch <- first_seed_batch[jobs_permutation]

# Compute target parameter for each parameters tuple id
target_parameters <- vector()
# for(i in 1:nrow(parameters_grid))
i <- 1
  target_parameters[i] <- compute_Psi_d_MC(type = parameters_grid[i, "type"],
                                           positivity_parameter = parameters_grid[i, "positivity_parameter"],
                                           alpha0 = parameters_grid[i, "alpha0"], 
                                           beta0 = parameters_grid[i, "beta0"],
                                           beta1 = parameters_grid[i, "beta1"],
                                           beta2 = parameters_grid[i, "beta2"],
                                           d0 = alwaysTreated0, M = 1e6)

# Save the parameters' grid
write.table(parameters_grid, file = "parameters_grid.csv", append = F, row.names=F, col.names=T,  sep=",")

# Perform the jobs in parallel
# library(Rmpi); library(doMPI)
# 
# cl <- startMPIcluster(72)
# registerDoMPI(cl)
# library(foreach); library(doParallel)
# cl <- makeCluster(getOption("cl.cores", 2), outfile = "")
# registerDoParallel(cl)

# results <- foreach(i = 1:length(jobs)) %dopar% { #job is a parameter_tuple_idS
# # for(i in 1:length(jobs)){
#   job <- jobs[i]
#   results_batch <- matrix(0, nrow = batch_size, ncol = 4)
#   colnames(results_batch) <- c("parameters_tuple_id", "EYd", "seed","Utgtd-untr")
#   for(j in 1:batch_size){
#     seed <- first_seed_batch[i] + j - 1; set.seed(seed)
#     observed_data <- generate_data(type = parameters_grid[job,]$type, 
#                                    positivity_parameter = parameters_grid[job,]$positivity_parameter, 
#                                    alpha0 = parameters_grid[job,]$alpha0,
#                                    beta0 = parameters_grid[job,]$beta0, beta1 = parameters_grid[job,]$beta1, 
#                                    beta2 = parameters_grid[job,]$beta2, n = parameters_grid[job,]$n)
#     
#     Utgtd_untr_Psi_n <- TMLE_truncated_target(observed_data, alwaysTreated0, 0, Q_misspecified = F)
# #     print(Utgtd_untr_Psi_n)
#     results_batch[j, ] <- c(job, target_parameters[job], seed, Utgtd_untr_Psi_n)
#   }
# #   print(results_batch)
#   
#   if(!file.exists("truncated_TMLE_results.csv")){
#     write.table(results_batch, file="truncated_TMLE_results.csv", append=T, row.names=F, col.names=T,  sep=",")
#   }else{
#     write.table(results_batch, file="truncated_TMLE_results.csv", append=T, row.names=F, col.names=F,  sep=",")
#   }
#   results_batch
# }
# 
# save(results, parameters_grid, file = "truncated_TMLE_results")
# 
# closeCluster(cl)
# mpi.quit()

job <- 23

observed_data <- generate_data(type = 'L0_unif', 
                               positivity_parameter = 4, 
                               alpha0 = parameters_grid[job,]$alpha0,
                               beta0 = parameters_grid[job,]$beta0, beta1 = parameters_grid[job,]$beta1, 
                               beta2 = parameters_grid[job,]$beta2, n = parameters_grid[job,]$n)

Utgtd_untr_Psi_n <- TMLE_truncated_target(observed_data, alwaysTreated0, 0, Q_misspecified = F)
ATE <- TMLE_truncated_target(observed_data, alwaysTreated0, 0, Q_misspecified = F) -
  TMLE_truncated_target(observed_data, neverTreated0, 0, Q_misspecified = F)

print(target_parameters[1])
print(Utgtd_untr_Psi_n)
print(ATE)
print(IPTW_ATE(observed_data))
print(tmle(observed_data$L1, observed_data$A0, as.matrix(observed_data$L0)))
