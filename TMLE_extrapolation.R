# Treatment rule(s)
alwaysTreated0 <- function(L0){
  1
}

# Functions to be used later
logit <- function(x){
  log(x/(1-x))
}

expit<-function(x){
  result <- exp(x)/(1+exp(x))
  result[is.nan(result)] <- 1
  result
}

g_to_g_delta<-function(delta, g){
  (g<delta) * delta + (g>=delta) * g
}

loss <- function(coefficients, outcome, covariates, boolean_subset, offset=0){
  Q_k <- as.vector(expit(covariates %*% coefficients + offset))
  sum(-boolean_subset * (outcome * log(Q_k) + (1 - outcome) * log(1 - Q_k)))
}

# Data generation
generate_data <- function(type = "L0_unif", positivity_parameter, alpha0, beta0, beta1, beta2, n){
  if(type == "L0_unif") 
    L0 <- runif(n, min= -positivity_parameter, max= positivity_parameter)
  else 
    L0 <- rexp(n, rate = 1 / positivity_parameter) * (1 - 2 * rbinom(n, 1, prob = 0.5))
  L0 <- runif(n, min = -positivity_parameter, max = positivity_parameter)
  g00 <- expit(alpha0*L0)
  A0 <- rbinom(n, 1, g00)
  PL1givenA0L0 <- expit(beta0+beta1*A0+beta2*L0)
  L1 <- rbinom(n, 1, PL1givenA0L0)
  list(L0 = L0, A0 = A0, L1 = L1)
}

# Compute a(delta0) (as defined in write up)
compute_a_delta0 <- function(delta0, order, n_points = 9, diff_step = NULL, verbose = F){
  
  #   cat("delta0 = ", delta0, " and order = ", order, "\n")
  
  if(order <= 0) return(list(a_delta0 = 1, deltas = delta0))
  
  if(n_points %% 2 == 0) n_points <- n_points + 1
  
  if(is.null(diff_step)){
    if(delta0 - (n_points-1)/2*1e-3 > 0){ diff_step=1e-3 }else{ diff_step = delta0/(n_points-1) }
  }
  bw <- diff_step * 2
  deltas <- delta0 + (1:n_points-1-(n_points-1) / 2)*diff_step
  weights <- exp(-(deltas-delta0)^2 / (2*bw^2)) / sqrt(2*pi*bw^2)
  
  X <- outer(deltas - delta0, 0:order, "^")
  
  A <- apply(diag(nrow(X)), 2, function(C) lm.wfit(X, C, weights)$coefficients)
  a_delta0 <- (-delta0)^(0:order) %*% A
  list(a_delta0 = a_delta0, differentiator = A, deltas = deltas)
}

# debug(compute_a_delta0)

# TMLE of truncated target parameter
TMLE_truncated_target <- function(observed_data, d0, delta, Q_misspecified = F){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  
  # 0. Fit models for g_{n,k=0}
  initial_model_for_A0 <- glm(A0 ~ 1 + L0, family=binomial)
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- as.vector(predict(initial_model_for_A0, type="response"))
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  if(Q_misspecified == FALSE){
    coeffs_Q1d_bar_0n <- optim(par=c(0,0,0), fn=loss, outcome=L1, 
                               covariates=cbind(1,L0,A0), 
                               boolean_subset = (A0==d0(L0)))$par
    offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0, d0(L0)) %*% coeffs_Q1d_bar_0n)
  }else{
    offset_vals_Q1d_bar_0n <- rep(logit(mean(L1[A0==d0(L0)])), n)
  }
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Compute clever covariate
  gn0_delta <- g_to_g_delta(delta, gn0)
  H_delta <- (A0 == d0(L0)) / gn0_delta
  H_delta_setting_A_to_d <- 1 / gn0_delta
  
  # Fit parametric submodel to training set
  epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = offset_vals_Q1d_bar_0n,
                 subset = which(A0 == d0(L0)))$coefficients[1]
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Return estimator
  mean(gn0 / gn0_delta * Q1d_bar_star_n)
}

# Untargeted extrapolation
TMLE_extrapolation_bis <- function(observed_data, d0, order, delta0, Q_misspecified = F,
                                   n_points = 11, diff_step = NULL){
  # Compute a_delta0 as defined in write-up
  result_compute_a_delta0 <- compute_a_delta0(delta0, order = order, n_points, diff_step, verbose=F)
  a_delta0 <- result_compute_a_delta0$a_delta0
  deltas <- result_compute_a_delta0$deltas
  
  # Get the Psi(delta) for each deltas
  Psi_deltas <- sapply(deltas, function(delta) TMLE_truncated_target(observed_data, d0, delta, Q_misspecified = F))
  
  # Extrapolate
  a_delta0 %*% Psi_deltas
}

# TMLE of extrapolation
TMLE_extrapolation <- function(observed_data, training_set, test_set, d0, order, delta0, Q_misspecified = F, 
                               n_points = 11, diff_step = NULL){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  
  # Fit models for g_{n,k=0}
  initial_model_for_A0 <- glm(A0 ~ 1 + L0, family=binomial, subset = training_set)
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- as.vector(predict(initial_model_for_A0, type="response"))
  
  # Compute a_delta0 as defined in write-up
  result_compute_a_delta0 <- compute_a_delta0(delta0, order = order, n_points, diff_step, verbose=F)
  a_delta0 <- result_compute_a_delta0$a_delta0
  deltas <- result_compute_a_delta0$deltas
  
  # Compute clever covariate
  gn0_delta0 <- g_to_g_delta(delta0, gn0)
  gn0_deltas <- sapply(deltas, g_to_g_delta, g=gn0)
  H_delta <- (A0==d0(L0)) / gn0 * (outer(gn0, rep(1, length(deltas))) / gn0_deltas) %*% t(a_delta0)
  H_delta_setting_A_to_d <- 1 / gn0 * (outer(gn0, rep(1, length(deltas))) / gn0_deltas) %*% t(a_delta0)
  
  # Fit initial model Q^1_{d,n} of Q^1_{0,d}
  if(Q_misspecified == FALSE){
    coeffs_Q1d_bar_0n <- optim(par = c(0,0,0), fn=loss, outcome=L1, 
                               covariates = cbind(1,L0,A0), 
                               boolean_subset = intersect(which(A0 == d0(L0)), training_set))$par
    offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0, d0(L0)) %*% coeffs_Q1d_bar_0n)
  }else{
    offset_vals_Q1d_bar_0n <- rep(logit(mean(L1[intersect(which(A0 == d0(L0)), training_set)])), n)
  }
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Fit parametric submodel to training set
  epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = offset_vals_Q1d_bar_0n,
                 subset = intersect(which(A0 == d0(L0)), training_set))$coefficients[1]
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Compute estimator and influence curve
  Psi_n <- mean(Q1d_bar_star_n * (outer(gn0, rep(1, length(deltas))) / gn0_deltas) %*% t(a_delta0))
  D_star_n <- H_delta + (1 / gn0_deltas) %*% t(a_delta0) * gn0 * Q1d_bar_star_n #It's actually D_star_n plus its mean
  var_D_star_n <- var(D_star_n)
  var_D_star_n_test <- var(D_star_n[test_set])
  
  # Return estimator, and Q_bar
  list(Psi_n = Psi_n, Q_bar = Q1d_bar_star_n, var_D_star_n = var_D_star_n, var_D_star_n_test = var_D_star_n_test)
}

# debug(C_TMLE_truncation)

# Simulations -------------------------------------------------------------
# Compute true value of EY^d
compute_Psi_d_MC <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, M){
  # Monte-Carlo estimation of the true value of mean of Yd
  if(type == "L0_unif") 
    L0_MC <- runif(M, min= -positivity_parameter, max= positivity_parameter)
  else 
    L0_MC <- rexp(M, rate = 1 / positivity_parameter) * (1 - 2 * rbinom(M, 1, prob = 0.5))
  A0_MC <- d0(L0_MC)
  g0_MC <- expit(alpha0 * L0_MC)
  PL1givenA0L0_MC <- expit(beta0 + beta1 * A0_MC + beta2 * L0_MC)
  mean(PL1givenA0L0_MC)
}

# Compute true value of truncation induced target parameter by numerical integration
compute_Psi_0_delta <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, delta){
  #g0_dw_w is g0(d(w)|w)
  g0_dw_w <- Vectorize(function(w) d0(w) * expit(alpha0 * w) + (1 - d0(w)) * (1 - expit(alpha0 * w)))
  
  # Q0_dw_w is \bar{Q}_0(d(w)| w)
  Q0_dw_w <- function(w) expit(beta0 + beta1 * d0(w) + beta2 * w)
  
  # q_w is q_w(w)
  if(type == "L0_exp"){
    q_w <- function(w) 1 / 2 * 1 / positivity_parameter * exp(-abs(w) / positivity_parameter)
  }else{
    q_w <- Vectorize(function(w) 1 / (2 * positivity_parameter) * (abs(w) <= positivity_parameter))
  }
  
  # Define integrand such that \int integrand(w) dw = Psi_0(\delta) and integrate it
  integrand <- Vectorize(function(w) q_w(w) * g0_dw_w(w) / max(delta, g0_dw_w(w)) * Q0_dw_w(w))
  integrate(integrand, lower = -5 * positivity_parameter, upper = 5 * positivity_parameter)$value
}


# Specify the jobs. A job is the computation of a batch. 
# It is fully characterized by the parameters_tuple_id that the batch corresponds to.
# ns <- c((1:9)*100, c(1:9)*1000, 2*c(1:5)*1e4)
ns <- c(50, 100, 150, 200, 250, (3:10) * 100)
delta0s <- 10^(seq(from = -8, to = log(0.4)/log(10), length = 50))
orders <- 0:4
parameters_grid <- rbind(expand.grid(type = "L0_unif", positivity_parameter = c(2, 4), 
                                     alpha0 = 2, beta0 = -3, beta1 = +1.5, beta2 = 1, n = ns, delta0 = delta0s, order = orders),
                         expand.grid(type = "L0_exp", positivity_parameter = c(2, 4), 
                                     alpha0 = 2, beta0 = -3, beta1 = +1.5, beta2 = 1, n = ns, delta0 = delta0s, order = orders))
# parameters_grid <- expand.grid(type = "L0_unif", positivity_parameter = 4,
#                                alpha0 = 1, beta0 = -3, beta1 = +1.5, beta2 = 1, n = ns, orders_set_id = 1:3)

batch_size <- 1000; nb_batchs <- 10

jobs <- kronecker(1:nrow(parameters_grid), rep(1, nb_batchs))
first_seed_batch <- 1:length(jobs) * batch_size
jobs_permutation <- sample(1:length(jobs))
jobs <- jobs[jobs_permutation]
first_seed_batch <- first_seed_batch[jobs_permutation]


# # Compute target parameter for each parameters tuple id
target_parameters <- vector()
tp_params <- unique(parameters_grid[,c("type", "positivity_parameter", "alpha0", "beta0", "beta1", "beta2")])
for(i in 1:nrow(tp_params)) #target_parameters[i] <- NA
  target_parameters[i] <- compute_Psi_0_delta(type = tp_params[i, "type"],
                                              positivity_parameter = tp_params[i, "positivity_parameter"],
                                              alpha0 = tp_params[i, "alpha0"], 
                                              beta0 = tp_params[i, "beta0"],
                                              beta1 = tp_params[i, "beta1"],
                                              beta2 = tp_params[i, "beta2"],
                                              d0 = alwaysTreated0,
                                              delta = tp_params[i, "delta0"])
target_parameters <- cbind(tp_params, target_parameters)

# Save the parameters' grid
write.table(parameters_grid, file = "parameters_grid_extrapolations.csv", append = F, row.names=F, col.names=T,  sep=",")

# Perform the jobs in parallel
# library(foreach); library(doParallel)
# cl <- makeCluster(getOption("cl.cores", 2), outfile = "")
# registerDoParallel(cl)
library(Rmpi); library(doMPI)
cl <- startMPIcluster(72)
registerDoMPI(cl)
# 
results <- foreach(i = 1:length(jobs)) %dopar% { #job is a parameter_tuple_idS
  # for(i in 1:length(jobs)){
  #   job <- 1
  job <- jobs[i]
  results_batch <- matrix(0, nrow = batch_size, ncol = 4)
  colnames(results_batch) <- c("parameters_tuple_id", "EYd", "TMLE_extrapolation", "TMLE_extrapolation_bis")
  system.time(for(j in 1:batch_size){
    seed <- first_seed_batch[i] + j - 1; #set.seed(seed)
    observed_data <- generate_data(type = parameters_grid[job,]$type, 
                                   positivity_parameter = parameters_grid[job,]$positivity_parameter, 
                                   alpha0 = parameters_grid[job,]$alpha0,
                                   beta0 = parameters_grid[job,]$beta0, beta1 = parameters_grid[job,]$beta1, 
                                   beta2 = parameters_grid[job,]$beta2, n = parameters_grid[job,]$n)
    
    result_extrapolations <- list(result_TMLE_extrapolation = NA, result_TMLE_extrapolation_bis = NA)
    training_set <- 1:parameters_grid[job,]$n
    test_set <- NULL
    try(result_TMLE_extrapolation <- TMLE_extrapolation(observed_data, training_set, 
                                                        test_set, alwaysTreated0, 
                                                        parameters_grid[job,]$order,
                                                        parameters_grid[job,]$delta0,
                                                        Q_misspecified = F)$Psi_n)
    try(result_TMLE_extrapolation_bis <- TMLE_extrapolation_bis(observed_data,
                                                                alwaysTreated0,
                                                                parameters_grid[job,]$order,
                                                                parameters_grid[job,]$delta0,
                                                                Q_misspecified = F))
    #cat("TMLE_extrapolation=", result_TMLE_extrapolation, "\n")
    #cat("TMLE_extrapolation_bis=", result_TMLE_extrapolation_bis, "\n")
    
    target_parameter <- as.numeric(subset(target_parameters, type == parameters_grid[job,]$type
                                          & positivity_parameter ==  parameters_grid[job,]$positivity_parameter
                                          & alpha0 == parameters_grid[job,]$alpha0
                                          & beta0 == parameters_grid[job,]$beta0
                                          & beta1 == parameters_grid[job,]$beta1
                                          & beta2 == parameters_grid[job,]$beta2, select = target_parameters))
    
    results_batch[j, ] <- c(job, target_parameter, result_TMLE_extrapolation, result_TMLE_extrapolation_bis)
  })
  
  if(!file.exists("TMLE_extrapolations_intermediate_results.csv")){
    write.table(results_batch, file="TMLE_extrapolations_intermediate_results.csv", append=T, row.names=F, col.names=T,  sep=",")
  }else{
    write.table(results_batch, file="TMLE_extrapolations_intermediate_results.csv", append=T, row.names=F, col.names=F,  sep=",")
  }
  results_batch
}

save(results, parameters_grid, file = "TMLE_extrapolations_results")

closeCluster(cl)
mpi.quit()