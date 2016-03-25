# Treatment rule(s)
alwaysTreated0 <- function(L0){
  1
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
  (g<delta) * delta + (g>=delta) * g
}

loss <- function(coefficients, outcome, covariates, boolean_subset, offset=0){
  Q_k <- as.vector(expit(covariates %*% coefficients + offset))
  sum(-boolean_subset * (outcome * log(Q_k) + (1 - outcome) * log(1 - Q_k)))
}

# Data generation
generate_data <- function(R, alpha0, beta0, beta1, beta2, n){
  L0 <- runif(n, min=-R, max=R)
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
                               boolean_subset=(A0==d0(L0)))$par
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
  epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = Q1d_bar_0n,
                 subset = which(A0 == d0(L0)))$coefficients[1]
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Return estimator
  mean(gn0 / gn0_delta * Q1d_bar_star_n)
}

# Untargeted extrapolation
untargeted_extrapolation <- function(observed_data, d0, order, delta0, Q_misspecified = F,
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
TMLE_extrapolation <- function(observed_data, d0, order, delta0, Q_misspecified = F, n_points = 11, diff_step = NULL){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  
  # Fit models for g_{n,k=0}
  initial_model_for_A0 <- glm(A0 ~ 1 + L0, family=binomial)
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
    coeffs_Q1d_bar_0n <- optim(par=c(0,0,0), fn=loss, outcome=L1, 
                               covariates=cbind(1,L0,A0), 
                               boolean_subset=(A0==d0(L0)))$par
    offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0, d0(L0)) %*% coeffs_Q1d_bar_0n)
  }else{
    offset_vals_Q1d_bar_0n <- rep(logit(mean(L1[A0==d0(L0)])), n)
  }
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Fit parametric submodel to training set
  epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = Q1d_bar_0n,
                 subset = which(A0 == d0(L0)))$coefficients[1]
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Return estimator
  mean(Q1d_bar_star_n * (outer(gn0, rep(1, length(deltas))) / gn0_deltas) %*% t(a_delta0))
}


# Simulations -------------------------------------------------------------
set.seed(0)
delta0s <- c(1e-4, 5e-4, (1:9)*1e-3, (1:9)*1e-2, (1:4)*1e-1)
# delta0s <- 1e-4
# orders <- 0:10
orders <- 9

# Compute true target parameter (EYd)
compute_Psi_d_MC <- function(R, alpha0, beta0, beta1, beta2, d0, M){
  # Monte-Carlo estimation of the true value of mean of Yd
  L0_MC <- runif(M, min=-R, max=R)
  A0_MC <- d0(L0_MC)
  g0_MC <- expit(alpha0 * L0_MC)
  PL1givenA0L0_MC <- expit(beta0 + beta1 * A0_MC + beta2 * L0_MC)
  mean(PL1givenA0L0_MC)
}

Psi_d0 <- compute_Psi_d_MC(R = 4, alpha0 = 2, beta0 = -3, beta1 = +1.5, beta2 = 1, alwaysTreated0, M = 1e6)

# Specify the jobs. A job is the computation of a batch. 
# It is fully characterized by the parameters_tuple_id that the batch corresponds to.
# ns <- c((1:9)*100, c(1:9)*1000, 2*c(1:5)*1e4)
# ns <- c(50, 100, 500, 1000, 5000, 1e4)
ns <- c(1e4, 2e4, 3e4)
delta0 <- 1e-1
order <- 3

parameters_grid <- expand.grid(R = 4, alpha0 = 2, beta0 = -3, beta1 = +1.5, beta2 = 1, n = ns)
batch_size <- 10; nb_batchs <- 5
jobs <- kronecker(1:nrow(parameters_grid), rep(1, nb_batchs))
jobs <- sample(jobs)

untargeted_Psi_n <- vector(); targeted_Psi_n <- vector()
for(job in jobs){
  observed_data <- generate_data(R = parameters_grid[job,]$R, alpha0 = parameters_grid[job,]$alpha0, 
                                 beta0 = parameters_grid[job,]$beta0, beta1 = parameters_grid[job,]$beta1, 
                                 beta2 = parameters_grid[job,]$beta2, n = parameters_grid[job,]$n)
  
  untargeted_Psi_n <- c(untargeted_Psi_n, untargeted_extrapolation(observed_data, alwaysTreated0, order, delta0))
  targeted_Psi_n <- c(targeted_Psi_n, TMLE_extrapolation(observed_data, alwaysTreated0, order, delta0, Q_misspecified = F))
}

plot(untargeted_Psi_n, targeted_Psi_n, xlim=c(0,.5))
