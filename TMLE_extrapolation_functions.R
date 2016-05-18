# Compute a_delta0 as defined in write-up
compute_a_delta0 <- function(delta0, order){
  if(order <= 0) return(list(a_delta0 = 1, deltas = delta0))
  
  if(delta0 / 2 < 5e-3){
    h <- delta0 / 2
  }else{
    h <- 5e-3
  }
  
  if(order == 1){
    diff1 <- c(-1 / (2 * h), 0, 1 / (2 * h))
    a_delta0 <- t(c(0, 1, 0) - delta0 * diff1)
    return(list(a_delta0 = a_delta0, deltas = c(delta0 - h, delta0, delta0 + h)))
  }
  if(order == 2){
    diff1 <- c(-1 / (2 * h), 0, 1 / (2 * h))
    diff2 <- c(1 / h^2, -2 / h^2, 1 / h^2)
    a_delta_0 <- t(c(0, 1, 0) - delta0 * diff1 + delta0^2 / 2 * diff2)
    return(list(a_delta0 = a_delta_0, deltas = c(delta0 - h, delta0, delta0 + h)))
  }
}

# Compute a(delta0) (as defined in write up)
compute_a_delta0 <- function(delta0, order, n_points = 9, diff_step = NULL, verbose = F){
  
  #   cat("delta0 = ", delta0, " and order = ", order, "\n")
  
  if(order <= 0) return(list(a_delta0 = 1, deltas = delta0))
  
  if(n_points %% 2 == 0) n_points <- n_points + 1
  
  if(is.null(diff_step)){
    if(delta0 - (n_points-1)/2*1e-3 > 0){ diff_step=1e-3 }else{ diff_step = delta0 / (n_points-1) }
  }
  bw <- diff_step * 2
  deltas <- delta0 + (1:n_points-1-(n_points-1) / 2)*diff_step
  weights <- exp(-(deltas-delta0)^2 / (2*bw^2)) / sqrt(2*pi*bw^2)
  
  X <- outer(deltas - delta0, 0:order, "^")
  
  A <- apply(diag(nrow(X)), 2, function(C) lm.wfit(X, C, weights)$coefficients)
  a_delta0 <- (-delta0)^(0:order) %*% A
  list(a_delta0 = a_delta0, differentiator = A, deltas = deltas)
}

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
TMLE_extrapolation <- function(observed_data, d0, order, delta0, Q_misspecified = F,
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