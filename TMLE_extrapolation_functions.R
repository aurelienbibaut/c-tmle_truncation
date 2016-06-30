# Manually define loss function (I guess at some point glm wasn't working fine,
# but I might not have been using it right)
loss <- function(coefficients, outcome, covariates, boolean_subset, offset=0){
  Q_k <- as.vector(expit(covariates %*% coefficients + offset))
  sum(-boolean_subset * (outcome * log(Q_k) + (1 - outcome) * log(1 - Q_k)))
}

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

# TMLE of truncated target parameter
TMLE_EY1 <- function(observed_data, delta, verbose = F){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  
  # 0. Fit models for g_{n,k=0}
  initial_model_for_A0 <- glm(A0 ~ 1 + L0, family=binomial)
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- as.vector(predict(initial_model_for_A0, type="response"))
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  coeffs_Q1d_bar_0n <- glm(L1 ~ L0, family = binomial, subset = which(A0 == 1))$coefficients
  offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0) %*% coeffs_Q1d_bar_0n)
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Compute clever covariate
  gn0_delta <- g_to_g_delta(delta, gn0)
  H_delta <- (A0 == 1) / gn0_delta
  H_delta_setting_A_to_d <- 1 / gn0_delta
  
  # Fit parametric submodel to training set
  epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = offset_vals_Q1d_bar_0n,
                 subset = which(A0 == 1))$coefficients[1]
  if(verbose) cat("delta = ", delta, ", epsilon = ", epsilon, "\n")
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Return estimator
  mean(gn0 / gn0_delta * Q1d_bar_star_n)
}

# Targeting step for TMLE of EY1
targeting_TMLE_EY1 <- function(data, indices){
  replicate_data <- data[indices, ]
  
  L0 <- replicate_data$L0; A0 <- replicate_data$A0; L1 <- replicate_data$L1
  offset_vals_Q1d_bar_0n <- replicate_data$offset_vals_Q1d_bar_0n
  gn0 <- replicate_data$gn0; gn0_delta <- replicate_data$gn0_delta
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Compute clever covariate
  H_delta <- (A0 == 1) / gn0_delta
  H_delta_setting_A_to_d <- 1 / gn0_delta
  
  # Fit parametric submodel
  epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = offset_vals_Q1d_bar_0n,
                 subset = which(A0 == 1))$coefficients[1]
  
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Return estimator
  mean(gn0 / gn0_delta  * Q1d_bar_star_n)
}

# Targeting step for TMLE of EY1 with speedglm
targeting_TMLE_EY1_speedglm <- function(data, indices){
  replicate_data <- data[indices, ]
  
  L0 <- replicate_data$L0; A0 <- replicate_data$A0; L1 <- replicate_data$L1
  offset_vals_Q1d_bar_0n <- replicate_data$offset
  gn0 <- replicate_data$gn0; gn0_delta <- replicate_data$gn0_delta
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Compute clever covariate
  H_delta <- (A0 == 1) / gn0_delta
  H_delta_setting_A_to_d <- 1 / gn0_delta
  
  # Fit parametric submodel
  epsilon <- speedglm(L1 ~ offset(offset) + H_delta - 1,
                      subset(data, A0 == 1), family = binomial())$coefficients[1]
  # if(verbose) cat("delta = ", delta, ", epsilon = ", epsilon, "\n")
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Influence curve and it's variance
  IC <- (as.numeric(A0 == 1)) / gn0_delta * (L1 - Q1d_bar_star_n) +
    gn0 / gn0_delta * Q1d_bar_star_n
  sigma_n <- var(IC)
  
  # Return estimator
  list(Psi_n = mean(gn0 / gn0_delta  * Q1d_bar_star_n), sigma_n = sigma_n)
}

# TMLE of truncated target parameter, with speedglm
TMLE_EY1_speedglm <- function(observed_data, delta, verbose = F){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  data <- data.frame(L0 = L0, A0 = A0, L1 = L1)
  n <- length(L0)
  
  # 0. Fit models for g_{n,k=0}
  initial_model_for_A0 <- speedglm(A0 ~ 1 + L0, data, family = binomial())
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- expit(cbind(1, L0) %*% initial_model_for_A0$coefficients)
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  coeffs_Q1d_bar_0n <- speedglm(L1 ~ L0, subset(data, A0 == 1), family = binomial())$coefficients
  offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0) %*% coeffs_Q1d_bar_0n)
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Compute clever covariate
  gn0_delta <- g_to_g_delta(delta, gn0)
  H_delta <- (A0 == 1) / gn0_delta
  H_delta_setting_A_to_d <- 1 / gn0_delta
  
  data <- cbind(data, H_delta = H_delta, offset = offset_vals_Q1d_bar_0n)
  
  # Fit parametric submodel to training set
  epsilon <- speedglm(L1 ~ offset(offset) + H_delta - 1,
                      subset(data, A0 == 1), family = binomial())$coefficients[1]
  # if(verbose) cat("delta = ", delta, ", epsilon = ", epsilon, "\n")
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Influence curve and its variance
  IC <- (as.numeric(A0 == 1)) / gn0_delta * (L1 - Q1d_bar_star_n) +
    gn0 / gn0_delta * Q1d_bar_star_n
  sigma_n <- var(IC)
  
  # Quantiles of gns
  gns.quantiles <- quantile(gn0, probs = c(0.1, 0.5))
  
  # Return estimator
  list(Psi_n = mean(gn0 / gn0_delta * Q1d_bar_star_n),
       sigma_n = sigma_n,
       gns.quantile_0.025 = gns.quantiles[1],
       gns.quantile_0.975 = gns.quantiles[2])
}

# TMLE of truncated target parameter with 
# bootstrapping of the targeting step, with speedglm
TMLE_EY1_bootstrap_speedglm <- function(observed_data, delta, nb_boostrap_samples = 1000){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  data <- data.frame(L0 = L0, A0 = A0, L1 = L1)
  n <- length(L0)
  
  # 0. Fit models for g_{n,k=0}
  initial_model_for_A0 <- speedglm(A0 ~ 1 + L0, data, family = binomial())
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- expit(cbind(1, L0) %*% initial_model_for_A0$coefficients)
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  coeffs_Q1d_bar_0n <- speedglm(L1 ~ L0, subset(data, A0 == 1), family = binomial())$coefficients
  offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0) %*% coeffs_Q1d_bar_0n)
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Compute clever covariate
  gn0_delta <- g_to_g_delta(delta, gn0)
  H_delta <- (A0 == 1) / gn0_delta
  H_delta_setting_A_to_d <- 1 / gn0_delta
  
  # Prepare data for targeting step
  data <- cbind(data, H_delta = H_delta, offset = offset_vals_Q1d_bar_0n,
                gn0, gn0_delta)
  # Fit parametric submodel to training set
  epsilon <- speedglm(L1 ~ offset(offset) + H_delta - 1,
                      subset(data, A0 == 1), family = binomial())$coefficients[1]
  # if(verbose) cat("delta = ", delta, ", epsilon = ", epsilon, "\n")
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  
  # Bootstrap targeting step + compute estimator on full data
  full_data_Psi_n <- targeting_TMLE_EY1_speedglm(data, 1:n)
  cat('For n=', n, ' and delta =', delta, ', Psi_n = ', full_data_Psi_n, '\n')
  bootstrapped_Psis_n <- boot(data, targeting_TMLE_EY1_speedglm, 
                              nb_boostrap_samples, sim = "ordinary")$t
  
  # Perform Shapiro Wilk's test for normality of the boostrapped estimates
  shapiro.p_value <- shapiro.test(bootstrapped_Psis_n)$p.value
  cat('For n = ', n, ' and delta = ', delta,
      ', Shapiro Wilk p-value = ', shapiro.p_value, '\n')
  
  list(Psi_n = full_data_Psi_n, shapiro.p_value = shapiro.p_value)
}

# TMLE of truncated target parameter with 
# bootstrapping of the targeting step
TMLE_EY1_bootstrap <- function(observed_data, delta, nb_boostrap_samples = 1000){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  
  # 0. Fit models for g_{n,k=0}
  initial_model_for_A0 <- glm(A0 ~ 1 + L0, family=binomial)
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- as.vector(predict(initial_model_for_A0, type="response"))
  gn0_delta <- g_to_g_delta(delta, gn0)
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  coeffs_Q1d_bar_0n <- glm(L1 ~ L0, family = binomial, subset = which(A0 == 1))$coefficients
  offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0) %*% coeffs_Q1d_bar_0n)
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Prepare data frame for targeting step
  targeting_step.df <- data.frame(L0, A0, L1, offset_vals_Q1d_bar_0n,
                                  gn0, gn0_delta)
  
  # Bootstrap targeting step + compute estimator on full data
  full_data_Psi_n <- targeting_TMLE_EY1(targeting_step.df, 1:n)
  cat('For n=', n, ' and delta =', delta, ', Psi_n = ', full_data_Psi_n, '\n')
  bootstrapped_Psis_n <- boot(targeting_step.df, targeting_TMLE_EY1, 
                              nb_boostrap_samples, sim = "ordinary")$t
  
  # Perform Shapiro Wilk's test for normality of the boostrapped estimates
  shapiro.p_value <- shapiro.test(bootstrapped_Psis_n)$p.value
  cat('For n = ', n, ' and delta = ', delta,
      ', Shapiro Wilk p-value = ', shapiro.p_value, '\n')
  
  list(Psi_n = full_data_Psi_n, shapiro.p_value = shapiro.p_value)
}

# Targeting step for TMLE of the finite difference b0(delta + Delta) - b0(delta)
# where b0 the bias function for the truncation induced parameter wrt EY1
targeting_TMLE_fin_diff_EY1 <- function(data, indices, delta, Delta){
  replicate_data <- data[indices, ]
  
  L0 <- replicate_data$L0; A0 <- replicate_data$A0; L1 <- replicate_data$L1
  offset_vals_Q1d_bar_0n <- replicate_data$offset_vals_Q1d_bar_0n
  gn0 <- replicate_data$gn0; gn0_delta <- replicate_data$gn0_delta
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Compute clever covariate
  H_delta <- (A0 == 1) * ((gn0 < delta + Delta) / (delta + Delta) - (gn0 < delta) / delta)
  H_delta_setting_A_to_d <- (gn0 < delta + Delta) / (delta + Delta) - (gn0 < delta) / delta
  
  # Fit parametric submodel
  epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = offset_vals_Q1d_bar_0n,
                 subset = which(A0 == 1))$coefficients[1]
  if(is.na(epsilon)) epsilon <- 0
  
  Q1d_bar_star_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
  
  # Return estimator
  mean(H_delta_setting_A_to_d * gn0 * Q1d_bar_star_n)
}

# TMLE, with bootstrap of the targeting step, of the 
# finite difference b0(delta + Delta) - b0(delta),
# where b0 the bias function for the truncation induced parameter wrt EY1
TMLE_fin_diff_EY1_bootstrap <- function(observed_data, delta, 
                                        Delta, nb_boostrap_samples = 1000){
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  
  # 0. Fit models for g_{n,k=0}
  initial_model_for_A0 <- glm(A0 ~ 1 + L0, family=binomial)
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- as.vector(predict(initial_model_for_A0, type="response"))
  gn0_delta <- g_to_g_delta(delta, gn0)
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  coeffs_Q1d_bar_0n <- glm(L1 ~ L0, family = binomial, subset = which(A0 == 1))$coefficients
  offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0) %*% coeffs_Q1d_bar_0n)
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Prepare data frame for targeting step
  targeting_step.df <- data.frame(L0, A0, L1, offset_vals_Q1d_bar_0n,
                                  gn0, gn0_delta)
  
  # Bootstrap targeting step + compute estimator on full data
  full_data_Psi_n <- targeting_TMLE_fin_diff_EY1(targeting_step.df, 1:n, delta, Delta)
  cat('For n=', n, ' and delta =', delta, ', Psi_n = ', full_data_Psi_n, '\n')
  bootstrapped_Psis_n <- boot(targeting_step.df, targeting_TMLE_fin_diff_EY1, 
                              nb_boostrap_samples, sim = "ordinary", 
                              delta = delta, Delta = Delta)$t
  
  # Perform Shapiro Wilk's test for normality of the boostrapped estimates
  shapiro.p_value <- 1
  try(shapiro.p_value <- shapiro.test(bootstrapped_Psis_n)$p.value)
  cat('For n = ', n, ' and delta = ', delta,
      ', Shapiro Wilk p-value = ', shapiro.p_value, '\n')
  
  list(Psi_n = full_data_Psi_n, shapiro.p_value = shapiro.p_value)
}

# TMLE of the finite difference Psi(delta + Delta) - Psi(delta), where Psi(0) = EY1
# with subsampling (without replacement) and comparison with full data.
# Return the median finite diff estimator on subsampled replicates. 
TMLE_fin_diff_EY1_subsampling <- function(observed_data, delta, 
                                          Delta, nb_subsamples = 1000){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  
  # 0. Fit models for g_{n,k=0}
  initial_model_for_A0 <- glm(A0 ~ 1 + L0, family=binomial)
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- as.vector(predict(initial_model_for_A0, type="response"))
  gn0_delta <- g_to_g_delta(delta, gn0)
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  coeffs_Q1d_bar_0n <- glm(L1 ~ L0, family = binomial, subset = which(A0 == 1))$coefficients
  offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0) %*% coeffs_Q1d_bar_0n)
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Prepare data frame for targeting step
  targeting_step.df <- data.frame(L0, A0, L1, offset_vals_Q1d_bar_0n,
                                  gn0, gn0_delta)
  
  # Bootstrap targeting step + compute estimator on full data
  full_data_Psi_n <- targeting_TMLE_fin_diff_EY1(targeting_step.df, 1:n)
  cat('For n=', n, ' and delta =', delta, ', Psi_n = ', full_data_Psi_n, '\n')
  
  subsampled_Psis_n <- replicate(nb_subsamples,
                                 targeting_TMLE_fin_diff_EY1(targeting_step.df[sample(1:n, floor(1 / sqrt(10)), replace = F), ]))
  
  # Compute median difference between 
  
  # Perform Shapiro Wilk's test for normality of the boostrapped estimates
  shapiro.p_value <- shapiro.test(subsampled_Psis_n)$p.value
  cat('For n = ', n, ' and delta = ', delta,
      ', Shapiro Wilk p-value = ', shapiro.p_value, '\n')
  
  
  list(Psi_n = full_data_Psi_n, shapiro.p_value = shapiro.p_value)
}

# debug(TMLE_EY1)

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