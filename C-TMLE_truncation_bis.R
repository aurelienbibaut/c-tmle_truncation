source('utilities.R')
source('generate_data.R')

library(speedglm)


# 
fit_propensity_score <- function(observed_data.df, d0, deg){
  attach(observed_data.df)
  # Fit models for g
  gn.fit <- speedglm(A0 ~ 1 + poly(L0, degree = deg, raw = TRUE), observed_data.df, family = binomial())
  gn <- predict(gn.fit, newdata = observed_data.df, type = 'response')
  (A0 == 1) * gn + (A0 == 0) * (1 - gn) # gn0_d = P(A0 = d0(L0)| L0)
}

# 
fit_Q1_n_d <- function(observed_data.df, d0, intercept_only = F){
  attach(observed_data.df)
  # Fit initial model for Q
  if(intercept_only){
    Q1_n_d.coeffs <- speedglm(L1 ~ A0, subset(observed_data.df, A0 == 1), 
                              family = binomial())$coefficients
    Q1_n_d.coeffs[is.na(Q1_n_d.coeffs)] <- 0
    return(expit(as.vector(cbind(1, A0) %*% Q1_n_d.coeffs)))
  }else{
    Q1_n_d.coeffs <- speedglm(L1 ~ L0 + A0, subset(observed_data.df, A0 == 1), 
                              family = binomial())$coefficients
    Q1_n_d.coeffs[is.na(Q1_n_d.coeffs)] <- 0
    return(expit(as.vector(cbind(1, L0, A0) %*% Q1_n_d.coeffs)))
  }
  
}

# Build a sequence of TMLEs with increasing fits
C_TMLE.build_TMLEs_sequence <- function(observed_data.df, d0, gns_d, Q1_n_d){
  attach(observed_data.df)
  
  fluctuation_models <- list()
  # First targeting step
  C_d <- 1 / gns_d[, 1] # clever covariate, setting A0 to d0(A0)
  Q1_n_d.logit <- logit(Q1_n_d)
  observed_data.df <- cbind(observed_data.df, C_d = C_d,
                            offset = Q1_n_d.logit)
  epsilon <- speedglm(L1 ~ offset(offset) + C_d - 1,
                      subset(observed_data.df, A0 == 1), 
                      family = binomial())$coefficients[1]
  Q1_n_d_star.logit <- Q1_n_d.logit + epsilon * C_d
  
  # Subsequent targeting steps
  Qns_d.logit <- list(Q1_n_d.logit)
  Qn_d_stars.logit <- list(Q1_n_d_star.logit)
  fluctuation_models[[1]] <- as.matrix(t(c(clever_covariate.id = 1, 
                                           epsilon = epsilon)),
                                       nrow = 1)
  
  for(k in 2:ncol(gns_d)){
    # Perform targeting step based on gn_d_delta, fluctuating Qn_{k-1}
    observed_data.df$offset <- Qns_d.logit[[length(Qns_d.logit)]]
    
    C_d <- 1 / gns_d[, k] # clever covariate, setting A0 to d0(A0)
    epsilon <- speedglm(L1 ~ offset(offset) + C_d - 1,
                        subset(observed_data.df, A0 == 1), 
                        family = binomial())$coefficients[1]
    candidate_Q_n_d_star.logit <- Qns_d.logit[[length(Qns_d.logit)]] + epsilon * C_d
    
    # Does the candidate Q_n_star have a better goodness of fit than the current latest Q_n_star?
    candidate_Q_n_d_star <- expit(candidate_Q_n_d_star.logit)
    candidate_Q_n_d_star.loss <- 
      sum((A0 == d0(L0)) * (L1 * log(candidate_Q_n_d_star) + 
                              (1 - L1) * log(1 - candidate_Q_n_d_star)))
    current_Q_n_d_star <- expit(Qns_d.logit[[length(Qns_d.logit)]])
    current_Q_n_d_star.loss <- 
      sum((A0 == d0(L0)) * (L1 * log(current_Q_n_d_star) + 
                              (1 - L1) * log(1 - current_Q_n_d_star)))
    
    if(candidate_Q_n_d_star.loss < current_Q_n_d_star.loss){
      cat("k = ", k, " fit improves directly\n")
      Qns_d.logit[[k]] <- Qns_d.logit[[k - 1]]
      Qn_d_stars.logit[[k]] <- candidate_Q_n_d_star.logit
      # The new fluctuation is obtained from the previous one by replacing 
      # the last clever covariate by the current one
      fluctuation_models[[k]] <- fluctuation_models[[k - 1]]
      fluctuation_models[[k]][nrow(fluctuation_models[[k]]), ] <- c(clever_covariate.id = k, 
                                                                    epsilon = epsilon)
    }else{
      cat("k =", k, " fluctuating Qn_d_star_{k-1}\n")
      # Fluctuate Q_n_star_{k-1} based on g_n_delta
      observed_data.df$offset <- Qn_d_stars.logit[[k - 1]]
      epsilon <- speedglm(L1 ~ offset(offset) + C_d - 1,
                          subset(observed_data.df, A0 == 1),
                          family = binomial())$coefficients[1]
      Qn_d_stars.logit[[k]] <- Qn_d_stars.logit[[k - 1]] + epsilon * C_d
      # The new fluctation model is obtained from the previous one by 
      # adding to it the current clever covariate
      fluctuation_models[[k]] <- rbind(fluctuation_models[[k - 1]],
                                       c(k, epsilon))
    }
  }
  fluctuation_models
}


# Test the shit
# Generate data
observed_data <- generate_data('L0_exp', 2, 4, 1, 2, 0.5, 1e4, alpha1 = -0.01, alpha2 = 1,
                               beta3 = 1)
d0 <- alwaysTreated0
observed_data.df <- data.frame(L0 = observed_data$L0, A0 = observed_data$A0,
                               L1 = observed_data$L1)

gns_d <- sapply(1:3, function(deg) fit_propensity_score(observed_data.df, d0, deg))
gns_d <- cbind(g_to_g_delta(0.2, gns_d[, 1]), g_to_g_delta(0.1, gns_d[, 1]), g_to_g_delta(0.01, gns_d[, 1]), gns_d)
Q1_d_n <- fit_Q1_n_d(observed_data.df, d0, T)
C_TMLE.build_TMLEs_sequence(observed_data.df, d0, gns_d, Q1_d_n)

# Compute the C-TMLE for given initial estimator, sequence of increasingly 
# unbiased propensity scores, and fluctuation model
C_TMLE.predict <- function(Q1_d_n, gns, fluctuation_model){
  
}

# Compute the cross validated log likelihood of a Q_d fit
CV_log_likelihood <- function(Qn_d_star.coeffs, epsilon,
                              observed_data.df, gn_d.coeffs, delta,
                              test_set.indices, d0){
  test_set.df <- observed_data.df[test_set.indices]
  L0 <- test_set.df$L0; A0 <- test_set.df$A0; L1 <-test_set$L1
  Qn_d_star.test_set <- expit(as.vector(cbind(1, L0, A0) %*% Qn_d_star.coeffs))
  sum((A0 == d0(L0)) * (L1 * log(Qn_d_star.test_set) + 
                          (1 - L1) * log(1 - Qn_d_star.test_set)))
}