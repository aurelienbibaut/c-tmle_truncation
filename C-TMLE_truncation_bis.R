source('utilities.R')
source('generate_data.R')

library(speedglm)

deltas <- sort(10^seq(from = -3, to = -1, by = 0.25), decreasing = T)

# Generate data
observed_data <- generate_data('L0_exp', 2, 4, 1, 2, 0.5, 1e3)
d0 <- alwaysTreated0

C_TMLE <- function(observed_data, d0, deltas){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1<- observed_data$L1
  observed_data.df <- data.frame(L0 = L0, A0 = A0, L1 = L1)
  
  # Fit models for g
  gn.fit <- speedglm(A0 ~ 1 + L0, observed_data.df, family = binomial())
  gn.fit$coefficients[is.na(gn.fit$coefficients)] <- 0
  gn <- expit(cbind(1, L0) %*% gn.fit$coefficients) # gn0 = P(A0 = 1 | L0)
  gn_d <- (A0 == 1) * gn + (A0 == 0) * (1 - gn) # gn0_d = P(A0 = d0(L0)| L0)
  
  # Fit initial model for Q
  Q1_n_d.coeffs <- speedglm(L1 ~ L0 + A0, subset(observed_data.df, A0 == 1), 
                          family = binomial())$coefficients
  Q1_n_d.coeffs[is.na(Q1_n_d.coeffs)] <- 0
  Q1_n_d.lin_pred <- as.vector(cbind(1, L0, A0) %*% Q1_n_d.coeffs)

  # First targeting step
  delta <- deltas[1]
  gn_d_delta <- g_to_g_delta(delta, gn_d)
  C_delta_d <- 1 / gn_d_delta # clever covariate, setting A0 to d0(A0)
  observed_data.df <- cbind(observed_data.df, C_delta_d = C_delta_d,
                            offset = Q1_n_d.lin_pred)
  epsilon <- speedglm(L1 ~ offset(offset) + C_delta_d - 1,
                      subset(observed_data.df, A0 == 1), 
                      family = binomial())$coefficients[1]
  Q1_n_d_star.lin_pred <- Q1_n_d.lin_pred + epsilon * C_delta_d
  
  # Build a sequence of TMLEs with increasing fit
  Qns_d.lin_pred <- list(Q1_n_d.lin_pred)
  Qn_d_stars.lin_pred <- list(Q1_n_d_star.lin_pred)
  gn_d_deltas <- list(gn_d_delta)
  
  for(k in 2:length(deltas)){
    # Perform targeting step based on gn_d_delta, fluctuating Qn_{k-1}
    delta <- deltas[k]
    gn_d_delta <- g_to_g_delta(delta, gn_d_delta)
    observed_data.df$offset <- Qns_d.lin_pred[[length(Qns_d.lin_pred)]]
    
    C_delta_d <- 1 / gn_d_delta # clever covariate, setting A0 to d0(A0)
    epsilon <- speedglm(L1 ~ offset(offset) + C_delta_d - 1,
                        subset(observed_data.df, A0 == 1), 
                        family = binomial())$coefficients[1]
    candidate_Q_n_d_star.lin_pred <- Qns_d.lin_pred[[length(Qns_d.lin_pred)]] + epsilon * C_delta_d
    
    # Does the candidate Q_n_star have a better goodness of fit than the current latest Q_n_star
    candidate_Q_n_d_star <- expit(candidate_Q_n_d_star.lin_pred)
    candidate_Q_n_d_star.loss <- sum((A0 == d0(L0)) * (L1 * log(candidate_Q_n_d_star) + (1 - L1) * log(1 - candidate_Q_n_d_star)))
    current_Q_n_d_star <- expit(Qns_d.lin_pred[[length(Qns_d.lin_pred)]])
    current_Q_n_d_star.loss <- sum((A0 == d0(L0)) * (L1 * log(current_Q_n_d_star) + (1 - L1) * log(1 - current_Q_n_d_star)))
    
    if(candidate_Q_n_d_star.loss < current_Q_n_d_star.loss){
      cat("k = ", k, " fit improves directly\n")
      Qns_d.lin_pred[[k]] <- Qns_d.lin_pred[[k - 1]]
      Qn_d_stars.lin_pred[[k]] <- candidate_Q_n_d_star.lin_pred
      gn_d_deltas[[k]] <- gn_d_delta
    }else{
      cat("k =", k, " fluctuating Qn_d_star_{k-1}\n")
      # Fluctuate Q_n_star_{k-1} based on g_n_delta
      observed_data.df$offset <- Qn_d_stars.lin_pred[[k - 1]]
      epsilon <- speedglm(L1 ~ offset(offset) + C_delta_d - 1,
                          subset(observed_data.df, A0 == 1),
                          family = binomial())$coefficients[1]
      Qn_d_stars.lin_pred[[k]] <- Qn_d_stars.lin_pred[[k - 1]] + epsilon * C_delta_d
    }
  }
  
  
  # Use cross validation to pick the best fit
  
}