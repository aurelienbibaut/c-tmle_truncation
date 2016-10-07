source('utilities.R')
source('generate_data.R')

library(speedglm); library(caret)

# Initial model fit for the Q part of the likelihood
fit_Q1_n_d <- function(observed_data.df, d0, intercept_only = F){
  if(intercept_only){
    Q1_n_d.model_fit <- speedglm(L1 ~ A0, subset(observed_data.df, A0 == 1), 
                                 family = binomial())
  }else{
    Q1_n_d.model_fit <- speedglm(L1 ~ A0, subset(observed_data.df, A0 == 1), 
                                 family = binomial())
  }
  Q1_n_d.fitted_values <- predict(Q1_n_d.model_fit, newdata = observed_data.df,
                                  type = 'response')
  list(model_fit = Q1_n_d.model_fit, 
       fitted_values = Q1_n_d.fitted_values)
}

# Build a sequence of TMLEs with increasing fits 
# and return the corresponding fluctuation models
C_TMLE.build_TMLEs_sequence <- function(observed_data.df, d0, gns_d, Q1_n_d, max_its = Inf){

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
  if(max_its == 1) return(fluctuation_models)
  for(k in 2:min(ncol(gns_d), max_its)){
    # Perform targeting step based on gn_d_delta, fluctuating Qn_{k-1}
    observed_data.df$offset <- Qns_d.logit[[length(Qns_d.logit)]]
    
    C_d <- 1 / gns_d[, k] # clever covariate, setting A0 to d0(A0)
    epsilon <- speedglm(L1 ~ offset(offset) + C_d - 1,
                        subset(observed_data.df, A0 == 1), 
                        family = binomial())$coefficients[1]
    candidate_Q_n_d_star.logit <- Qns_d.logit[[length(Qns_d.logit)]] + epsilon * C_d
    
    # Does the candidate Q_n_star have a better goodness of fit than the current latest Q_n_star?
    candidate_Q_n_d_star <- expit(candidate_Q_n_d_star.logit)
    candidate_Q_n_d_star.loss.summand <- 
      - (observed_data.df$A0 == d0(observed_data.df$L0)) * (observed_data.df$L1 * log(candidate_Q_n_d_star) + 
                                (1 - observed_data.df$L1) * log(1 - candidate_Q_n_d_star))
    current_Q_n_d_star <- expit(Qns_d.logit[[length(Qns_d.logit)]])
    current_Q_n_d_star.loss.summand <- 
      - (observed_data.df$A0 == d0(L0)) * (observed_data.df$L1 * log(current_Q_n_d_star) + 
                                (1 - observed_data.df$L1) * log(1 - current_Q_n_d_star))
    # It's okay if Q(Y | A, W) = 1 and Y = 1, or Q(Y | A, W) = 0 and Y = 0. Thus I take the corresponding NaN's out
    candidate_Q_n_d_star.loss.summand[is.nan(candidate_Q_n_d_star.loss.summand)] <- 0
    current_Q_n_d_star.loss.summand[is.nan(current_Q_n_d_star.loss.summand)] <- 0
    candidate_Q_n_d_star.loss <- sum(candidate_Q_n_d_star.loss.summand)
    current_Q_n_d_star.loss<- sum(current_Q_n_d_star.loss.summand)
    
    if(candidate_Q_n_d_star.loss < current_Q_n_d_star.loss){
      cat("k = ", k, ", fluctuating Q_n_d_{k-1} improves the fit directly\n")
      Qns_d.logit[[k]] <- Qns_d.logit[[k - 1]]
      Qn_d_stars.logit[[k]] <- candidate_Q_n_d_star.logit
      # The new fluctuation is obtained from the previous one by replacing 
      # the last clever covariate by the current one
      fluctuation_models[[k]] <- fluctuation_models[[k - 1]]
      fluctuation_models[[k]][nrow(fluctuation_models[[k]]), ] <- c(clever_covariate.id = k, 
                                                                    epsilon = epsilon)
    }else{
      cat("k =", k, ", fluctuating Qn_d_star_{k-1}\n")
      # Fluctuate Q_n_star_{k-1} based on g_n_delta
      observed_data.df$offset <- Qn_d_stars.logit[[k - 1]]
      Qns_d.logit[[k]] <- Qn_d_stars.logit[[k - 1]]
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

# Compute the C-TMLE for given initial estimator, sequence of increasingly 
# unbiased propensity scores, and fluctuation model
C_TMLE.predict <- function(Q1_d_n, gns_d, fluctuation_model){
  C_ds <- 1 / gns_d
  expit(logit(Q1_d_n) + 
          as.matrix(C_ds[, fluctuation_model[, "clever_covariate.id"]]) %*% fluctuation_model[, "epsilon.C_d"])
}

# Compute the cross validated log likelihood of a Q_d fit
test_set.log_likelihood <- function(Qn_d_star.test_set, test_set.df){
  L0 <- test_set.df$L0; A0 <- test_set.df$A0; L1 <-test_set.df$L1
  summand <- -((A0 == d0(L0)) * (L1 * log(Qn_d_star.test_set) + 
                          (1 - L1) * log(1 - Qn_d_star.test_set)))
  # It's okay if Q(Y | A, W) = 1 and Y = 1, or Q(Y | A, W) = 0 and Y = 0. Thus I take the corresponding NaN's out
  summand[is.nan(summand)] <- 0
  sum(summand)
}

# Compute increasingly unbiased sequence of propensity scores
# To be modified for one's particular application
fit_propensity_scores_sequence <- function(training_set.df, d0, test_set.df = NULL, max_degree = 3){
  gns_d.training_set <- vector()
  if(!is.null(test_set.df)){
    gns_d.test_set <- vector()
  }else{
    gns_d.test_set <- NULL
  }
  for(deg in 1:max_degree){
    gn.model_fit <- speedglm(A0 ~ 1 + poly(L0, degree = deg, raw = TRUE), 
                             data = training_set.df, family = binomial())
    gn.training_set <- predict(gn.model_fit, newdata = training_set.df,
                               type = 'response')
    if(!is.null(test_set.df)) gn.test_set <- predict(gn.model_fit, newdata = test_set.df,
                                                   type = 'response')
    gns_d.training_set <- cbind(gns_d.training_set,
                                (training_set.df$A0 == d0(training_set.df$L0)) * gn.training_set +
                                  (training_set.df$A0 != d0(training_set.df$L0)) * (1 - gn.training_set))
    if(!is.null(test_set.df)) gns_d.test_set <- cbind(gns_d.test_set,
                                         (test_set.df$A0 == d0(test.df$L0)) * gn.test_set +
                                           (test_set.df$A0 != d0(test.df$L0)) * (1 - gn.test_set))
  }
  
  gns_d.training_set <- cbind(g_to_g_delta(0.2, gns_d.training_set[, 1]), 
                              g_to_g_delta(0.1, gns_d.training_set[, 1]), 
                              g_to_g_delta(0.01, gns_d.training_set[, 1]), gns_d.training_set)
  if(!is.null(test_set.df)) gns_d.test_set <- cbind(g_to_g_delta(0.2, gns_d.test_set[, 1]), 
                                       g_to_g_delta(0.1, gns_d.test_set[, 1]), 
                                       g_to_g_delta(0.01, gns_d.test_set[, 1]), gns_d.test_set)
  list(gns_d.training_set = gns_d.training_set,
       gns_d.test_set = gns_d.test_set)
}

# Main C-TMLE function
C_TMLE <- function(observed_data.df, d0, V = 5, Q.intercept_only = F){
  # Make V folds of the observed dataset
  n <- nrow(observed_data.df)
  folds <- createFolds(1:n, V)
  test_sets.loss <- matrix(NA, nrow = 6, ncol = V) # loss of each TMLE in the sequence on each test set
  
  # Compute the sequence of TMLEs for each fold
  # and its log likelihood loss on the test set
  for(fold_id in 1:length(folds)){
    cat("Fold ", fold_id, "\n")
    test_set.df <- observed_data.df[folds[[fold_id]], ]
    training_set.df <- observed_data.df[setdiff(1:n, folds[[fold_id]]), ]
    # Fit the sequence of TMLEs on the training set
    Q1_d_n.fit <- fit_Q1_n_d(training_set.df, d0, intercept_only = Q.intercept_only)
    gns_d <- fit_propensity_scores_sequence(training_set.df, d0, test_set.df, 3)
    
    fluctuation_models <- C_TMLE.build_TMLEs_sequence(training_set.df, d0, 
                                                      gns_d$gns_d.training_set, Q1_d_n.fit$fitted_values)
    # Compute the loss of each TMLE on the test set
    Q1_d_n.test_set <- predict(Q1_d_n.fit$model_fit, newdata = test_set.df,
                               type = 'response')
    Q_d_n_stars.test_set <- lapply(fluctuation_models,
                                   function(fluctuation_model) C_TMLE.predict(Q1_d_n.test_set,
                                                                              gns_d$gns_d.test_set,
                                                                              fluctuation_model))
    test_sets.loss[, fold_id] <- unlist(lapply(Q_d_n_stars.test_set, 
                                        function(Q_d_n_star.test_set) test_set.log_likelihood(Q_d_n_star.test_set,
                                                                                              test_set.df)))
  }
  
  # Find  best TMLE in the sequence based on the cross validated losses
  CV_losses <- apply(test_sets.loss, 1, sum)
  best_index <- which.min(CV_losses)
  cat("Best TMLE's index:", best_index, "\n")
  # Return the corresponding TMLE based on the full data
  Q1_d_n.full_data <- fit_Q1_n_d(observed_data.df, d0, intercept_only = Q.intercept_only)$fitted_values
  gns_d.full_data <- fit_propensity_scores_sequence(observed_data.df, d0, test_set.df, 3)
  fluctuation_models.full_data <- C_TMLE.build_TMLEs_sequence(training_set.df, d0, 
                                                              gns_d$gns_d.training_set, 
                                                              Q1_d_n.fit$fitted_values,
                                                              max_its = best_index)
  C_TMLE.predict(Q1_d_n.test_set,
                 gns_d$gns_d.test_set, fluctuation_models.full_data[[best_index]])
}

# Test the C-TMLE procedure
observed_data <- generate_data('L0_exp', 2, 4, 1, 2, 0.5, 1e4, 
                               alpha1 = -0.01, alpha2 = 1,
                               beta3 = 1)
d0 <- alwaysTreated0
observed_data.df <- data.frame(L0 = observed_data$L0, 
                               A0 = observed_data$A0,
                               L1 = observed_data$L1)

C_TMLE.Q_d_n_star <- C_TMLE(observed_data.df, d0, V = 5, Q.intercept_only = F)
