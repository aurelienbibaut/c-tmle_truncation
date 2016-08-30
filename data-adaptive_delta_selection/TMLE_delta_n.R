running_environment <- 'AWS'
source('./find_gamma-functions.R')
source('./find_beta-functions-no_bootstrap.R')

library(foreach); library(doParallel)

# Load deeplearning fit for rate regression
# load(file = 'rate.deeplearning_fit.RData')

# Sample data-generating distribution's parameters
sample_datagen_dist.parameters <- function(alpha0_max){
  # Define the vertices of a polytope of parameters (alpha0, lambda^-2, abs(beta2))
  # for which the target parameter is weakly identifiable
  vertices <- t(rbind(c(0, 0, 0),
                      c(alpha0_max, 0, 0),
                      c(alpha0_max, alpha0_max, 0),
                      c(alpha0_max, 0, alpha0_max))) # One column per vertex
  # Pick uniformly the parameter vector (alpha0, lambda^-1, abs(beta2))
  # in the polytope defined by the above defined vertices
  unormalized_weights <- runif(4)
  weights <- unormalized_weights / sum(unormalized_weights)
  main_params <- as.vector(vertices %*% weights)
  alpha0 <- main_params[1]; lambda <- 1 / main_params[2]
  beta2 <- (1 - 2 * rbinom(1, 1, 0.5)) * main_params[3]
  
  beta0 <- runif(1, min = -2, max = 2)
  beta1 <- runif(1, min = -2, max = 2)
  gamma <- (alpha0 - 1 / lambda - abs(beta2)) / (2 * alpha0)
  
  list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1,
       beta2 = beta2, gamma = gamma)
}

TMLE_delta_n <- function(observed_data, C0 = 0.05, n0 = 100){
  # Extract features
  # Compute finite differences
  Delta.delta_rates <- c(0.8, 1, 1.1, 1.375, 1.5)
  fin_diffs_df <- compute_finite_difference(observed_data, Delta.delta_rates, n)
  
  # Compute the empirical variances of the IC at different values of the truncation level
  var_IC_df <- compute_variances(observed_data)
  var_IC.plot <- NULL
  
  # Extract gamma features
  gamma_features <- NULL
  cat('About to call extract_gamma_features')
  try(gamma_features <- extract_gamma_features(var_IC_df, gamma, plotting = F, 
                                               var_IC.plot = var_IC.plot))
  
  # Extract bete features
  cat('About to call beta_features')
  beta_features <- NULL
  try(beta_features <- extract_beta_features(fin_diffs_df, beta, plotting = F))
  
  # Generate features vector
  if(!is.null(gamma_features) & !is.null(beta_features)){
    features <- cbind(gamma_features, beta_features)
  }else if(!is.null(beta_features) & is.null(gamma_features)){
    features <- matrix(NA, nrow = 1, ncol = 76 + 278)
    features[, 77:(76 + 278)] <- as.matrix(beta_features)
  }else if(!is.null(gamma_features) & is.null(beta_features)){
    features <- matrix(NA, nrow = 1, ncol = 76 + 278)
    features[, 0:76] <- as.matrix(gamma_features)
  }else{
    stop("Could extract neither beta nor gamma features")
  }
  features <- matrix(as.numeric(as.character(cbind(dataset_id = 1, features))), nrow = 1)
  colnames(features) <- colnames(test_set)[1:355]
  
  # Estimate rate based on deeplearning fit
  h2o.rate_regression_fit <- h2o.loadModel("/home/rstudio/c-tmle_truncation/data-adaptive_delta_selection/rate_regression.deeplearning_fit/DeepLearning_model_R_1472474363291_1")
  features.h2o <- as.h2o(as.data.frame(features))
  predicted_delta_rate <- as.vector(h2o.predict(h2o.rate_regression_fit, features.h2o))
  
  # Compute delta_n based on predicted rate
  if(predicted_delta_rate < 0) stop('Negative predicted rate. This does not make sense')
  delta_n <- C0 * (n / n0)^(-predicted_delta_rate)
  
  # Returns TMLE of target parameter with truncation level delta_n
  list(Psi_n = TMLE_EY1_speedglm(observed_data, delta_n)$Psi_n, delta_n = delta_n)
}

# debug(TMLE_delta_n)
# debug(extract_beta_features)
# Simulation
options(warn = -1)
# Set up cluster
cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()), outfile = '')
registerDoParallel(cl)


simulate_and_estimate_once <-function(){
  # Sample a data generating distribution
  current_data_generating_distributions.parameters <- sample_datagen_dist.parameters(runif(1, min = 2, max = 10))
  beta <- (current_data_generating_distributions.parameters$alpha0 - 1 / current_data_generating_distributions.parameters$lambda -
             max(0, current_data_generating_distributions.parameters$beta2)) / current_data_generating_distributions.parameters$alpha0
  gamma <- current_data_generating_distributions.parameters$gamma
  n <- floor(10^runif(1, min = 3, max = 4.8))
  
  # Compute the true value of EY1
  EY1 <- compute_true_Psi0_delta("L0_exp", current_data_generating_distributions.parameters$lambda, 
                                 current_data_generating_distributions.parameters$alpha0, 
                                 current_data_generating_distributions.parameters$beta0, 
                                 current_data_generating_distributions.parameters$beta1, 
                                 current_data_generating_distributions.parameters$beta2, 
                                 alwaysTreated0, 0)
  
  # Sample a dataset from the above sampled data-generating distribution
  observed_data <- generate_data("L0_exp", current_data_generating_distributions.parameters$lambda, 
                                 current_data_generating_distributions.parameters$alpha0, 
                                 current_data_generating_distributions.parameters$beta0, 
                                 current_data_generating_distributions.parameters$beta1, 
                                 current_data_generating_distributions.parameters$beta2, 
                                 n)
  
  
  baseline_delta_n <- 0.01 * (n / 100)^(-1/3)
  baseline_Psi_n <- TMLE_EY1_speedglm(observed_data, baseline_delta_n)
  TMLE_delta_n.result <- TMLE_delta_n(observed_data)
  cat('EY1 = ', EY1, '\n',
      'baseline delta_n = ', baseline_delta_n, ' TMLE at baseline delta_n = ', baseline_Psi_n$Psi_n, '\n',
      'delta_n = ', TMLE_delta_n.result$delta_n, ', TMLE at delta_n = ',  TMLE_delta_n.result$Psi_n, '\n',
      'loss baseline ', (baseline_Psi_n$Psi_n - EY1)^2, ', loss data adaptive selection', (TMLE_delta_n.result$Psi_n - EY1)^2, '\n')
  c(EY1, baseline_delta_n, baseline_Psi_n$Psi_n, TMLE_delta_n.result$delta_n,  TMLE_delta_n.result$Psi_n)
}

# Repeat the dataset simulation and estimation process a bunch of times
results <- vector()
for(i in 1:100){
  try(results <- rbind(results,
                       simulate_and_estimate_once()))
}

stopCluster(cl)