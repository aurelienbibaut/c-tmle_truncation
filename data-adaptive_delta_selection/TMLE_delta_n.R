running_environment <- 'AWS'
source('./find_gamma-functions.R')
source('./find_beta-functions-no_bootstrap.R')

library(foreach); library(doParallel)
library(h2o)
library(ggplot2); library(gridExtra)

useEnsemble <- T
if(useEnsemble) library(h2oEnsemble)

localH2O <- h2o.init(nthreads = -1)

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

estimate_optimal_rate <- function(observed_data){
  # Extract features
  # Compute finite differences
  Delta.delta_rates <- c(0.8, 1, 1.1, 1.375, 1.5)
  n <- length(observed_data$L0)
  fin_diffs_df <- compute_finite_difference(observed_data, Delta.delta_rates, n)
  
  # Compute the empirical variances of the IC at different values of the truncation level
  var_IC_df <- compute_variances(observed_data)
  var_IC.plot <- NULL
  
  # Extract gamma features
  gamma_features <- NULL
  cat('About to call extract_gamma_features')
  try(gamma_features <- extract_gamma_features(var_IC_df, 0.5, plotting = F, 
                                               var_IC.plot = var_IC.plot))
  
  # Extract bete features
  cat('About to call beta_features')
  beta_features <- NULL
  try(beta_features <- extract_beta_features(fin_diffs_df, 0.5, plotting = F))
  
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
  features.h2o <- as.h2o(as.data.frame(features))
  if(useEnsemble){
    h2o.rate_regression_fit <- h2o.load_ensemble("/home/rstudio/c-tmle_truncation/data-adaptive_delta_selection/rate_regression.ensemble_fit/", 
                                                 import_levelone = F)
    predicted_delta_rate <- as.vector(predict(h2o.rate_regression_fit, features.h2o)$pred)
  }else{
    h2o.rate_regression_fit <- h2o.loadModel("/home/rstudio/c-tmle_truncation/data-adaptive_delta_selection/rate_regression.deeplearning_fit/DeepLearning_model_R_1472825287709_4")
    predicted_delta_rate <- as.vector(h2o.predict(h2o.rate_regression_fit, features.h2o))
  }
  
  
  
  # Compute delta_n based on predicted rate
  if(predicted_delta_rate < 0) stop('Negative predicted rate. This does not make sense')
  
  predicted_delta_rate
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
  #   current_data_generating_distributions.parameters <- sample_datagen_dist.parameters(runif(1, min = 2, max = 10))
  #   beta <- (current_data_generating_distributions.parameters$alpha0 - 1 / current_data_generating_distributions.parameters$lambda -
  #              max(0, current_data_generating_distributions.parameters$beta2)) / current_data_generating_distributions.parameters$alpha0
  #   gamma <- current_data_generating_distributions.parameters$gamma
  # Set of parameters 3
  current_data_generating_distributions.parameters <- list(lambda = 2, alpha0 = 4, beta2 =-3, beta0 = -1, beta1 = 1,
                                                           beta = 7/8, gamma = 1/16)
  
  # Compute true quantities
  lambda <- current_data_generating_distributions.parameters$lambda
  alpha0 <- current_data_generating_distributions.parameters$alpha0
  beta0 <- current_data_generating_distributions.parameters$beta0
  beta1 <- current_data_generating_distributions.parameters$beta1
  beta2 <- current_data_generating_distributions.parameters$beta2
  # Compute true optimal constant
  bias.kappa <- (alpha0 + max(beta2, 0) + lambda^-1) / alpha0
  var.kappa <- (alpha0 + abs(beta2) + lambda^-1) / alpha0
  beta <- 2 - bias.kappa
  gamma <- 1 - var.kappa / 2
  
  C_b <- lambda^- 1 /  2 * exp((beta0 + beta1) * (beta2 > 0)) * 1 / (lambda^-1 + max(beta2, 0) + alpha0)
  C_sigma2 <- lambda^-1 / 2 * exp(sign(beta2) * (beta0 + beta1)) / (lambda^-1 + alpha0 + abs(beta2))
  
  abs_optimal_rate <- 1 / (2 * (gamma + 1 - beta))
  optimal_const <- (C_sigma2 / C_b^2 * gamma / (1 - beta))^abs_optimal_rate
  
  # n <- floor(10^runif(1, min = 3, max = 4.8))
  n <- ns[which(rmultinom(1, 1, (1:length(ns)) / length(ns)) == 1)]
  
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
  
  
  baseline_delta_n.01 <- 0.01 * (n / 100)^(-1/3)
  optimal_delta_n.01 <- 0.01 * (n / 100)^(-abs_optimal_rate)
  baseline_delta_n.05 <- 0.05 * (n / 100)^(-1/3)
  optimal_delta_n.05 <- 0.05 * (n / 100)^(-abs_optimal_rate)
  optimal_delta_n <- optimal_const * n^(-abs_optimal_rate)
  
  
  baseline_Psi_n.0.01 <- TMLE_EY1_speedglm(observed_data, baseline_delta_n.01)
  optimal_Psi_n.0.01 <- TMLE_EY1_speedglm(observed_data, optimal_delta_n.01)
  optimal_Psi_n <- TMLE_EY1_speedglm(observed_data, optimal_delta_n)
  baseline_Psi_n.0.05 <- TMLE_EY1_speedglm(observed_data, baseline_delta_n.05)
  optimal_Psi_n.0.05 <- TMLE_EY1_speedglm(observed_data, optimal_delta_n.05)
  
  # Estimate optimal rate
  try.estimate_rate.result <- try(estimated_optimal_rate <- estimate_optimal_rate(observed_data))
  if(class(try.estimate_rate.result) == "try-error") stop('Could not estimate rate')
  delta_n.0.01 <- 0.01 * (n / 100)^(-estimated_optimal_rate)
  delta_n.0.05 <- 0.05 * (n / 100)^(-estimated_optimal_rate)
  delta_n <- true_optimal_const * n^(-estimated_optimal_rate)
  Psi_n.0.01 <- TMLE_EY1_speedglm(observed_data, delta_n.0.01)
  Psi_n.0.05 <- TMLE_EY1_speedglm(observed_data, delta_n.0.05)
  Psi_n <- TMLE_EY1_speedglm(observed_data, delta_n)
  #   cat('n = ', n, ', EY1 = ', EY1, '\n',
  #       'baseline delta_n 0.01 = ', baseline_delta_n.01, ' TMLE at baseline delta_n 0.01 = ', baseline_Psi_n.0.01$Psi_n, '\n',
  #       'baseline delta_n 0.05 = ', baseline_delta_n.05, ' TMLE at baseline delta_n 0.05= ', baseline_Psi_n.0.05$Psi_n, '\n',
  #       'delta_n.0.01 = ', TMLE_delta_n.result.0.01$delta_n, ', TMLE at delta_n.0.01 = ',  TMLE_delta_n.result.0.01$Psi_n, '\n',
  #       'delta_n.0.05 = ', TMLE_delta_n.result.0.05$delta_n, ', TMLE at delta_n.0.05 = ',  TMLE_delta_n.result.0.05$Psi_n, '\n',
  #       'loss baseline 0.01', (baseline_Psi_n.0.01$Psi_n - EY1)^2,
  #       ', loss baseline 0.05', (baseline_Psi_n.0.05$Psi_n - EY1)^2,
  #       ', loss data adaptive selection', (TMLE_delta_n.result.0.01$Psi_n - EY1)^2, '\n')
  
  list(estimates = rbind(c(EY1 = EY1, n = n, estimator = 'baseline0.01', delta = baseline_delta_n.01, Psi_n = baseline_Psi_n.0.01$Psi_n),
                         c(EY1 = EY1, n = n, estimator = 'baseline0.05', delta = baseline_delta_n.05, Psi_n = baseline_Psi_n.0.05$Psi_n),
                         c(EY1 = EY1, n = n, estimator = 'optimal0.01', delta = optimal_delta_n.01, Psi_n = optimal_Psi_n.0.01$Psi_n),
                         c(EY1 = EY1, n = n, estimator = 'optimal0.05', delta = optimal_delta_n.05, Psi_n = optimal_Psi_n.0.05$Psi_n),
                         c(EY1 = EY1, n = n, estimator = 'optimal', delta = optimal_delta_n, Psi_n = optimal_Psi_n$Psi_n),
                         c(EY1 = EY1, n = n, estimator = 'data_adaptive.0.01', delta = delta_n.0.01, Psi_n = Psi_n.0.01$Psi_n),
                         c(EY1 = EY1, n = n, estimator = 'data_adaptive.0.05', delta = delta_n.0.05, Psi_n = Psi_n.0.05$Psi_n),
                         c(EY1 = EY1, n = n, estimator = 'data_adaptive', delta = delta_n, Psi_n = Psi_n$Psi_n)),
       rate.estimate = estimated_optimal_rate,
       true_rate = optimal_rate, n = n)
}

debug(simulate_and_estimate_once)

# Repeat the dataset simulation and estimation process a bunch of times

plot_MSEs <- function(results.estimates, results.rates){
  results_df <- transform(as.data.frame(results.estimates),
                          EY1 = as.numeric(as.character(EY1)),
                          n = as.numeric(as.character(n)),
                          delta = as.numeric(as.character(delta)),
                          Psi_n = as.numeric(as.character(Psi_n)),
                          estimator =  as.character(estimator))
  # MSEs
  MSEs_df <- vector()
  estimation_tasks <- unique(results_df[, 2:3])
  for(i in 1:nrow(estimation_tasks)){
    estimation_task.df <- subset(results_df, n == estimation_tasks[i, "n"] &
                                   estimator == estimation_tasks[i, "estimator"])
    MSEs_df <- rbind(MSEs_df,
                     c(n = estimation_tasks[i, "n"], estimator = as.character(estimation_tasks[i, "estimator"]),
                       MSE = mean((estimation_task.df$Psi_n - estimation_task.df$EY1)^2)))
  }
  
  MSEs_df <- transform(MSEs_df, n = as.numeric(as.character(n)),
                       MSE = as.numeric(as.character(MSE)))
  
  MSE_plot <- ggplot(MSEs_df, aes(x = n, y = MSE, colour = estimator)) + geom_line() + geom_point()
  rates_plot <- qplot(results.rates[,3], geom="histogram") + geom_vline(xintercept = results.rates[1, 2])
  grid.arrange(MSE_plot, rates_plot, nrow = 2, ncol = 1)
}

results.estimates <- vector()
results.rates <- vector()

# debug(estimate_optimal_rate)

ns <- c(1e3, 5e3, 1e4, 2e4, 5e4)
# Repeat the dataset simulation and estimation process a bunch of times
for(i in 1:1e4){
  cat('Iteration ', i, '\n')
  try.result <- try(iteration.result <- simulate_and_estimate_once())
  if(class(try.result) != "try-error"){
    results.estimates <- rbind(results.estimates,
                               iteration.result$estimates)
    results.rates <- rbind(results.rates,
                           c(n = iteration.result$n, true_rate = iteration.result$true_rate, estimated_rate = iteration.result$rate.estimate))
    if((i %% 5) == 0) plot_MSEs(results.estimates, results.rates)
  }
}

stopCluster(cl)