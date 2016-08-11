source('./find_beta-functions-no_bootstrap.R')
library(R.utils)

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


# Parameters
generate_beta_datapoint <- function(plotting = F){
  Delta.delta_rates <- c(0.8, 1, 1.1, 1.375, 1.5)
  
  # Sample a data generating distribution
  current_data_generating_distributions.parameters <- sample_datagen_dist.parameters(runif(1, min = 2, max = 10))
  beta <- (current_data_generating_distributions.parameters$alpha0 - 1 / current_data_generating_distributions.parameters$lambda -
             max(0, current_data_generating_distributions.parameters$beta2)) / current_data_generating_distributions.parameters$alpha0
  n <- floor(10^runif(1, min = 3, max = 4.8))
  
  fin_diffs.results <- compute_finite_difference("L0_exp", current_data_generating_distributions.parameters$lambda, 
                                                 current_data_generating_distributions.parameters$alpha0, 
                                                 current_data_generating_distributions.parameters$beta0, 
                                                 current_data_generating_distributions.parameters$beta1, 
                                                 current_data_generating_distributions.parameters$beta2, 
                                                 n, Delta.delta_rates)
  fin_diffs.results <- subset(fin_diffs.results, fin_diff != 0)
  
  if(plotting){
    fin_diffs.plot <- plot_finite_differences(fin_diffs.results, beta)
    print(fin_diffs.plot)
  }
  
  LTS_fits <- vector()
  for(Delta.delta_rate in Delta.delta_rates){
    LTS_fits <- rbind(LTS_fits,
                      c(lts_regression(fin_diffs.results, Delta.delta_rate), Delta.delta_rate = Delta.delta_rate))
  }
  
  if(plotting) fin_diffs.plot <- fin_diffs.plot + geom_abline(data = as.data.frame(LTS_fits),
                                                              aes(intercept = intercept, slope = slope, 
                                                                  colour = factor(Delta.delta_rate)))
  
  broken_lines.results <-  matrix(NA, nrow = 20, ncol = 16)
  colnames(broken_lines.results) <- c("x.start", "x.end", "y.start",
                                      "y.end", "Delta.delta_rate", "RSS", "r_squared", "nb_points", "leftmost",
                                      "segment.squared_length", "slope", "loss", "is_best", "relative_segment_squared_length",
                                      "linearity", "relative_linearity")
  for(nb_breakpoints in 1:5){
    try(broken_lines.results[sum(0:nb_breakpoints):(sum(0:(nb_breakpoints + 1)) - 1), ] <- 
          as.matrix(fit_broken_line(results_df = subset(fin_diffs.results, fin_diff !=0 & Delta.delta_rate == 1.1),
                                    nb_breakpoints = nb_breakpoints, beta, 1.1, T)$broken_line_points_df))
  }
  broken_lines.results <- cbind(id = 1, segment_id = 1:nrow(broken_lines.results), broken_lines.results)
  broken_lines.results.wide <- reshape(as.data.frame(broken_lines.results), v.names = setdiff(colnames(broken_lines.results), 
                                                                                              c("x.start", "x.end", "y.start", "y.end", "id")), 
                                       timevar = "segment_id", idvar = "id", direction = "wide", 
                                       drop = c("x.start", "x.end", "y.start", "y.end"))
  
  
  LTS_fits.wide <- reshape(as.data.frame(cbind(id = 1, LTS_fits)), v.names = c("intercept", "slope", "r.squared"),
                           timevar = "Delta.delta_rate", idvar = "id", direction = "wide")
  
  cbind(broken_lines.results.wide, LTS_fits.wide, true_beta = beta)
  
}

# Set up cluster
cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()), outfile = '')
registerDoParallel(cl)

# debug(generate_beta_datapoint)
# for(i in 1:100){
find_beta.results <- vector()
for(i in 1:100){
  find_beta.results <- 
    find_beta.results <- rbind(find_beta.results, 
                               generate_beta_datapoint(T))
}
stopCluster(cl)
