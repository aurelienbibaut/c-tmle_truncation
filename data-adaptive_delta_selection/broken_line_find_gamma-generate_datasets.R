source('./broken_line_find_gamma-functions.R')

data_generating_distributions.parameters <- list()
# Set of parameters 1
data_generating_distributions.parameters[[1]] <- list(lambda = 2, alpha0 = 2, beta0 = -3, beta1 = 1.5, beta2 = 1,
                                                      beta = 0.25, gamma = 0.125)
# Set of parameters 2
data_generating_distributions.parameters[[2]] <- list(lambda = 2, alpha0 = 4, beta0 = -3, beta1 = 1.5, beta2 = 0.5,
                                                      beta = 2 - 5/4, gamma =  1 - 5/8)
# Set of parameters 3
data_generating_distributions.parameters[[3]] <- list(lambda <= 2, alpha0 = 4, beta2 =-3, beta0 = -1, beta1 = 1,
                                                      beta = 7/8, gamma = 5/16)
ns <- c(10^3.5, 10^4.5, 10^5.5)

# Define jobs
nb_repeats <- 20
parameters_grid <- expand.grid(distribution_id = 1:length(data_generating_distributions.parameters), n = ns)
jobs <- sample(kronecker(1:nrow(parameters_grid), rep(1, nb_repeats)))

broken_line.results <- vector()
for(job in jobs){
  current_broken_line.result <- NULL
  try(current_broken_line.result <- generate_data_and_gamma_broken_line('L0_exp', parameters_grid[job, ]$lambda, 
                                                                        parameters_grid[job, ]$alpha0, 
                                                                        parameters_grid[job, ]$beta0, 
                                                                        parameters_grid[job, ]$beta1, 
                                                                        parameters_grid[job, ]$beta2, 
                                                                        parameters_grid[job, ]$n, 
                                                                        parameters_grid[job, ]$gamma)$broken_line_points_df)
  if(!is.null(current_broken_line.result)){
    cat('For job')
    print(parameters_grid[job, ])
    cat('Results:\n')
    print(current_broken_line.result)
    if(!file.exists("broken_lines.results.csv")){
      write.table(current_broken_line.result, file = "broken_lines.results.csv", append = T, row.names = F, col.names = T,  sep = ",")
    }else{
      write.table(current_broken_line.result, file = "broken_lines.results.csv", append = T, row.names = F, col.names = F,  sep = ",")
    }
    broken_line.results <- rbind(broken_line.results, current_broken_line.result)
  }
}