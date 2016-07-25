source('./broken_line_find_gamma-functions-mpi.R')
library(R.methodsS3, lib.loc = '~/Rlibs')
library(R.oo, lib.loc = '~/Rlibs')
library(R.utils, lib.loc = '~/Rlibs')

data_generating_distributions.parameters <- list()
# Set of parameters 1
data_generating_distributions.parameters[[1]] <- list(lambda = 2, alpha0 = 2, beta0 = -3, beta1 = 1.5, beta2 = 1,
                                                      beta = 0.25, gamma = 0.125)
# Set of parameters 2
data_generating_distributions.parameters[[2]] <- list(lambda = 2, alpha0 = 4, beta0 = -3, beta1 = 1.5, beta2 = 0.5,
                                                      beta = 2 - 5/4, gamma =  1 - 5/8)
# Set of parameters 3
data_generating_distributions.parameters[[3]] <- list(lambda = 2, alpha0 = 4, beta2 =-3, beta0 = -1, beta1 = 1,
                                                      beta = 7/8, gamma = 5/16)
#ns <- floor(c(10^3.5, 10^4, 10^4.5))
ns <- floor(c(10^3, 10^3.5))

# Define jobs
nb_repeats <- 20
parameters_grid <- expand.grid(distribution_id = 1:length(data_generating_distributions.parameters), n = ns)
jobs <- sample(kronecker(1:nrow(parameters_grid), rep(1, nb_repeats)))
jobs_completed <- vector()

broken_line.results <- vector()
for(job in jobs){
  current_broken_line.result <- NULL
  current_data_generating_distributions.parameters <- data_generating_distributions.parameters[[parameters_grid[job, ]$distribution_id]]
  
  cat('About to perform task ', job, ': n = ', parameters_grid[job, ]$n,
      ' and distribution id ', parameters_grid[job, ]$distribution_id, '\n')
  
  try(current_broken_line.result <- generate_data_and_gamma_broken_line('L0_exp',
                                                                        current_data_generating_distributions.parameters$lambda, 
                                                                        current_data_generating_distributions.parameters$alpha0,
                                                                        current_data_generating_distributions.parameters$beta0, 
                                                                        current_data_generating_distributions.parameters$beta1, 
                                                                        current_data_generating_distributions.parameters$beta2, 
                                                                        parameters_grid[job, ]$n, 
                                                                        current_data_generating_distributions.parameters$gamma,
                                                                        plotting = F))
  try(closeCluster(cl)())
  if(!is.null(current_broken_line.result)){
    cat('For job')
    print(parameters_grid[job, ])
    cat('Results:\n')
    print(current_broken_line.result)
    if(!file.exists("broken_lines.results.csv")){
      current_broken_line.result <- cbind(dataset_id = 1, current_broken_line.result)
      cat("The results file does not exist yet. About to write:\n")
      print(current_broken_line.result)
      write.table(current_broken_line.result, file = "broken_lines.results.csv", append = T, row.names = F, col.names = T,  sep = ",")
    }else{
      n_lines <- countLines("broken_lines.results.csv")[1]
      last_dataset_id <- read.csv("broken_lines.results.csv", skip = n_lines - 2)[1]
      current_broken_line.result <- cbind(dataset_id = as.numeric(last_dataset_id + 1), current_broken_line.result)
      cat("The results file already exists. About to write:\n")
      print(current_broken_line.result)
      write.table(current_broken_line.result, file = "broken_lines.results.csv", append = T, row.names = F, col.names = F,  sep = ",")
    }
    broken_line.results <- rbind(broken_line.results, current_broken_line.result)
    
    jobs_completed <- c(jobs_completed, job)
  }
}
