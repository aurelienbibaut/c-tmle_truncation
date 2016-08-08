source('./broken_line_find_gamma-functions-mpi.R')
library(R.utils)
library(Rmpi); library(doMPI)

# Set up cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Define jobs
nb_repeats <- 1e4
jobs <- 1:nb_repeats
jobs_completed <- vector()

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

broken_line.results <- vector()
for(job in jobs){
  current_broken_line.result <- NULL
  current_data_generating_distributions.parameters <- sample_datagen_dist.parameters(runif(1, min = 2, max = 10))
  n <- floor(10^runif(1, min = 3, max = 4.8))
  
  cat('About to perform task ', job, ', n = ', n, ', and data-generating distribution\' parameters:\n')
  print(current_data_generating_distributions.parameters)
  
  try(current_broken_line.result <- generate_data_and_gamma_broken_line('L0_exp',
                                                                        current_data_generating_distributions.parameters$lambda, 
                                                                        current_data_generating_distributions.parameters$alpha0,
                                                                        current_data_generating_distributions.parameters$beta0, 
                                                                        current_data_generating_distributions.parameters$beta1, 
                                                                        current_data_generating_distributions.parameters$beta2, 
                                                                        n, 
                                                                        current_data_generating_distributions.parameters$gamma,
                                                                        plotting = F))
  # Figure out to which file to write
  if(!outfile_name.defined){
    file_number <- NULL
    while(file.exists(paste("broken_lines.results", file_number, ".csv", sep = ''))){
      if(is.null(file_number)){
        file_number <- 1
      }else{
        file_number <- file_number + 1
      }
    }
    outfile <- paste("broken_lines.results", file_number, ".csv", sep = '')
    outfile_name.defined <- T
    cat("We'll write results in ", outfile, "\n")
  }
  
  # Write the results to outfile
  if(!is.null(current_broken_line.result)){
    cat('Results:\n')
    print(current_broken_line.result)
    if(!file.exists("broken_lines.results.csv")){
      current_broken_line.result <- cbind(dataset_id = 1, current_broken_line.result)
      cat("The results file does not exist yet. About to write:\n")
      print(current_broken_line.result)
      write.table(current_broken_line.result, file = outfile, append = T, row.names = F, col.names = T,  sep = ",")
    }else{
      n_lines <- countLines("broken_lines.results.csv")[1]
      last_dataset_id <- read.csv("broken_lines.results.csv", skip = n_lines - 2)[1]
      current_broken_line.result <- cbind(dataset_id = as.numeric(last_dataset_id + 1), current_broken_line.result)
      cat("The results file already exists. About to write:\n")
      print(current_broken_line.result)
      write.table(current_broken_line.result, file = outfile, append = T, row.names = F, col.names = F,  sep = ",")
    }
    broken_line.results <- rbind(broken_line.results, current_broken_line.result)
    
    jobs_completed <- c(jobs_completed, job)
  }
}

closeCluster(cl)