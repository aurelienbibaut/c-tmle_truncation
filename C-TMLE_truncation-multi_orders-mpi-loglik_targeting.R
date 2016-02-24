# Treatment rule(s)
alwaysTreated0 <- function(L0){
  1
}

# Functions to be used later
logit <- function(x){
  log(x/(1-x))
}

expit<-function(x){
  result<-exp(x)/(1+exp(x))
  result[is.nan(result)]<-1
  result
}

g_to_g_delta<-function(delta, g){
  (g<delta)*delta+(g>=delta)*g
}

loss <- function(coefficients, outcome, covariates, boolean_subset, offset=0){
  Q_k <- as.vector(expit(covariates %*% coefficients + offset))
  sum(-boolean_subset * (outcome * log(Q_k) + (1 - outcome) * log(1 - Q_k)))
}

# Data generation
generate_data <- function(R, alpha0, beta0, beta1, beta2, n){
  L0 <- runif(n, min=-R, max=R)
  g00 <- expit(alpha0*L0)
  A0 <- rbinom(n, 1, g00)
  PL1givenA0L0 <- expit(beta0+beta1*A0+beta2*L0)
  L1 <- rbinom(n, 1, PL1givenA0L0)
  list(L0 = L0, A0 = A0, L1 = L1)
}

# Compute a(delta0) (as defined in write up)
compute_a_delta0 <- function(delta0, order, n_points = 9, diff_step = NULL, verbose = F){
  
#   cat("delta0 = ", delta0, " and order = ", order, "\n")
  
  if(order <= 0) return(list(a_delta0 = 1, deltas = delta0))
  
  if(n_points %% 2 == 0) n_points <- n_points + 1
  
  if(is.null(diff_step)){
    if(delta0 - (n_points-1)/2*1e-3 > 0){ diff_step=1e-3 }else{ diff_step = delta0/(n_points-1) }
  }
  bw <- diff_step*2
  deltas <- delta0 + (1:n_points-1-(n_points-1)/2)*diff_step
  weights <- exp(-(deltas-delta0)^2/(2*bw^2))/sqrt(2*pi*bw^2)
  
  X <- outer(deltas - delta0, 0:order, "^")
  
  
  A <- apply(diag(nrow(X)), 2, function(C) lm.wfit(X, C, weights)$coefficients)
  a_delta0 <- (-delta0)^(0:order) %*% A
  list(a_delta0 = a_delta0, deltas = deltas)
}

# debug(compute_a_delta0)

# C-TMLE_truncation
C_TMLE_truncation <- function(observed_data, d0, orders, delta0s, Q_misspecified = F, 
                              n_points = 11, diff_step = NULL){
  
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
                                                  boolean_subset=(A0==d0(L0)))$par
    offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0, d0(L0)) %*% coeffs_Q1d_bar_0n)
  }else{
    offset_vals_Q1d_bar_0n <- rep(logit(mean(L1[A0==d0(L0)])), n)
  }
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Cross validated targeting: compute the loss for \bar{Q} for 
  # every (order, delta0) pair.
  test_set_indices <- sample(1:n, floor(n/5), replace = F)
  training_set_indices <- setdiff(1:n, test_set_indices)
  best_loss <- mean(((L1 - Q1d_bar_0n)^2 * (A0 == d0(L0)))[test_set_indices])
  best_Q1d_bar_n <- Q1d_bar_0n
  best_indices <- list(order = -1, delta0 = 0) #(-1,0) stands for no targeting
  
  for(order in orders){
    for(delta0 in delta0s){
#       cat("Starting targeting step for (", order, ", ", delta0, ")\n")
      
      # Define test set and training set
      test_set_indices <- sample(1:n, floor(n/5), replace = F)
      training_set_indices <- setdiff(1:n, test_set_indices)
      
      # Compute a_delta0 as defined in write-up
      result_compute_a_delta0 <- compute_a_delta0(delta0, order = order, n_points, diff_step, verbose=F)
      a_delta0 <- result_compute_a_delta0$a_delta0
      deltas <- result_compute_a_delta0$deltas
      
      # Compute clever covariate
      gn0_delta0 <- g_to_g_delta(delta0, gn0)
      gn0_deltas <- sapply(deltas, g_to_g_delta, g=gn0)
      H_delta <- (A0==d0(L0)) / gn0 * (outer(gn0, rep(1, length(deltas))) / gn0_deltas) %*% t(a_delta0)
      H_delta_setting_A_to_d <- 1 / gn0 * (outer(gn0, rep(1, length(deltas))) / gn0_deltas) %*% t(a_delta0)
      
      # Fit parametric submodel to training set
      epsilon <- glm(L1 ~ H_delta - 1, family = binomial, offset = Q1d_bar_0n,
                     subset = intersect(which(A0 == d0(L0)), training_set_indices))$coefficients[1]
      candidate_Q1d_bar_n <- expit(logit(Q1d_bar_0n) + epsilon * H_delta_setting_A_to_d)
      
      # Compute loss of the candidate on training set. Update best loss and best candidate
      # if needed
      loss_evaluation_support <- intersect(which(A0 == d0(L0)), test_set_indices)
      
      candidate_loss <- - sum(log(candidate_Q1d_bar_n)[intersect(loss_evaluation_support,
                                                                which(L1 == 1))]) -
                          sum(log(1 - candidate_Q1d_bar_n)[intersect(loss_evaluation_support,
                                               which(L1 == 0))])
      
      if(is.nan(candidate_loss)){
        cat("order = ", order, "and delta0 = ", delta0, "\n")
        cat("candidate_loss = ", candidate_loss, " and best_loss = ", best_loss, "\n")
        #browser()
      }
      if(candidate_loss < best_loss){
        best_loss <- candidate_loss
        best_Q1d_bar_n <- candidate_Q1d_bar_n
        best_indices <- list(order = order, delta0 = delta0)
      }
    }
  }
  # End of cross validated targeting
  
  # Compute estimate. 
  # First recompute deltas, a_delta0, deltas and gn0_deltas for the pair of best indices
  result_compute_a_delta0 <- compute_a_delta0(best_indices$delta0, order = best_indices$order, 
                                              n_points, diff_step, verbose = F)
  a_delta0 <- result_compute_a_delta0$a_delta0
  deltas <- result_compute_a_delta0$deltas
  gn0_deltas <- sapply(deltas, g_to_g_delta, g=gn0)
  # Then compute the TMLE Psi_n(order, delta0)
  Psi_n <- mean(best_Q1d_bar_n * (outer(gn0, rep(1, length(deltas))) / gn0_deltas) %*% t(a_delta0))
  
  # For the record, compute the untargeted G-computation formula, and return outputs
  Utgtd_Psi_n <- mean(Q1d_bar_0n)
  list(Utgtd_Psi_n = Utgtd_Psi_n, Psi_n = Psi_n, tp_indices = best_indices)
}


# Simulations -------------------------------------------------------------
set.seed(0)
delta0s <- c(1e-4, 5e-4, (1:9)*1e-3, (1:9)*1e-2, (1:4)*1e-1)
# delta0s <- 1e-4
orders <- 0:10
# orders <- 9

# Compute true value of EY^d
compute_Psi_d_MC <- function(R, alpha0, beta0, beta1, beta2, d0, M){
  # Monte-Carlo estimation of the true value of mean of Yd
  L0_MC <- runif(M, min=-R, max=R)
  A0_MC <- d0(L0_MC)
  g0_MC <- expit(alpha0 * L0_MC)
  PL1givenA0L0_MC <- expit(beta0 + beta1 * A0_MC + beta2 * L0_MC)
  mean(PL1givenA0L0_MC)
}

Psi_d0 <- compute_Psi_d_MC(R = 4, alpha0 = 2, beta0 = -3, beta1 = -1.5, beta2 = -2, alwaysTreated0, M = 1e6)

# Specify the jobs. A job is the computation of a batch. 
# It is fully characterized by the parameters_tuple_id that the batch corresponds to.
ns <- c((1:9)*100, c(1:9)*1000, 2*c(1:5)*1e4)
parameters_grid <- expand.grid(R = 4, alpha0 = 2, beta0 = -3, beta1 = -1.5, beta2 = -2, n = ns)
batch_size <- 1; nb_batchs <- 200
jobs <- kronecker(1:nrow(parameters_grid), rep(1, nb_batchs))
jobs <- sample(jobs)


# # Perform the jobs in parallel
library(Rmpi); library(doMPI)

cl <- startMPIcluster(32)
registerDoMPI(cl)

results <- foreach(i=1:length(jobs)) %dopar% { #job is a parameter_tuple_id
  i <- 1
  job <- jobs[i]
  results_batch <- matrix(0, nrow = batch_size, ncol = 5)
  colnames(results_batch) <- c("parameters_tuple_id", "Utgtd", "C-TMLE", "order", "delta0")
  for(i in 1:batch_size){
    observed_data <- generate_data(R = parameters_grid[job,]$R, alpha0 = parameters_grid[job,]$alpha0, 
                                   beta0 = parameters_grid[job,]$beta0, beta1 = parameters_grid[job,]$beta1, 
                                   beta2 = parameters_grid[job,]$beta2, n = parameters_grid[job,]$n)
    result_C_TMLE <- list(Utgtd_Psi_n = NA, Psi_n = NA)
    try(result_C_TMLE <- C_TMLE_truncation(observed_data, alwaysTreated0, orders, 
                                           delta0s, Q_misspecified = F))
    print(result_C_TMLE)
    results_batch[i, ] <- c(job, result_C_TMLE$Utgtd_Psi_n, result_C_TMLE$Psi_n, 
                            result_C_TMLE$tp_indices$order, result_C_TMLE$tp_indices$delta0)
  }

  if(!file.exists("C-TMLE_multi_orders_intermediate_results.csv")){
    write.table(results_batch, file="C-TMLE_multi_orders_intermediate_results.csv", append=T, row.names=F, col.names=T,  sep=",")
  }else{
    write.table(results_batch, file="C-TMLE_multi_orders_intermediate_results.csv", append=T, row.names=F, col.names=F,  sep=",")
  }
  results_batch
}

save(results, parameters_grid, file = "C-TMLE_multi_orders_results")

closeCluster(cl)
mpi.quit()


# # Combine the results together
# full_results_matrix <- results[[1]]
# for(i in 2:length(results)){
#   full_results_matrix <- rbind(full_results_matrix, results[[i]])
# }
# 
# # Compute the MSEs for each parameter tuple id
# MSEs <- matrix(NA, nrow = nrow(parameters_grid), ncol = 4)
# for(i in 1:nrow(parameters_grid)){
#   MSE_utgtd <- mean((full_results_matrix[full_results_matrix[,"parameters_tuple_id"] == i, "Utgtd"] - Psi_d0)^2)
#   MSE_C_TMLE <- mean((full_results_matrix[full_results_matrix[,"parameters_tuple_id"] == i, "C-TMLE"] - Psi_d0)^2)
#   MSEs[i, ] <- c(i, parameters_grid[i,]$n, MSE_utgtd, MSE_C_TMLE)
# }
# colnames(MSEs) <- c("parameter_tuple_id", "n", "MSE_utgtd", "MSE_C-TMLE")
# 
# plot(MSEs[, "n"],  MSEs[, "n"]* MSEs[,"MSE_utgtd"])
# lines(MSEs[, "n"],  MSEs[, "n"]* MSEs[,"MSE_C-TMLE"])
