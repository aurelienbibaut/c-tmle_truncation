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

# C-TMLE_truncation
C_TMLE_truncation <- function(observed_data, d0, deltas){
  
  L0 <- observed_data$L0; A0 <- observed_data$A0; L1 <- observed_data$L1
  n <- length(L0)
  deltas <- sort(deltas, decreasing  = T)
  
  # 0. Fit models for g_{n,k=0}
  initial_model_for_A0 <- glm(A0 ~ 1 + L0, family=binomial)
  initial_model_for_A0$coefficients[is.na(initial_model_for_A0$coefficients)] <- 0
  gn0 <- as.vector(predict(initial_model_for_A0, type="response"))
  
  # 1.a Fit initial model Q^1_{d,n} of Q^1_{0,d}
  coeffs_Q1d_bar_0n <- optim(par=c(0,0,0), fn=loss, outcome=L1, 
                                                covariates=cbind(1,L0,A0), 
                                                boolean_subset=(A0==d0(L0)))$par
  offset_vals_Q1d_bar_0n <- as.vector(cbind(1, L0, d0(L0)) %*% coeffs_Q1d_bar_0n)
  Q1d_bar_0n <- expit(offset_vals_Q1d_bar_0n)
  
  # Targeting steps: iterate targeting until either no likelihood gain is obtained,
  # or we run out of deltas.
  deltas_indices_used <- vector() # Sequence of deltas that will index the TMLES
  current_Q1d_bar_n <- Q1d_bar_0n
  
  test_set_indices <- sample(1:n, floor(n/5), replace = F)
  training_set_indices <- setdiff(1:n, test_set_indices)
  current_loss <- mean(((L1 - current_Q1d_bar_n)^2 * (A0 == d0(L0)))[test_set_indices])
  
  while(length(deltas_indices_used) < length(deltas)){
    
    test_set_indices <- sample(1:n, floor(n/5), replace = F)
    training_set_indices <- setdiff(1:n, test_set_indices)
    remaining_deltas_indices <- setdiff(1:length(deltas), deltas_indices_used)
    
    best_loss <- Inf
    for(i in remaining_deltas_indices){
      gn0_delta <- g_to_g_delta(deltas[i], gn0)
      H_delta <- (A0 == d0(L0)) / gn0_delta
      H_delta_setting_A_to_d <- 1 / gn0_delta
      
      # Fit parametric submodel to training set
      epsilon <- sum((((L1 - current_Q1d_bar_n) * H_delta) * (A0 == d0(L0)))[training_set_indices]) /
                           sum((H_delta^2 * (A0 == d0(L0)))[training_set_indices])
      candidate_Q1d_bar_n <- current_Q1d_bar_n + epsilon * H_delta
      # Compute loss of the candidate on training set
      candidate_loss <- mean(((L1 - candidate_Q1d_bar_n)^2 * (A0 == d0(L0)))[test_set_indices])
      
      if(candidate_loss < best_loss){
        best_loss <- candidate_loss
        index_best_candidate <- i
        best_candidate_Q1d_bar_n <- candidate_Q1d_bar_n
      }
    }
    
    if(best_loss < current_loss){
      deltas_indices_used <- c(deltas_indices_used, index_best_candidate)
      current_Q1d_bar_n <- best_candidate_Q1d_bar_n
      current_loss <- best_loss
    }else{
      break
    }
  }
  # End of iterative targeting
  
  # Compute estimate
  Psi_n <- mean(current_Q1d_bar_n)
  list(Psi_n = Psi_n, delta_sequence = deltas[deltas_indices_used])
}


# Simulations -------------------------------------------------------------
observed_data <- generate_data(R = 4, alpha0 = 2, beta0 = -3, beta1 = -1.5, beta2 = -2, n = 1000)
deltas <- c(1e-4, 5e-4, (1:9)*1e-3, (1:9)*1e-2, (1:4)*1e-1)
print(C_TMLE_truncation(observed_data, alwaysTreated0, deltas))

compute_Psi_d_MC <- function(R, alpha0, beta0, beta1, beta2, d0, M){
  # Monte-Carlo estimation of the true value of mean of Yd
  L0_MC <- runif(M, min=-R, max=R)
  A0_MC <- d0(L0_MC)
  g0_MC <- expit(alpha0*L0_MC)
  PL1givenA0L0_MC <- expit(beta0+beta1*A0_MC+beta2*L0_MC)
  mean(PL1givenA0L0_MC)
}

Psi_d0 <- compute_Psi_d_MC(R = 4, alpha0 = 2, beta0 = -3, beta1 = -1.5, beta2 = -2, alwaysTreated0, M = 1e6)

