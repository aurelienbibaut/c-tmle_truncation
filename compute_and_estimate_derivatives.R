# Treatment rule(s)
alwaysTreated0 <- function(L0){
  1
}

# Functions to be used later
logit <- function(x){
  log(x/(1-x))
}

expit <- function(x){
  result <- exp(x)/(1+exp(x))
  result[is.nan(result)] <- 1
  result
}

g_to_g_delta<-function(delta, g){
  (g<delta) * delta + (g>=delta) * g
}

# Compute Psi_0(delta) and derivatives and ICs
compute_tps_derivatives_var_ICs <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, delta, lambda = 1){
  #g0_dw_w is g0(d(w)|w)
  g0_dw_w <- Vectorize(function(w) d0(w) * expit(alpha0 * w) + (1 - d0(w)) * (1 - expit(alpha0 * w)))
  
  # Q0_dw_w is \bar{Q}_0(d(w)| w)
  Q0_dw_w <- function(w) expit(beta0 + beta1 * d0(w) + beta2 * w)
  
  # q_w is q_w(w)
  if(type == "L0_exp"){
    q_w <- function(w) 1 / 2 * 1 / positivity_parameter * exp(-abs(w) / positivity_parameter)
  }else{
    q_w <- Vectorize(function(w) 1 / (2 * positivity_parameter) * (abs(w) <= positivity_parameter))
  }
  
  # Define integrand such that \int integrand(w) dw = Psi_0(\delta) and integrate it
  integrand <- Vectorize(function(w) q_w(w) * g0_dw_w(w) / max(delta, g0_dw_w(w)) * Q0_dw_w(w))
  Psi0_delta <- integrate(integrand, lower = -5 * positivity_parameter, upper = 5 * positivity_parameter)$value
  
  # First derivative and second. Might need to fix it for negative values of alpha0
  integrand <- Vectorize(function(w) q_w(w) * g0_dw_w(w) * Q0_dw_w(w))
  if(delta <= 0 || logit(delta) == -Inf){
    first_derivative <- 0 
  }else{
    first_derivative <- -1 / delta^2 * integrate(integrand, lower = min(-10, -3 * logit(delta)), 
                                                 upper = logit(delta) / alpha0)$value
    
    second_derivative <- 2 / delta^3 * integrate(integrand, lower = min(-10, -3 * logit(delta)), 
                                                 upper = logit(delta) / alpha0)$value -
      1 / (delta^3 * (1 - delta)) * 1 / alpha0 * q_w(logit(delta) / alpha0) * 
      g0_dw_w(logit(delta) / alpha0) * Q0_dw_w(logit(delta) / alpha0)
  }
  

  
  # Variance of the influence curve of the plain trucation induced target parameter Psi0(delta)
  integrand1 <- Vectorize(function(w) q_w(w) * g0_dw_w(w) / max(delta, g0_dw_w(w))^2 * (Q0_dw_w(w) - Q0_dw_w(w)^2))
  integrand2 <- Vectorize(function(w) q_w(w) * g0_dw_w(w)^2 / max(delta, g0_dw_w(w))^2 * Q0_dw_w(w)^2)
  var_IC0 <- integrate(integrand1, lower = -5 * positivity_parameter, upper = 5 * positivity_parameter)$value +
    integrate(integrand2, lower = -5 * positivity_parameter, upper = 5 * positivity_parameter)$value - Psi0_delta^2
  
  # Variance of the influence curve of the first order Taylor expansion
  if(delta > 0){
    integrand1 <-  Vectorize(function(w) q_w(w) * g0_dw_w(w) * (1 / max(delta, g0_dw_w(w)) + 
                                                                  lambda * as.numeric(g0_dw_w(w) < delta) / delta)^2 *
                               (Q0_dw_w(w) - Q0_dw_w(w)^2))
    integrand2 <- Vectorize(function(w) q_w(w) * g0_dw_w(w)^2 * (1 / max(delta, g0_dw_w(w)) + 
                                                                   lambda * as.numeric(g0_dw_w(w) < delta) / delta)^2 *
                              Q0_dw_w(w)^2)
  }
  Psi1_delta <- Psi0_delta - first_derivative * delta
  var_IC1 <- integrate(integrand1, lower = -5 * positivity_parameter, upper = 5 * positivity_parameter)$value +
    integrate(integrand2, lower = -5 * positivity_parameter, upper = 5 * positivity_parameter)$value - Psi1_delta^2
  
  list(Psi0_delta = Psi0_delta, first_derivative = first_derivative, 
       second_derivative = second_derivative,
       var_IC0 = var_IC0, var_IC1 = var_IC1)
}

# Compute a(delta0) (as defined in write up)
compute_a_delta0 <- function(delta0, order, n_points = 9, diff_step = NULL, verbose = F){
  if(order <= 0) return(list(a_delta0 = 1, deltas = delta0))
  
  if(n_points %% 2 == 0) n_points <- n_points + 1
  
  # if(is.null(diff_step)){
  #   if(delta0 - (n_points-1)/2*5e-3 > 0){ diff_step=5e-3 }else{ diff_step = delta0 / (n_points-1) }
  # }
  diff_step <- 2 / (n - 1) * delta0^1.4
  bw <- diff_step * 2
  deltas <- delta0 + (1:n_points - (n_points + 1) / 2)*diff_step
  weights <- exp(-(deltas-delta0)^2 / (2*bw^2)) / sqrt(2*pi*bw^2)
  
  X <- outer(deltas - delta0, 0:order, "^")
  
  A <- apply(diag(nrow(X)), 2, function(C) lm.wfit(X, C, weights)$coefficients)
  a_delta0 <- (-delta0)^(0:order) %*% A
  list(a_delta0 = a_delta0, differentiator = A, deltas = deltas)
}

compute_extrapolation <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, 
                                  order, delta0){
  
  result_compute_a0 <- compute_a_delta0(delta0, order)
  deltas <- result_compute_a0$deltas
  a_delta0 <- result_compute_a0$a_delta0
  
  g0_dw_w <- Vectorize(function(w) d0(w) * expit(alpha0 * w) + (1 - d0(w)) * (1 - expit(alpha0 * w)))
  
  # Q0_dw_w is \bar{Q}_0(d(w)| w)
  Q0_dw_w <- function(w) expit(beta0 + beta1 * d0(w) + beta2 * w)
  
  # q_w is q_w(w)
  if(type == "L0_exp"){
    q_w <- function(w) 1 / 2 * 1 / positivity_parameter * exp(-abs(w) / positivity_parameter)
  }else{
    q_w <- Vectorize(function(w) 1 / (2 * positivity_parameter) * (abs(w) <= positivity_parameter))
  }
  
  # Define integrand such that \int integrand(w) dw = Psi_0(\delta) and integrate it
  Psi0_deltas <- vector()
  for(delta in deltas){
    integrand <- Vectorize(function(w) q_w(w) * g0_dw_w(w) / max(delta, g0_dw_w(w)) * Q0_dw_w(w))
    Psi0_deltas <- c(Psi0_deltas, 
                     integrate(integrand, lower = -5 * positivity_parameter, upper = 5 * positivity_parameter)$value)
  }
  Psi0_deltas %*% t(a_delta0)
}

# debug(compute_tps_derivatives_var_ICs)

Psi0 <- function(delta) compute_tps_derivatives_var_ICs("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta)$Psi0_delta
first_derivative <- function(delta) compute_tps_derivatives_var_ICs("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta)$first_derivative
second_derivative <- function(delta) compute_tps_derivatives_var_ICs("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta)$second_derivative
WLS1 <- function(delta) compute_extrapolation("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, 1, delta)
WLS2 <- function(delta) compute_extrapolation("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, 2, delta)
sigma0 <- function(delta) sqrt(compute_tps_derivatives_var_ICs("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta)$var_IC0)
sigma1 <- function(delta) sqrt(compute_tps_derivatives_var_ICs("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta)$var_IC1)
sigma_half <- function(delta) sqrt(compute_tps_derivatives_var_ICs("L0_exp", 2, 2, -3, 1.5, 1, alwaysTreated0, delta, 0.5)$var_IC1)

library(numDeriv)


# Plot derivatives
deltas <- seq(from = 1e-5, to = 0.4, length = 200)
Psi0_vals <- sapply(deltas, Psi0)
WLS1_vals <- sapply(deltas, WLS1)
WLS2_vals <- sapply(deltas, WLS2)
# first_derivative_vals_num_diff <- sapply(deltas, function(x) grad(Psi0, x))
sigma0_vals <- sapply(deltas, function(x) sigma0(x))
sigma1_vals <- sapply(deltas, function(x) sigma1(x))
sigma_half_vals <- sapply(deltas, function(x) sigma_half(x))

first_derivative_vals <- sapply(deltas, first_derivative)
second_derivative_vals <- sapply(deltas, second_derivative)
# second_derivative_vals <- sapply(deltas, function(x) hessian(Psi0, x))
first_order_extr <- Psi0_vals - deltas * first_derivative_vals
middle_combination <- Psi0_vals - deltas / 2 * first_derivative_vals
second_order_extr <- Psi0_vals - deltas * first_derivative_vals + deltas^2 / 2 * second_derivative_vals

plot(deltas, Psi0_vals, ylim=c(0.28, 0.3))
lines(deltas, first_order_extr)
lines(deltas, middle_combination)
abline(a = Psi0_vals[1], b = 0)

lines(deltas, Psi0_vals + 1.96 * sigma0_vals / sqrt(1e5), col = "green")
lines(deltas, Psi0_vals - 1.96 * sigma0_vals / sqrt(1e5), col = "green")

lines(deltas, first_order_extr + 1.96 * sigma1_vals / sqrt(1e5), col = "blue")
lines(deltas, first_order_extr - 1.96 * sigma1_vals / sqrt(1e5), col = "blue")

lines(deltas, first_order_extr + 1.96 * sigma1_vals / sqrt(1e5), col = "blue")
lines(deltas, first_order_extr - 1.96 * sigma1_vals / sqrt(1e5), col = "blue")

lines(deltas, middle_combination + 1.96 * sigma_half_vals / sqrt(1e5), col = "red")
lines(deltas, middle_combination - 1.96 * sigma_half_vals / sqrt(1e5), col = "red")

# lines(deltas, WLS1_vals)
# lines(deltas, WLS2_vals)

lines(deltas, second_order_extr)
