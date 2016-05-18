# Generate densities (and conditional expectations) used to define 
# the target parameters and their ICs
generate_densities <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0){
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
  
  list(g0_dw_w = g0_dw_w, Q0_dw_w = Q0_dw_w, q_w = q_w)
}

# Compute Psi_0(delta), the truncation induced target parameter at truncation level delta
compute_true_Psi0_delta <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, delta){
  # Generate and get densities and conditional expectations functions
  densities_and_cond_expectation <- generate_densities(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0)
  g0_dw_w <- densities_and_cond_expectation$g0_dw_w
  Q0_dw_w <- densities_and_cond_expectation$Q0_dw_w
  q_w <- densities_and_cond_expectation$q_w
  
  # Define integrand such that \int integrand(w) dw = Psi_0(\delta) and integrate it
  integrand <- Vectorize(function(w) q_w(w) * g0_dw_w(w) / max(delta, g0_dw_w(w)) * Q0_dw_w(w))
  integrate(integrand, lower = -5 * positivity_parameter, upper = 5 * positivity_parameter)$value
}

# Compute the true variance of the influence curve of Psi0(delta) and of the first order Taylor expansion
compute_true_var_IC_Psi0_delta <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, delta, lambda = 1){
  # Generate and get densities and conditional expectations functions
  densities_and_cond_expectation <- generate_densities(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0)
  g0_dw_w <- densities_and_cond_expectation$g0_dw_w
  Q0_dw_w <- densities_and_cond_expectation$Q0_dw_w
  q_w <- densities_and_cond_expectation$q_w
  
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
  
  list(var_IC0 = var_IC0, var_IC1 = var_IC1)
}

# Compute the first and second derivatives of Psi_0 at delta
compute_true_1st_and_2nd_derivatives_Psi0_delta <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, delta, lambda = 1){
  # Generate and get densities and conditional expectations functions
  densities_and_cond_expectation <- generate_densities(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0)
  g0_dw_w <- densities_and_cond_expectation$g0_dw_w
  Q0_dw_w <- densities_and_cond_expectation$Q0_dw_w
  q_w <- densities_and_cond_expectation$q_w
  
  # First derivative and second derivatives. Might need to fix it for negative values of alpha0
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
  
  list(first_derivative = first_derivative, second_derivative = second_derivative)
}
