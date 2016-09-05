source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')

# Set of parameters 1
lambda <- 2; alpha0 <- 2; beta0 <- -3; beta1 <- 1.5; beta2 <- 1
# beta <- 0.25
# gamma <- 0.125
# kappa <- 1 / (2 * (gamma + 1 - beta))
# Set of parameters 2
# lambda <- 2; alpha0 <- 4; beta0 <- -3; beta1 <- 1.5; beta2 <- 0.5
# kappa <- 5 / 4
# beta <- 2 - kappa
# gamma <- 1 - kappa / 2


# Compute key integral
compute_bias.key_integral <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, delta){
  # Generate and get densities and conditional expectations functions
  densities_and_cond_expectation <- generate_densities(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0)
  g0_dw_w <- densities_and_cond_expectation$g0_dw_w
  Q0_dw_w <- densities_and_cond_expectation$Q0_dw_w
  q_w <- densities_and_cond_expectation$q_w
  
  # Define integrand such that \int integrand(w) dw = Psi_0(\delta) and integrate it
  integrand <- Vectorize(function(w) q_w(w) * g0_dw_w(w) * Q0_dw_w(w))
  integrate(integrand, lower =  -Inf, upper = logit(delta) / alpha0)$value
}

compute_var.key_integral <- function(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0, delta){
  # Generate and get densities and conditional expectations functions
  densities_and_cond_expectation <- generate_densities(type, positivity_parameter, alpha0, beta0, beta1, beta2, d0)
  g0_dw_w <- densities_and_cond_expectation$g0_dw_w
  Q0_dw_w <- densities_and_cond_expectation$Q0_dw_w
  q_w <- densities_and_cond_expectation$q_w
  
  # Define integrand such that \int integrand(w) dw = Psi_0(\delta) and integrate it
  integrand <- Vectorize(function(w) q_w(w) * g0_dw_w(w) * Q0_dw_w(w) * (1 - Q0_dw_w(w)))
  integrate(integrand, lower =  -Inf, upper = logit(delta) / alpha0)$value
}

deltas <- 10^seq(from = -8, to = -1, length = 1000)

bias.key_integrals <- sapply(deltas, function(delta) compute_bias.key_integral("L0_exp", lambda, alpha0, beta0, beta1, beta2,
                                                                       alwaysTreated0, delta))
var.key_integrals <- sapply(deltas, function(delta) compute_var.key_integral("L0_exp", lambda, alpha0, beta0, beta1, beta2,
                                                                               alwaysTreated0, delta))

# fin_diffs <- Psi_0_deltas[2:length(deltas)] - Psi_0_deltas[1:(length(deltas) - 1)] / 
#   (deltas[2:length(deltas)] - deltas[1:(length(deltas) - 1)])
bias.kappa <- (alpha0 + max(beta2, 0) + lambda^-1) / alpha0
var.kappa <- (alpha0 + abs(beta2) + lambda^-1) / alpha0
beta <- 2 - bias.kappa
gamma <- 1 - var.kappa / 2

C_b <- lambda^- 1 /  2 * exp((beta0 + beta1) * (beta2 > 0)) * 1 / (lambda^-1 + max(beta2, 0) + alpha0)
C_sigma2 <- lambda^-1 / 2 * exp(sign(beta2) * (beta0 + beta1)) / (lambda^-1 + alpha0 + abs(beta2))

# Optimal constant
abs_optimal_rate <- 1 / (gamma + 1 - beta)
optimal_const <- (C_sigma2 / C_b^2 * gamma / (1 - beta))^abs_optimal_rate
delta_1000 <- optimal_const * 1000^-abs_optimal_rate
cat('Optimal n rate: ', -abs_optimal_rate,
    '\nOptimal constant: ', optimal_const,
    '\nThus delta_1000 = ', delta_1000, '\n')

par(mfrow = c(1, 2))
plot(log(deltas) / log(10), log(abs(bias.key_integrals)) / log(10), main = "Bias' key integral")
abline(a = log(C_b) / log(10), 
       b = bias.kappa)

plot(log(deltas) / log(10), log(abs(var.key_integrals)) / log(10), main = "Variance' key integral")
abline(a = log(C_sigma2) / log(10),
       b = var.kappa)
