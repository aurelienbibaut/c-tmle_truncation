source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)


# Specify data-generating distribution ------------------------------------

# Set of parameters 1
# lambda <- 2; alpha0 <- 2; beta0 <- -3; beta1 <- 1.5; beta2 <- 1
# beta <- 0.25
# gamma <- 0.125
# kappa <- 1 / (2 * (gamma + 1 - beta))
# Set of parameters 2
lambda <- 2; alpha0 <- 4; beta0 <- -3; beta1 <- 1.5; beta2 <- 0.5
kappa <- 5 / 4
beta <- 2 - kappa
gamma <- 1 - kappa / 2


# Define inference functions ----------------------------------------------

# Define finite difference function
finite_difference <- function(data, delta, n){
  Delta <- n^(-0.25) * delta^((beta + 1 - gamma) / 2)
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(data, delta + Delta)$Psi_n
  Psi_n_delta <- TMLE_EY1_speedglm(data, delta)$Psi_n
  
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}

# Define find_beta, a function that infers beta based on finite differences
find_beta <- function(observed_data, ns, etas){
  jobs <- expand.grid(n = ns, eta = etas)
  
  # Set up cluster
  cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
  registerDoParallel(cl)
  
  results <- foreach(i=1:nrow(jobs), .combine = rbind, .export = c('finite_difference', 'beta', 'gamma', 'kappa',
                                                                   'lambda', 'alpha0', 'beta0', 'beta1', 'beta2',
                                                                   'TMLE_EY1_speedglm', 'expit', 'logit', 'g_to_g_delta'),
                     .packages = c('speedglm', 'robustbase', 'ggplot2'), .verbose = T, .inorder = T) %dopar% {
                       
                       n <- jobs[i, ]$n; eta <- jobs[i, ]$eta
                       delta <- n^(-1 / (2 * eta * (gamma + 1 - beta)))
                       
                       if(n == max(ns)){
                         indices <- t(replicate(1000, sample(1:max(ns), n, replace = T)))
                       }else if(n < max(ns)){
                         indices <- t(replicate(1000, sample(1:max(ns), n, replace = F)))
                       }else{
                         indices <- t(replicate(1000, sample(1:max(ns), n, replace = T)))
                       }
                       
                       cat('Job ', i, ', dim(indices): ')
                       print(dim(indices))
                       
                       fin_diffs <- apply(indices, 1, function(y) finite_difference(lapply(observed_data, 
                                                                                           function(x) x[y]), delta, n))
                       fin_diff <- median(fin_diffs)
                       
                       shapiro.p_value <- shapiro.test(fin_diffs)$p.value
                       
                       cat('Job ', i, ', n = ', n, ', eta = ', eta, 
                           ', delta = ', delta, ', fin_diff = ', fin_diff, '\n')
                       c(n, eta, delta, fin_diff, shapiro.p_value)
                     }
  
  
  stopCluster(cl)
  
  row.names(results) <- NULL
  colnames(results) <- c('n', 'eta', 'delta', 'fin_diff', 'shapiro.p_value')
  
  results_df <- as.data.frame(results)
  results_df <- transform(results_df, n = as.numeric(as.character(n)), 
                          eta = as.numeric(as.character(eta)),
                          delta = as.numeric(as.character(delta)),
                          fin_diff = as.numeric(as.character(fin_diff)),
                          shapiro.p_value = as.numeric(as.character(shapiro.p_value)))
  
  regression_df <- as.data.frame(cbind(log_fin_diff = log(abs(results_df$fin_diff))/log(10), 
                                       log_delta = log(results_df$delta) / log(10),
                                       shapiro.p_value = results_df$shapiro.p_value,
                                       eta = results_df$eta,
                                       n = results_df$n))
  line_fit <- lm(formula = log_fin_diff ~ log_delta, regression_df)
  lts.line_fit <- ltsReg(regression_df$log_delta, regression_df$log_fin_diff)
  
  trimmed_regression_df <- subset(regression_df, 
                                  shapiro.p_value  >= max(quantile(regression_df$shapiro.p_value)[2], 1e-5) & 
                                    log_delta <= quantile(regression_df$log_delta)[5])
  lts.p_val_trimmed.line_fit <- NULL
  try(lts.p_val_trimmed.line_fit <-ltsReg(trimmed_regression_df$log_delta, trimmed_regression_df$log_fin_diff))
  
  if(!is.null(lts.p_val_trimmed.line_fit)){
    if(max(ns) < 4e4){
      beta_hat <- max(-lts.line_fit$raw.coefficients[2],
                      lts.p_val_trimmed.line_fit$raw.coefficients[2])
    }else{
      beta_hat <- -line_fit$coefficients["log_delta"]
    }
  }else{
    if(max(ns) < 4e4){
      beta_hat <- -lts.line_fit$raw.coefficients[2]
    }else{
      beta_hat <- -line_fit$coefficients["log_delta"]
    }
  }
  
  
  fin_diffs.all_etas.plot <- ggplot(results_df, 
                                    aes(x = log(delta) / log(10), 
                                        y = log(abs(fin_diff)) / log(10),
                                        colour = factor(eta),
                                        label = round(log(n) / log(10), 1))) + 
    geom_line() + geom_text() + geom_point(aes(size = shapiro.p_value)) +
    geom_point(data = trimmed_regression_df, aes(x = log_delta, y = log_fin_diff), shape = 0, size = 10) +
    geom_abline(intercept = -2, slope = -beta, linetype = 'dashed') +
    geom_abline(intercept = -2.1, slope = -beta, linetype = 'dashed') +
    geom_abline(intercept = -2.2, slope = -beta, linetype = 'dashed') +
    geom_abline(intercept = line_fit$coefficients["(Intercept)"],
                slope = line_fit$coefficients["log_delta"], colour = 'blue', size = 1) +
    geom_abline(intercept = lts.line_fit$raw.coefficients[1],
                slope = lts.line_fit$raw.coefficients[2], colour = 'green', size = 1) +
    geom_abline(intercept = lts.p_val_trimmed.line_fit$raw.coefficients[1],
                slope = lts.p_val_trimmed.line_fit$raw.coefficients[2], colour = 'red', size = 1) +
    annotate("text", x = quantile(regression_df$log_delta)[2],
             y = quantile(regression_df$log_fin_diff)[2],
             label = paste(c("beta_hat = ", round(beta_hat, 3)), collapse = '')) + 
    xlab(expression(log[10](delta))) + 
    ylab(expression(log[10](widehat(Delta * b[n])(delta)))) +
    ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                         group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                       list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))
  
  list(beta_hat = beta_hat, beta_plot = fin_diffs.all_etas.plot)
}

# Find gamma from empirical variance
find_gamma <- function(observed_data, deltas){
  # Set up cluster
  cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
  registerDoParallel(cl)
  
  results <- foreach(delta = deltas, .combine = rbind, .export = c('beta', 'gamma', 'kappa',
                                                                   'lambda', 'alpha0', 'beta0', 'beta1', 'beta2',
                                                                   'TMLE_EY1_speedglm', 'expit', 'logit', 'g_to_g_delta'),
                     .packages = c("speedglm", "robustbase"), .verbose = T, .inorder = T) %dopar% {
                       
                       cat('In find gamma, delta = ', delta, '\n')
                       
                       TMLE_delta.results <- TMLE_EY1_speedglm(observed_data, delta)
                       var_IC_delta <- TMLE_delta.results$sigma_n
                       c(delta, var_IC_delta)
                     }
  stopCluster(cl)
  row.names(results) <- NULL
  colnames(results) <- c('delta', 'var_IC_delta')
  
  log_var_IC_delta <- log(results[, 'var_IC_delta']) / log(10)
  log_var_min <- 0.9 * min(log_var_IC_delta) + 0.1 * max(log_var_IC_delta)
  log_var_max <- 0.1 * min(log_var_IC_delta) + 0.9 * max(log_var_IC_delta)
  delta_min <- max(results[log_var_IC_delta >= log_var_max, 'delta'])
  delta_max <- max(results[log_var_IC_delta >= log_var_min, 'delta'])
  
  results <- cbind(results, in_range = results_n[, 'delta'] <= delta_max & results_n[, 'delta'] >= delta_min)
  
  regression_df <- subset(data.frame(results), delta <= delta_max & delta >= delta_min)
  regression_df <- cbind(regression_df, log_delta = log(regression_df$delta) / log(10),
                         log_var_IC = log(regression_df$var_IC_delta) / log(10))
  lts_fit <- ltsReg(regression_df$log_delta, regression_df$log_var_IC)
  ols_fit <- lm(log_delta ~ log_var_IC, regression_df)
  results_df <- as.data.frame(results)
  
  sigmas_ns.plot <- ggplot(results_df, aes(x = log(delta) / log(10), 
                                           y = log(var_IC_delta) / log(10), 
                                           colour = factor(n))) + 
    geom_line(aes(size = factor(in_range))) + 
    geom_abline(intercept = -1, slope = - 2 * gamma) +
    geom_hline(yintercept =  log_var_max, colour = 'red') +
    geom_hline(yintercept =  log_var_min, colour = 'red') +
    geom_vline(xintercept = log(delta_min) / log(10), colour = 'red') +
    geom_vline(xintercept = log(delta_max) / log(10), colour = 'red') +
    geom_abline(intercept = lts_fit$raw.coefficients[1],
                slope = lts_fit$raw.coefficients[2], colour = 'green') +
    geom_abline(intercept = ols_fit$coefficients[1],
                slope = ols_fit$coefficients[2], colour = 'blue') +
    xlim(-6, 0) +
    xlab(expression(log[10](delta))) + 
    ylab(expression(log[10](sigma[n](delta)))) +
    annotate("text", x = -5, y = -0.2, label = paste(c("gamma_hat = ", 
                                                       round(-lts_fit$raw.coefficients[2] / 2, 3)), collapse = '')) + 
    ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                         group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                       list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))
  list(gamma_hat = -lts_fit$raw.coefficients[2] / 2, gamma_plot = sigmas_ns.plot)
}
# debug(find_gamma)

# Define algorithm tuning parameters --------------------------------------
etas <- seq(from = 2, to = 5, by = 0.5)
finite_diffs <- vector(); finite_diffs_bias <- vector(); true_finite_diffs_bias <- vector()
ns <- floor(10^seq(from = 2, to = 4.5, by = 0.5))
deltas <- 10^seq(from = -1, to = -7.5, by = -0.05)

# Run inference algorithm -------------------------------------------------
observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, max(ns))

find_beta.result <- find_beta(observed_data, ns, etas)
cat('beta_hat = ', find_beta.result$beta_hat, '\n')
print(find_beta.result$beta_plot)

find_gamma.result <- find_gamma(observed_data, deltas)
cat('gamma_hat = ', find_gamma.result$gamma_hat, '\n')
print(find_gamma.result$gamma_plot)

# Estimated rate
estimated_rate <- -1 / (2 * (find_gamma.result$gamma_hat + 1 - find_beta.result$beta_hat))
grid.arrange(find_beta.result$beta_plot, find_gamma.result$gamma_plot,
             nrow = 2, ncol = 1,
             bottom = textGrob(paste(c("Estimated rate: ", estimated_rate, 
                                       "\nTrue rate: ", -1 / (2 * (gamma + 1 - beta))), collapse = ''),
                               gp=gpar(fontsize=20,font=3)))

cat('Estimated rate = ', estimated_rate, '\n')