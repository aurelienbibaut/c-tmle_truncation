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

# Define finite difference functions
finite_difference_for_boot <- function(data, indices, delta, Delta){
  replicate.data <- data[indices, ]
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(replicate.data, delta + Delta)$Psi_n
  Psi_n_delta <- TMLE_EY1_speedglm(replicate.data, delta)$Psi_n
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}

# Define find_beta, a function that infers beta based on finite differences
find_beta <- function(observed_data, deltas, Delta.delta_rates){
  jobs <- expand.grid(delta = deltas, Delta.delta_rate = Delta.delta_rates)
  
  # Set up cluster
  cl <- makeCluster(getOption("cl.cores", 32), outfile = '')
  registerDoParallel(cl)
  
  results <- foreach(i=1:nrow(jobs), .combine = rbind, .export = c('finite_difference_for_boot', 'beta', 'gamma', 'kappa',
                                                                   'lambda', 'alpha0', 'beta0', 'beta1', 'beta2',
                                                                   'TMLE_EY1_speedglm', 'expit', 'logit', 'g_to_g_delta'),
                     .packages = c('speedglm', 'boot'), .verbose = T, .inorder = T) %dopar% {
                       
                       delta <- jobs[i, ]$delta; Delta <- n^(-0.25) * delta^jobs[i, ]$Delta.delta_rate
                       
                       observed_data.df <- data.frame(L0 = observed_data$L0,
                                                      A0 = observed_data$A0,
                                                      L1 = observed_data$L1)
                       fin_diff <- finite_difference_for_boot(observed_data.df, 1:n, delta, Delta)
                       fin_diff.bootstrap <- boot(data = observed_data.df,
                                                  statistic = finite_difference_for_boot,
                                                  R = 1000, sim = 'ordinary',
                                                  delta = delta, Delta = Delta)$t
                       shapiro.p_value <- 0
                       try(shapiro.p_value <- shapiro.test(fin_diff.bootstrap)$p.value)
                       if(is.null(shapiro.p_value)) shapiro.p_value <- 0
                       
                       c(delta, jobs[i, ]$Delta.delta_rate, fin_diff, shapiro.p_value)
                     }
  row.names(results) <- NULL
  colnames(results) <- c('delta', 'Delta.delta_rate', 'fin_diff', 'p_value')
  
  results_df <- as.data.frame(results)
  
  # Make plot of finite differences
  fin_diffs.plot <- ggplot(results_df, 
                           aes(x = log(delta) / log(10), 
                               y = log(abs(fin_diff)) / log(10),
                               colour = factor(Delta.delta_rate))) + 
    geom_line() + geom_point(aes(size = log(p_value) / log(10))) + 
    geom_abline(intercept = -0.5, slope = -beta, linetype = 'dotted') +
    geom_abline(intercept = -1, slope = -beta, linetype = 'dotted') +
    geom_abline(intercept = -1.5, slope = -beta, linetype = 'dotted') +
    xlab(expression(log[10](delta))) + 
    ylab(expression(log[10](widehat(Delta * b[n])(delta)))) +
    ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                         group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                       list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))
  
  print(fin_diffs.plot)
  
  # Find beta
  results_df <- cbind(results_df, log_delta = log(results_df$delta) / log(10),
                      log_p_value = log(results_df$p_value) / log(10))
  
  # Find the right window of log_deltas where to fit the linear regression and do it
  moving_windows.fits <- function(results_df){
    regression_df <- data.frame(log_delta = log(results_df$delta) / log(10), 
                                log_fin_diff = log(abs(results_df$fin_diff)) / log(10),
                                p_value = results_df$p_value)
    regression_df <- subset(regression_df, is.finite(log_fin_diff) & p_value != 0)
    
    regressions.results <- vector()
    sorted_log_delta <- sort(unique(regression_df$log_delta))
    for(id_lower_bound in 1:(length(sorted_log_delta) - 5)){
      for(id_upper_bound in (id_lower_bound + 5):length(sorted_log_delta)){
        lower_bound <- sorted_log_delta[id_lower_bound]
        upper_bound <- sorted_log_delta[id_upper_bound]
        in_window_subset <- subset(regression_df, log_delta >= lower_bound &
                                     log_delta <= upper_bound)
        lm.out <- lm(log_fin_diff ~ log_delta, in_window_subset)
        regressions.results <- rbind(regressions.results,
                                     c(lower_bound = lower_bound, 
                                       upper_bound = upper_bound,
                                       r_squared  = summary(lm.out)$r.squared,
                                       nb_points  = nrow(in_window_subset),
                                       avg_log_p_value = mean(log(in_window_subset$p_value) / log(10)),
                                       intercept = as.numeric(lm.out$coefficients[1]),
                                       slope = as.numeric(lm.out$coefficients[2]),
                                       true_beta = beta))
      }
    }
    regressions.results
  }
  
  
  
  # Fit p-value cutoff
  lambda <- 0.01
  log_p_value_cutoff <- (1 - lambda) * max(results_df$log_p_value) + lambda * min(results_df$log_p_value)
  
  while(nrow(subset(results_df, log_p_value > log_p_value_cutoff)) < 10 * length(Delta.delta_rates) &
        log_p_value_cutoff > -15){
    lambda <- lambda + 0.01
    log_p_value_cutoff <- (1 - lambda) * max(results_df$log_p_value) + lambda * min(results_df$log_p_value)
    lower_log_delta_bound <- min(subset(results_df, log_p_value > log_p_value_cutoff, select = log_delta))
    upper_log_delta_bound <- max(subset(results_df, log_p_value > log_p_value_cutoff, select = log_delta))
    if(lambda > 0.4) break
  }
  
  plot(results_df$log_delta, results_df$log_p_value)
  abline(h = log_p_value_cutoff)
  abline(v = lower_log_delta_bound); abline(v = upper_log_delta_bound)
  
  moving_windows.fits.results <- data.frame(moving_windows.fits(subset(results_df, 
                                                                       log_p_value > log_p_value_cutoff)))
  fits.scores <- vector()
  for(i in 1:nrow(moving_windows.fits.results)){
    fits.scores <- c(fits.scores,
                     rank(-moving_windows.fits.results$avg_log_p_value)[i]^2 +
                       rank(-moving_windows.fits.results$r_squared)[i]^2 +
                       rank(-moving_windows.fits.results$r_squared)[i]^2 +
                       0.5 * rank(-moving_windows.fits.results$nb_points)[i]^2)
  }
  moving_windows.fits.results <- cbind(moving_windows.fits.results, 
                                       fits.scores = fits.scores)
  beta_hat <- -as.numeric(subset(moving_windows.fits.results, 
                                 fits.scores == min(fits.scores), 
                                 select = slope)[1])
  intercept <- as.numeric(subset(moving_windows.fits.results, 
                                 fits.scores == min(fits.scores), 
                                 select = intercept)[1])
  
  list(beta_hat = beta_hat, 
       beta_plot = moving_windows.fits + geom_abline(intercept = intercept, slope = -beta_hat))
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
  
  results <- cbind(results, in_range = results[, 'delta'] <= delta_max & results[, 'delta'] >= delta_min)
  
  regression_df <- subset(data.frame(results), delta <= delta_max & delta >= delta_min)
  regression_df <- cbind(regression_df, log_delta = log(regression_df$delta) / log(10),
                         log_var_IC = log(regression_df$var_IC_delta) / log(10))
  lts_fit <- ltsReg(regression_df$log_delta, regression_df$log_var_IC)
  ols_fit <- lm(log_var_IC ~ log_delta, regression_df)
  results_df <- as.data.frame(results)
  
  sigmas_ns.plot <- ggplot(results_df, aes(x = log(delta) / log(10), 
                                           y = log(var_IC_delta) / log(10))) + 
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