source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)
library(boot)


# Specify data-generating distribution ------------------------------------

# Set of parameters 1
# lambda <- 2; alpha0 <- 2; beta0 <- -3; beta1 <- 1.5; beta2 <- 1
# beta <- 0.25
# gamma <- 0.125
# kappa <- 1 / (2 * (gamma + 1 - beta))
# Set of parameters 2
# lambda <- 2; alpha0 <- 4; beta0 <- -3; beta1 <- 1.5; beta2 <- 0.5
# kappa <- 5 / 4
# beta <- 2 - kappa
# gamma <- 1 - kappa / 2
# Set of parameters 3
lambda <- 2; alpha0 <- 4; beta2 <- -3; beta0 <- -1; beta1 <- 1
beta <- 7/8; gamma <- 5/16

# Define inference functions ----------------------------------------------

# Define finite difference functions
finite_difference_with_bootstrap_of_targeting_step <- function(data, delta, Delta){
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(data, delta + Delta)$Psi_n
  TMLE_Psi_0_delta.bootstrap <- TMLE_EY1_bootstrap_speedglm(data, delta, nb_boostrap_samples = 1000)
  Psi_n_delta <- TMLE_Psi_0_delta.bootstrap$Psi_n
  list(fin_diff = (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta,
       shapiro.p_value = TMLE_Psi_0_delta.bootstrap$shapiro.p_value)
}

finite_difference_for_boot <- function(data, indices, delta, Delta){
  replicate.data <- data[indices, ]
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(replicate.data, delta + Delta)$Psi_n
  Psi_n_delta <- TMLE_EY1_speedglm(replicate.data, delta)$Psi_n
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}

# Compute linear fits for moving windows
moving_windows.fits <- function(results_df, min_gap = 4){
  regression_df <- data.frame(log_delta = log(results_df$delta) / log(10), 
                              log_fin_diff = log(abs(results_df$fin_diff)) / log(10),
                              p_value = results_df$p_value)
  regression_df <- subset(regression_df, is.finite(log_fin_diff) & p_value != 0)
  
  regressions.results <- vector()
  sorted_log_delta <- sort(unique(regression_df$log_delta))
  min_gap <- min(4, length(sorted_log_delta) - 1)
  for(id_lower_bound in 1:(length(sorted_log_delta) - min_gap)){
    for(id_upper_bound in (id_lower_bound + min_gap):length(sorted_log_delta)){
      lower_bound <- sorted_log_delta[id_lower_bound]
      upper_bound <- sorted_log_delta[id_upper_bound]
      # cat('Lower_bound = ', lower_bound, ' and upper bound = ', upper_bound)
      in_window_subset <- subset(regression_df, log_delta >= lower_bound &
                                   log_delta <= upper_bound)
      # cat('Length of in window subset:', nrow(in_window_subset))
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

# Main directions fits: perform moving windows regressions and score them
# This includes nice plotting
main_directions_fits <- function(results_df, min_gap = 5){
  
  regression_df <- data.frame(log_delta = log(results_df$delta) / log(10), 
                              log_fin_diff = log(abs(results_df$fin_diff)) / log(10),
                              p_value = results_df$p_value)
  regression_df <- subset(regression_df, is.finite(log_fin_diff) & p_value != 0)
  
  bounds <- matrix(0, nrow  = length(unique(regression_df$log_delta)),
                   ncol = length(unique(regression_df$log_delta)))
  
  regressions.results <- vector()
  sorted_log_delta <- sort(unique(regression_df$log_delta))
  min_gap <- min(min_gap, length(sorted_log_delta) - 1)
  for(id_lower_bound in 1:(length(sorted_log_delta) - min_gap)){
    for(id_upper_bound in (id_lower_bound + min_gap):length(sorted_log_delta)){
      lower_bound <- sorted_log_delta[id_lower_bound]
      upper_bound <- sorted_log_delta[id_upper_bound]
      bounds[id_lower_bound, id_upper_bound] <- 1
      # cat('Lower_bound = ', lower_bound, ' and upper bound = ', upper_bound)
      in_window_subset <- subset(regression_df, log_delta >= lower_bound &
                                   log_delta <= upper_bound)
      # cat('Length of in window subset:', nrow(in_window_subset))
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
  regressions.results <- data.frame(regressions.results)
  
  # Scoring of fits
  fits.penalty <- vector()
  for(i in 1:nrow(regressions.results)){
    fits.penalty <- c(fits.penalty,
                      3 * rank(-regressions.results$avg_log_p_value)[i]^2 +
                        rank(-regressions.results$r_squared)[i]^2 +
                        500 * as.numeric(regressions.results$slope[i] < -1))
  }
  
  
  regressions.results <- cbind(regressions.results,
                               y.start = regressions.results$lower_bound * regressions.results$slope + 
                                 regressions.results$intercept,
                               y.end = regressions.results$upper_bound * regressions.results$slope + 
                                 regressions.results$intercept,
                               fits.penalty = fits.penalty)
  # Lines with true slope
  intercepts.range <- range(c(regressions.results$y.end, regressions.results$y.start))
  true_lines.df <- data.frame(intercept = seq(from = intercepts.range[1] - 2, 
                                              to = intercepts.range[2] + 2, 
                                              length = 15) + beta * min(regressions.results$lower_bound),
                              slope = rep(-beta, 15))
  
  beta_hat <- -subset(regressions.results, rank(fits.penalty, ties.method = 'first') <= 1, select = slope)[[1]]
  
  fits.plot <- ggplot(data = regressions.results, aes(x = lower_bound, y = y.start)) +
    geom_segment(aes(xend = upper_bound, yend = y.end, colour = avg_log_p_value,
                     alpha = r_squared)) +
    geom_segment(data = subset(regressions.results, rank(fits.penalty, ties.method = 'first') <= 10 ), 
                 aes(x = lower_bound, y = y.start, xend = upper_bound, yend = y.end), colour = 'yellow') +
    geom_segment(data = subset(regressions.results, rank(fits.penalty, ties.method = 'first') <= 1), 
                 aes(x = lower_bound, y = y.start, xend = upper_bound, yend = y.end), colour = 'red') +
    geom_abline(data = true_lines.df, mapping = aes(intercept = intercept,
                                                    slope = slope),
                colour = 'black', linetype = 'dotted') +
    annotate("text", x = -2, y = mean(intercepts.range), 
             label = paste('beta_hat = ', round(beta_hat, 4), '\nbeta = ', beta, sep = '')) +
    xlab(expression(log[10](delta))) + 
    ylab(expression(log[10](widehat(Delta * b[n])(delta)))) +
    ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                         group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                       list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))
  
  print(fits.plot)
  
  list(fits.plot = fits.plot, beta_hat = beta_hat)
}

# Set up cluster
cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()), outfile = '')
registerDoParallel(cl)

# Simulate data
n <- 10^4
observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, n)

# Set up tasks
deltas <- 10^seq(from = -5, to = -0.8, by = 0.1)
Delta.delta_rates <- c(0.8, 1, 1.1, 1.375, 1.5)
jobs <- expand.grid(delta = deltas, Delta.delta_rate = Delta.delta_rates)
outlier_correction = TRUE

results <- foreach(i=1:nrow(jobs), .combine = rbind,
                   .packages = c('speedglm', 'boot'), .verbose = T, .inorder = T) %dopar% {
                     
                     delta <- jobs[i, ]$delta; Delta <- n^(-0.25) * delta^jobs[i, ]$Delta.delta_rate
                     
                     observed_data.df <- data.frame(L0 = observed_data$L0,
                                                    A0 = observed_data$A0,
                                                    L1 = observed_data$L1)
                     fin_diff <- finite_difference_for_boot(observed_data.df, 1:n, delta, Delta)
                     fin_diff.bootstrap <- boot(data = observed_data.df,
                                                statistic = finite_difference_for_boot,
                                                R = 100, sim = 'ordinary',
                                                delta = delta, Delta = Delta)$t
                     
                     # Outlier correction
                     if(outlier_correction){
                       if(fin_diff < quantile(fin_diff.bootstrap, c(0.1, 0.9))[1] |
                          fin_diff > quantile(fin_diff.bootstrap, c(0.1, 0.9))[2])
                         fin_diff <- median(fin_diff.bootstrap)
                     }
                     
                     shapiro.p_value <- 0
                     try(shapiro.p_value <- shapiro.test(fin_diff.bootstrap)$p.value)
                     if(is.null(shapiro.p_value)) shapiro.p_value <- 0
                     
                     c(delta, jobs[i, ]$Delta.delta_rate, fin_diff, shapiro.p_value)
                   }
stopCluster(cl)

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

moving_windows.fits.results <- data.frame(moving_windows.fits(subset(results_df, log_p_value > log_p_value_cutoff)))
fits.scores <- vector()
for(i in 1:nrow(moving_windows.fits.results)){
  fits.scores <- c(fits.scores,
                   rank(-moving_windows.fits.results$avg_log_p_value)[i]^2 +
                     rank(-moving_windows.fits.results$r_squared)[i]^2 +
                     rank(-moving_windows.fits.results$r_squared)[i]^2 +
                     0.5 * rank(-moving_windows.fits.results$nb_points)[i]^2)
}
moving_windows.fits.results <- cbind(moving_windows.fits.results, fits.scores = fits.scores)
beta_hat <- -as.numeric(subset(moving_windows.fits.results, fits.scores == min(fits.scores), select = slope)[1])
intercept <- as.numeric(subset(moving_windows.fits.results, fits.scores == min(fits.scores), select = intercept)[1])

print(fin_diffs.plot + geom_abline(intercept = intercept, slope = -beta_hat) + 
        geom_abline(intercept = intercept - 1, slope = -beta, linetype = 'dotted'))

# # Is there a perfectly linear region with more than 15 points?
# cat('About to compute the linear fits for moving windows\n')
# fits <- data.frame(moving_windows.fits(results_df, min_gap = 10))
# cat('Done with computing linear fits for moving windows\n')
# if(any(fits$r_squared > 0.99)){
#   beta_hat <- -median(subset(fits, r_squared > 0.99 & upper_bound > -2, select = slope)[[1]])
#   intercept <- median(subset(fits, r_squared > 0.99 & upper_bound > -2, select = intercept)[[1]])
#   cat('Perfectly linear region with more than 10 points found')
# }else{
#   cat('No perfectly linear region with more than 10 points found')
#   # Fit p-value cutoff
#   lambda <- 0.01
#   log_p_value_cutoff <- (1 - lambda) * max(results_df$log_p_value) + lambda * min(results_df$log_p_value)
#   lower_log_delta_bound <- Inf; upper_log_delta_bound <- Inf
#   while(nrow(subset(results_df, log_p_value > log_p_value_cutoff)) < 10 * length(Delta.delta_rates) &
#         log_p_value_cutoff > -15 |
#         lower_log_delta_bound == upper_log_delta_bound){
#     lambda <- lambda + 0.01
#     log_p_value_cutoff <- (1 - lambda) * max(results_df$log_p_value) + lambda * min(results_df$log_p_value)
#     lower_log_delta_bound <- min(subset(results_df, log_p_value > log_p_value_cutoff, select = log_delta))
#     upper_log_delta_bound <- max(subset(results_df, log_p_value > log_p_value_cutoff, select = log_delta))
#     if(lambda > 0.4 & lower_log_delta_bound != upper_log_delta_bound) break
#   }
#   
#   plot(results_df$log_delta, results_df$log_p_value)
#   abline(h = log_p_value_cutoff)
#   abline(v = lower_log_delta_bound); abline(v = upper_log_delta_bound)
#   
#   moving_windows.fits.results <- data.frame(moving_windows.fits(subset(results_df, log_p_value > log_p_value_cutoff)))
#   fits.scores <- vector()
#   for(i in 1:nrow(moving_windows.fits.results)){
#     fits.scores <- c(fits.scores,
#                      rank(-moving_windows.fits.results$avg_log_p_value)[i]^2 +
#                        rank(-moving_windows.fits.results$r_squared)[i]^2 +
#                        rank(-moving_windows.fits.results$r_squared)[i]^2 +
#                        0.5 * rank(-moving_windows.fits.results$nb_points)[i]^2)
#   }
#   moving_windows.fits.results <- cbind(moving_windows.fits.results, fits.scores = fits.scores)
#   beta_hat <- -as.numeric(subset(moving_windows.fits.results, fits.scores == min(fits.scores), select = slope)[1])
#   intercept <- as.numeric(subset(moving_windows.fits.results, fits.scores == min(fits.scores), select = intercept)[1])
# }
# print(fin_diffs.plot + geom_abline(intercept = intercept, slope = -beta_hat) + 
#         geom_abline(intercept = intercept - 1, slope = -beta, linetype = 'dotted'))


# Find main direction of fits ---------------------------------------------

main_directions_fits.results <- main_directions_fits(results_df, min_gap = 7)
print(main_directions_fits$results.fits.plot)
cat('beta_hat =', main_directions_fits.results$beta_hat)

# Save the finite differences
save('results_df', file = paste('finite_differences_',
                                floor(n), '_', lambda, '_', alpha0, '_', beta0, '_', beta1, '_', beta2, '_',
                                as.integer(as.POSIXct( Sys.time() )),
                                '.results',
                                sep = ''))