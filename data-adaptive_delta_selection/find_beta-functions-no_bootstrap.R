source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)
library(boot); library(segmented)

# Define finite difference functions
finite_difference <- function(observed_ta, delta, Delta){
  Psi_n_delta_plus_Delta <- TMLE_EY1_speedglm(observed_data, delta + Delta)$Psi_n
  Psi_n_delta <- TMLE_EY1_speedglm(observed_data, delta)$Psi_n
  (Psi_n_delta_plus_Delta - Psi_n_delta) / Delta
}


# Compute finite differences for a different rates
compute_finite_difference <- function(type, lambda, alpha0, beta0, beta1, beta2, n, Delta.delta_rates){
  # Simulate data
  observed_data <- generate_data(type, lambda, alpha0, beta0, beta1, beta2, n)
  
  # Set up tasks
  deltas <- 10^seq(from = -5, to = -0.8, by = 0.1)
  
  jobs <- expand.grid(delta = deltas, Delta.delta_rate = Delta.delta_rates)
  outlier_correction = TRUE
  
  results <- foreach(i=1:nrow(jobs), .combine = rbind,
                     .packages = c('speedglm', 'boot'),
                     .export = c('TMLE_EY1_speedglm', 'expit', 'logit', 'g_to_g_delta', 'finite_difference'),
                     .verbose = T, .inorder = T) %dopar% {
                       
                       delta <- jobs[i, ]$delta; Delta <- n^(-0.25) * delta^jobs[i, ]$Delta.delta_rate
                       
                       observed_data.df <- data.frame(L0 = observed_data$L0,
                                                      A0 = observed_data$A0,
                                                      L1 = observed_data$L1)
                       fin_diff <- finite_difference(observed_data.df, delta, Delta)
                       
                       c(delta, jobs[i, ]$Delta.delta_rate, fin_diff)
                     }
  
  row.names(results) <- NULL
  colnames(results) <- c('delta', 'Delta.delta_rate', 'fin_diff')
  
  results_df <- as.data.frame(results)
  
  results_df <- cbind(results_df, log_delta = log(results_df$delta) / log(10),
                      log_fin_diff = log(abs(results_df$fin_diff)) / log(10))
  results_df
}


# Plot finite differences
plot_finite_differences <- function(results_df, beta){
  # Make plot of finite differences
  fin_diffs.plot <- ggplot(results_df, 
                           aes(x = log_delta, 
                               y = log_fin_diff,
                               colour = factor(Delta.delta_rate))) + 
    geom_line() + 
    geom_abline(intercept = -0.5, slope = -beta, linetype = 'dotted') +
    geom_abline(intercept = -1, slope = -beta, linetype = 'dotted') +
    geom_abline(intercept = -1.5, slope = -beta, linetype = 'dotted') +
    xlab(expression(log[10](delta))) + 
    ylab(expression(log[10](widehat(Delta * b[n])(delta)))) +
    ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                         group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                       list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))
  
  fin_diffs.plot
}

# Fit broken line on a log(fin_diff) against log(delta) plot. One Delta.delta_rate at a time
# Fit a broken line to the in-range area of the curve
fit_broken_line <- function(results_df, nb_breakpoints = 2, beta, Delta.delta_rate, plotting = F){
  
  initial_breakpoints <- quantile(results_df$log_delta, (1:nb_breakpoints) / (nb_breakpoints + 1))
  
  lm.out <- lm(log_fin_diff ~ log_delta, results_df)
  segmented.fit.sucess <- F; nb_its <- 0
  while(!segmented.fit.sucess){
    try.out <- try(segmented.out <- segmented.lm(lm.out, seg.Z=~log_delta, psi = initial_breakpoints))
    if(paste(class(try.out), collapse = '') != "try-error") segmented.fit.sucess <- T
    nb_its <- nb_its + 1
    if(nb_its > 15) stop('Could not fit broken line')
  }
  broken_line_df <- as.data.frame(cbind(intercept = as.vector(intercept(segmented.out)$log_delta),
                                        slope = as.vector(slope(segmented.out)$log_delta[, 1])))
  
  fitted_breakpoints <- c(min(results_df$log_delta), 
                          as.vector(segmented.out$psi[, 2]), 
                          max(results_df$log_delta))
  
  broken_line_points_df <- matrix(0, ncol = 13, nrow = length(fitted_breakpoints) - 1)
  colnames(broken_line_points_df) <- c('x.start', 'x.end', 'y.start', 'y.end', 'Delta.delta_rate',
                                       'RSS', 'r_squared',
                                       'nb_points', 'leftmost', 'segment.squared_length',
                                       'slope', 'loss', 'is_best')
  for(i in 1:(length(fitted_breakpoints) - 1)){
    broken_line_points_df[i, 'x.start'] <- fitted_breakpoints[i]
    broken_line_points_df[i, 'x.end'] <- fitted_breakpoints[i + 1]
    broken_line_points_df[i, 'y.start'] <- broken_line_df$intercept[i] + broken_line_df$slope[i] * fitted_breakpoints[i]
    broken_line_points_df[i, 'y.end'] <- broken_line_df$intercept[i] + broken_line_df$slope[i] * fitted_breakpoints[i + 1]
    broken_line_points_df[i, 'Delta.delta_rate'] <- Delta.delta_rate
    segment.subset <- subset(results_df, log_delta >= fitted_breakpoints[i] &
                               log_delta <= fitted_breakpoints[i + 1])
    lm_on_subset.out <- lm(log_fin_diff ~ log_delta, segment.subset)
    broken_line_points_df[i, 'RSS'] <- sum(lm_on_subset.out$residuals^2)
    broken_line_points_df[i, 'r_squared'] <- as.numeric(summary(lm_on_subset.out)$r.squared)
    broken_line_points_df[i, 'nb_points'] <- nrow(segment.subset)
    broken_line_points_df[i, 'leftmost'] <- i
    broken_line_points_df[i, 'slope'] <- broken_line_df$slope[i]
  }
  broken_line_points_df <- as.data.frame(broken_line_points_df)
  
  RSS <- sum(broken_line_points_df$RSS)
  
  # Figure out which segment is best
  broken_line_points_df$loss <- sapply(broken_line_points_df$slope, 
                                       function(x) (abs(x) - 2 * beta)^2 + 1e6 * as.numeric(x < -1 | x > 0))
  
  # Score the segments
  segments.squared_lengths <- (broken_line_points_df$x.end - broken_line_points_df$x.start)^2 +
    (broken_line_points_df$y.end - broken_line_points_df$y.start)^2
  
  broken_line_points_df$segment.squared_length <- segments.squared_lengths
  broken_line_points_df$relative_segment_squared_length <- segments.squared_lengths / max(segments.squared_lengths)
  broken_line_points_df$linearity <- -log(1 - broken_line_points_df$r_squared) / log(10)
  broken_line_points_df$relative_linearity <- broken_line_points_df$linearity / max(broken_line_points_df$linearity)
  
  print(broken_line_points_df)
  
  # Create a plotting layer for the broken line
  broken_line.plot_layer <- NULL
  if(plotting){
    broken_line.plot_layer <- list(geom_segment(data = broken_line_points_df,
                                                mapping = aes(x = x.start, y = y.start,
                                                              xend = x.end, yend = y.end),
                                                size = 2),
                                   geom_point(data = data.frame(x = c(broken_line_points_df$x.start, broken_line_points_df$x.end),
                                                                y = c(broken_line_points_df$y.start, broken_line_points_df$y.end),
                                                                Delta.delta_rate = Delta.delta_rate),
                                              aes(x = x, y = y), shape = 4, size = 5),
                                   geom_abline(intercept = as.numeric(broken_line_df$intercept[which.min(broken_line_points_df$loss)]),
                                               slope = as.numeric(broken_line_df$slope[which.min(broken_line_points_df$loss)]),
                                               linetype = 'dotted'))
  }
  
  list(broken_line.plot_layer = broken_line.plot_layer,
       RSS = RSS,
       broken_line_points_df = broken_line_points_df)
}

# Least trimmed squares regression
lts_regression <- function(results_df, Delta.delta_rate_){
  regression_df <- subset(results_df, fin_diff !=0 & Delta.delta_rate == Delta.delta_rate_)
  lts_fit <- ltsReg(regression_df$log_delta, regression_df$log_fin_diff)
  c(intercept = as.numeric(lts_fit$coefficients[1]), slope =as.numeric(lts_fit$coefficients[2]), 
    r.squared = lts_fit$rsquared)
}

# Try it out
# Delta.delta_rates <- c(0.8, 1, 1.1, 1.375, 1.5)
# # Set of parameters 1
# # lambda <- 2; alpha0 <- 2; beta0 <- -3; beta1 <- 1.5; beta2 <- 1
# # beta <- 0.25
# # gamma <- 0.125
# # kappa <- 1 / (2 * (gamma + 1 - beta))
# # Set of parameters 2
# lambda <- 2; alpha0 <- 4; beta0 <- -3; beta1 <- 1.5; beta2 <- 0.5
# kappa <- 5 / 4
# beta <- 2 - kappa
# gamma <- 1 - kappa / 2
# # Set of parameters 3
# # lambda <- 2; alpha0 <- 4; beta2 <- -3; beta0 <- -1; beta1 <- 1
# # beta <- 7/8; gamma <- 5/16
# # n <- 1e4
# 
# fin_diffs.results <- compute_finite_difference("L0_exp", lambda, alpha0, beta0, beta1, beta2, n, Delta.delta_rates)
# fin_diffs.plot <- plot_finite_differences(subset(fin_diffs.results, fin_diff != 0), beta)
# print(fin_diffs.plot)
# # debug(fit_broken_line)
# # LTS fits
# LTS_fits <- vector()
# for(Delta.delta_rate in Delta.delta_rates){
#   LTS_fits <- rbind(LTS_fits,
#                     c(lts_regression(results_df, Delta.delta_rate), Delta.delta_rate = Delta.delta_rate))
# }
# fin_diffs.plot <- fin_diffs.plot + geom_abline(data = as.data.frame(LTS_fits),
#                                                aes(intercept = intercept, slope = slope, 
#                                                    colour = factor(Delta.delta_rate)))
# 
# broken_line_fit <- fit_broken_line(results_df = subset(fin_diffs.results, fin_diff !=0 & Delta.delta_rate == 1.1),
#                 nb_breakpoints = 1, beta, 1.1, T)
# fin_diffs.plot <- fin_diffs.plot + broken_line_fit$broken_line.plot_layer[[1]] +
#   broken_line_fit$broken_line.plot_layer[[2]] +
#   broken_line_fit$broken_line.plot_layer[[3]]
# print(fin_diffs.plot)