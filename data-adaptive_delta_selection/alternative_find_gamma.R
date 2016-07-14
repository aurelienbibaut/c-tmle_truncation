source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)
library(boot); library(segmented)


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

# Wrapper for bootstrap of the variance of the IC of Psi_n(delta)
compute_variance_for_boot <- function(data, indices, delta){
  replicate.data <- data[indices, ]
  TMLE_EY1.result <- TMLE_EY1_speedglm(replicate.data, delta)
  c(TMLE_EY1.result$Psi_n, TMLE_EY1.result$var_IC)
}

# Main directions fits: perform moving windows regressions and score them
# This includes nice plotting
main_directions_fits <- function(results_df, min_gap = 5, plotting = F){
  
  regression_df <- subset(results_df, is.finite(log_var_IC_delta) & 
                            p_value != 0 & in_range == T)
  
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
      lm.out <- lm(log_var_IC_delta ~ log_delta, in_window_subset)
      
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
                      rank(-regressions.results$avg_log_p_value)[i]^2 +
                        rank(-regressions.results$r_squared)[i]^2 +
                        500 * as.numeric(regressions.results$slope[i] < -1) +
                        0.25 * rank(-regressions.results$nb_points)[i])
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
  
  gamma_hat <- -subset(regressions.results, rank(fits.penalty, ties.method = 'first') <= 1, select = slope)[[1]] / 2
  
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
    annotate("text", x = mean(range(regression_df$log_delta)), y = mean(intercepts.range), 
             label = paste('gamma_hat = ', round(gamma_hat, 4), '\ngamma = ', gamma, sep = '')) +
    xlab(expression(log[10](delta))) + 
    ylab(expression(log[10](widehat(Delta * b[n])(delta)))) +
    ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2]),")") ==
                         group("(",list(lambda, alpha0, beta0, beta1, beta2),")"),
                       list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2)))
  
  if(plotting) print(fits.plot)
  
  list(fits.plot = fits.plot, gamma_hat = gamma_hat)
}

# Fit a broken line to the in-range area of the curve
fit_broken_line <- function(results_df, nb_breakpoints = 1, plotting = F){
  regression_df <- subset(results_df, in_range == T)
  
  initial_breakpoints <- quantile(regression_df$log_delta, (1:nb_breakpoints) / (nb_breakpoints + 1))
  
  lm.out <- lm(log_var_IC_delta ~ log_delta, regression_df)
  segmented.out <- segmented.lm(lm.out, seg.Z=~log_delta, psi = initial_breakpoints)
  broken_line_df <- as.data.frame(cbind(intercept = as.vector(intercept(segmented.out)$log_delta),
                                        slope = as.vector(slope(segmented.out)$log_delta[, 1])))
  
  fitted_breakpoints <- c(log(delta_min) / log(10), 
                          as.vector(segmented.out$psi[,2]), 
                          log(delta_max) / log(10))
  broken_line_points_df <- matrix(0, ncol = 9, nrow = length(fitted_breakpoints) - 1)
  colnames(broken_line_points_df) <- c('x.start', 'x.end', 'y.start', 'y.end',
                                       'avg_log_p_value', 'RSS', 'r_squared',
                                       'nb_points', 'leftmost')
  for(i in 1:(length(fitted_breakpoints) - 1)){
    broken_line_points_df[i, 'x.start'] <- fitted_breakpoints[i]
    broken_line_points_df[i, 'x.end'] <- fitted_breakpoints[i + 1]
    broken_line_points_df[i, 'y.start'] <- broken_line_df$intercept[i] + broken_line_df$slope[i] * fitted_breakpoints[i]
    broken_line_points_df[i, 'y.end'] <- broken_line_df$intercept[i] + broken_line_df$slope[i] * fitted_breakpoints[i + 1]
    segment.subset <- subset(regression_df, log_delta >= fitted_breakpoints[i] &
                               log_delta <= fitted_breakpoints[i + 1])
    lm_on_subset.out <- lm(log_var_IC_delta ~ log_delta, segment.subset)
    broken_line_points_df[i, 'RSS'] <- sum(lm_on_subset.out$residuals^2)
    broken_line_points_df[i, 'r_squared'] <- as.numeric(summary(lm_on_subset.out)$r.squared)
    broken_line_points_df[i, 'avg_log_p_value'] <- mean(segment.subset$log_p_value)
    broken_line_points_df[i, 'nb_points'] <- nrow(segment.subset)
    broken_line_points_df[i, 'leftmost'] <- i
  }
  broken_line_points_df <- as.data.frame(broken_line_points_df)
  
  var_IC.plot_with_broken_segments <- NULL
  if(plotting){
    var_IC.plot_with_broken_segments <- var_IC.plot + geom_segment(data = broken_line_points_df,
                                                                   mapping = aes(x = x.start, y = y.start,
                                                                                 xend = x.end, yend = y.end,
                                                                                 alpha = avg_log_p_value,
                                                                                 colour = r_squared),
                                                                   size = 2) +
      geom_point(data = data.frame(x = c(broken_line_points_df$x.start, broken_line_points_df$x.end),
                                   y = c(broken_line_points_df$y.start, broken_line_points_df$y.end)),
                 aes(x = x, y = y), colour = 'red', shape = 4, size = 5)
  }
  RSS <- sum(broken_line_points_df$RSS)
  nb_parameters <- (nb_breakpoints + 1) * 3
  BIC <- -log(RSS) + 0.5 * nb_parameters * log(nrow(regression_df))
  
  # Score the segments
  segments.squared_lengths <- (broken_line_points_df$x.end - broken_line_points_df$x.start)^2 +
    (broken_line_points_df$y.end - broken_line_points_df$y.start)^2
  
  scores <- rank(-broken_line_points_df$r_squared)^2 +
    rank(-broken_line_points_df$avg_log_p_value)^2 + 
    broken_line_points_df$leftmost + 
    -(broken_line_points_df$nb_points / max(broken_line_points_df$nb_points) * nrow(broken_line_df))^2 +
    -2 * segments.squared_lengths / max(segments.squared_lengths) +
    0.25 * rank(broken_line_df$slope) + 
    1e6  * as.numeric(broken_line_df$slope > 0)
  cat('Segment scores, from left to right:')
  print(scores)
  # browser()
  gamma_hat <- -broken_line_df$slope[which.min(scores)] / 2
  cat('gamma_hat = ', gamma_hat, '\n')
  
  cat('broken_line_df:\n')
  print(broken_line_df)
  
  if(plotting){
    var_IC.plot_with_broken_segments <- var_IC.plot_with_broken_segments + 
      geom_abline(intercept = as.numeric(broken_line_df$intercept[which.min(scores)]),
                  slope = as.numeric(broken_line_df$slope[which.min(scores)]),
                  colour = 'red', linetype = 'dotted')
    print(var_IC.plot_with_broken_segments)
  }
  
  cat('With ', nb_breakpoints, ', RSS = ', RSS, ', and BIC = ', BIC, '\n')
  list(var_IC.plot = var_IC.plot_with_broken_segments, 
       RSS = RSS, BIC = BIC,
       gamma_hat = gamma_hat)
}

# Set up cluster
cat(detectCores(), 'cores detected\n')
cl <- makeCluster(getOption("cl.cores", detectCores()), outfile = '')
registerDoParallel(cl)

# Simulate data
n <- 10^4.5
observed_data <- generate_data("L0_exp", lambda, alpha0, beta0, beta1, beta2, n)

# Set up tasks
deltas <- 10^seq(from = -5, to = -0.8, by = 0.05)

results <- foreach(delta=deltas, .combine = rbind,
                   .packages = c('speedglm', 'boot'), .verbose = T, .inorder = T) %dopar% {
                     
                     observed_data.df <- data.frame(L0 = observed_data$L0,
                                                    A0 = observed_data$A0,
                                                    L1 = observed_data$L1)
                     var_IC <- compute_variance_for_boot(observed_data.df, 1:n, delta)[2]
                     var_IC.bootstrap <- boot(data = observed_data.df,
                                              statistic = compute_variance_for_boot,
                                              R = 1000, sim = 'ordinary',
                                              delta = delta)$t
                     shapiro.p_value <- 0
                     try(shapiro.p_value <- shapiro.test(var_IC.bootstrap[, 1])$p.value)
                     if(is.null(var_IC.bootstrap)) shapiro.p_value <- 0
                     
                     c(delta, var_IC, shapiro.p_value)
                   }
stopCluster(cl)

row.names(results) <- NULL
colnames(results) <- c('delta', 'var_IC', 'p_value')

results_df <- as.data.frame(results)
results_df <- cbind(results_df, log_delta = log(results_df$delta) / log(10),
                    log_var_IC_delta = log(results_df$var_IC) / log(10),
                    log_p_value = log(results_df$p_value) / log(10))

# Make plot of finite differences
var_IC.plot <- ggplot(results_df, 
                      aes(x = log_delta, 
                          y = log_var_IC_delta)) + 
  geom_line() + geom_point(aes(size = log_p_value)) + 
  geom_abline(intercept = -0.5, slope = -2 * gamma, linetype = 'dotted') +
  geom_abline(intercept = -1, slope = -2 * gamma, linetype = 'dotted') +
  geom_abline(intercept = -1.5, slope = -2 * gamma, linetype = 'dotted') +
  xlab(expression(log[10](delta))) + 
  ylab(expression(log[10](sigma[n]^2*(delta)))) +
  ggtitle(substitute(group("(", list(Lambda, alpha[0], beta[0], beta[1], beta[2], 'n'),")") ==
                       group("(",list(lambda, alpha0, beta0, beta1, beta2, n),")"),
                     list(lambda = lambda, alpha0 = alpha0, beta0 = beta0, beta1 = beta1, beta2 = beta2, n = n)))

# Find gamma --------------------------------------------------------------

# Restrict the domain of study by cropping the curve
# log_var_min <- 0.9 * min(results_df$log_var_IC_delta) + 0.1 * max(results_df$log_var_IC_delta)
log_var_min <- min(results_df$log_var_IC_delta)
log_var_max <- 0.1 * min(results_df$log_var_IC_delta) + 0.9 * max(results_df$log_var_IC_delta)
delta_min <- max(results_df$delta[results_df$log_var_IC_delta >= log_var_max])
delta_max <- max(results_df$delta[results_df$log_var_IC_delta >= log_var_min])

var_IC.plot <- var_IC.plot + geom_hline(yintercept = log_var_min) +
  geom_hline(yintercept = log_var_max) +
  geom_vline(xintercept = log(delta_min) / log(10)) + 
  geom_vline(xintercept = log(delta_max) / log(10))

print(var_IC.plot)

results_df <- cbind(results_df, in_range = results_df$delta <= delta_max & results_df$delta >= delta_min)

# Fit a segmented regression on the search range
broken_lines_fits.BIC <- sapply(1:4, function(nb_breakpoints){
  BIC <- Inf
  try(BIC <- fit_broken_line(results_df, nb_breakpoints, F)$BIC)
  BIC
})
cat(which.min(broken_lines_fits.BIC), ' breakpoints selected based on BIC\n')
broken_lines.final_fit <- fit_broken_line(results_df, which.min(broken_lines_fits.BIC), T)

# Find main direction of fits ---------------------------------------------
main_directions_fits.results <- main_directions_fits(results_df, min_gap = 5)
print(main_directions_fits.results$fits.plot)
cat('gamma_hat =', main_directions_fits.results$gamma_hat)

# Recapitulate results
cat('Broken line fitting method: gamma_hat = ', broken_lines.final_fit$gamma_hat, '\n')
cat('Moving windows method: gamma_hat = ', main_directions_fits.results$gamma_hat, '\n')
cat('True gamma = ', gamma, '\n')

# Save the finite differences
save('results_df', file = paste('finite_differences_',
                                floor(n), '_', lambda, '_', alpha0, '_', beta0, '_', beta1, '_', beta2, '_',
                                as.integer(as.POSIXct( Sys.time() )),
                                '.results',
                                sep = ''))