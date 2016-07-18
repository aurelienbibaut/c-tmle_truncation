source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

library(ggplot2); library(gridExtra); library(grid)
library(foreach); library(doParallel)
library(robustbase); library(speedglm)
library(boot); library(segmented)

# Define inference functions ----------------------------------------------

# Wrapper for bootstrap of the variance of the IC of Psi_n(delta)
compute_variance_for_boot <- function(data, indices, delta){
  replicate.data <- data[indices, ]
  TMLE_EY1.result <- TMLE_EY1_speedglm(replicate.data, delta)
  c(TMLE_EY1.result$Psi_n, TMLE_EY1.result$var_IC)
}

# Main directions fits: perform moving windows regressions and score them
# This includes nice plotting

# Fit a broken line to the in-range area of the curve
fit_broken_line <- function(results_df, nb_breakpoints = 2, delta_min, delta_max, plotting = F){
  regression_df <- subset(results_df, in_range == T)
  
  initial_breakpoints <- quantile(regression_df$log_delta, (1:nb_breakpoints) / (nb_breakpoints + 1))
  
  lm.out <- lm(log_var_IC_delta ~ log_delta, regression_df)
  segmented.fit.sucess <- F; nb_its <- 0
  while(!segmented.fit.sucess){
    try.out <- try(segmented.out <- segmented.lm(lm.out, seg.Z=~log_delta, psi = initial_breakpoints))
    if(paste(class(try.out), collapse = '') != "try-error") segmented.fit.sucess <- T
    nb_its <- nb_its + 1
    if(nb_its > 15) stop('Could not fit broken line')
  }
  broken_line_df <- as.data.frame(cbind(intercept = as.vector(intercept(segmented.out)$log_delta),
                                        slope = as.vector(slope(segmented.out)$log_delta[, 1])))
  
  fitted_breakpoints <- c(log(delta_min) / log(10), 
                          as.vector(segmented.out$psi[,2]), 
                          log(delta_max) / log(10))
  broken_line_points_df <- matrix(0, ncol = 13, nrow = length(fitted_breakpoints) - 1)
  colnames(broken_line_points_df) <- c('x.start', 'x.end', 'y.start', 'y.end',
                                       'avg_log_p_value', 'RSS', 'r_squared',
                                       'nb_points', 'leftmost', 'segment.squared_length',
                                       'slope', 'loss', 'is_best')
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
    broken_line_points_df[i, 'slope'] <- broken_line_df$slope[i]
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
  BIC <- -log(RSS) + 0.25 * nb_parameters * log(nrow(regression_df))
  
  # Figure out which segment is best
  broken_line_points_df$loss <- sapply(broken_line_points_df$slope, function(x) (abs(x) - beta)^2 + 1e6 * as.numeric(x < -1 | x > 0))
  broken_line_points_df$is_best[which.min(broken_line_points_df$loss)] <- 1
  
  # Score the segments
  segments.squared_lengths <- (broken_line_points_df$x.end - broken_line_points_df$x.start)^2 +
    (broken_line_points_df$y.end - broken_line_points_df$y.start)^2
  
  broken_line_points_df$segment.squared_length <- segments.squared_lengths
  broken_line_points_df$relative_segment_squared_length <- segments.squared_lengths / max(segments.squared_lengths)
  broken_line_points_df$linearity <- -log(1 - broken_line_points_df$r_squared) / log(10)
  broken_line_points_df$relative_linearity <- broken_line_points_df$linearity / max(broken_line_points_df$linearity)
  
  scores <- 2*rank(-broken_line_points_df$r_squared)^2 +
    -0.5 * log(1 - broken_line_points_df$r_squared) / log(10) * (broken_line_points_df$nb_points > 6)
  rank(-broken_line_points_df$avg_log_p_value)^2 + 
    -(broken_line_points_df$leftmost == 2) + 
    -(broken_line_points_df$nb_points / max(broken_line_points_df$nb_points) * nrow(broken_line_df))^2 +
    -2 * segments.squared_lengths / max(segments.squared_lengths) +
    0.25 * rank(broken_line_df$slope) + 
    1e6  * as.numeric(broken_line_df$slope > 0 | broken_line_df$slope < -1)
  
  cat('Segment scores, from left to right:')
  
  print(broken_line_points_df)
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
       gamma_hat = gamma_hat,
       broken_line_points_df = broken_line_points_df)
}

generate_data_and_gamma_broken_line <- function(type, lambda, alpha0, beta0, beta1, beta2, n, gamma){
  # Set up cluster
  cat(detectCores(), 'cores detected\n')
  cl <- makeCluster(getOption("cl.cores", detectCores()), outfile = '')
  registerDoParallel(cl)
  
  # Set up tasks
  deltas <- 10^seq(from = -5, to = -0.8, by = 0.05)
  
  # Simulate data
  observed_data <- generate_data(type, lambda, alpha0, beta0, beta1, beta2, n)
  
  results <- foreach(delta=deltas, .combine = rbind,
                     .packages = c('speedglm', 'boot'), .verbose = T,
                     .export = c('beta', 'gamma', 'kappa', 'n',
                                 'lambda', 'alpha0', 'beta0', 'beta1', 'beta2',
                                 'TMLE_EY1_speedglm', 'expit', 'logit', 'g_to_g_delta',
                                 'generate_data', 'compute_variance_for_boot',
                                 'observed_data'),
                     .inorder = T) %dopar% {
                       
                       observed_data.df <- data.frame(L0 = observed_data$L0,
                                                      A0 = observed_data$A0,
                                                      L1 = observed_data$L1)
                       
                       var_IC <- compute_variance_for_boot(observed_data.df, 1:n, delta)[2]
                       var_IC.bootstrap <- boot(data = observed_data.df,
                                                statistic = compute_variance_for_boot,
                                                R = 100, sim = 'ordinary',
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
  
  # Is there a discontinuity (probably due to numerical issues)
  nb_points <- length(results_df$log_var_IC_delta)
  first_diff <- abs(results_df$log_var_IC_delta[2:nb_points] - results_df$log_var_IC_delta[1:(nb_points - 1)])
  jumps <- abs(first_diff) > 20 * mean(abs(first_diff))
  if(any(jumps)){
    delta_first_jump <- max(results_df$delta[which(jumps) + 1])
    delta_min <- max(delta_first_jump, delta_min)
    log_var_max <- results_df$log_var_IC_delta[which(results_df$delta == delta_min)]
  }
  var_IC.plot <- var_IC.plot + geom_hline(yintercept = log_var_min) +
    geom_hline(yintercept = log_var_max) +
    geom_vline(xintercept = log(delta_min) / log(10)) + 
    geom_vline(xintercept = log(delta_max) / log(10))
  
  print(var_IC.plot)
  
  results_df <- cbind(results_df, in_range = results_df$delta <= delta_max & results_df$delta >= delta_min)
  
  broken_lines.final_fit <- fit_broken_line(results_df, 2, delta_min, delta_max)
}
