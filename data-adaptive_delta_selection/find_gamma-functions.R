source('../utilities.R')
source('../true_target_parameters_derivatives_and_ICs.R')
source('../generate_data.R')
source('../TMLE_extrapolation_functions.R')

if(running_environment != 'SAVIO2'){
  library(speedglm)
  library(boot)
  library(segmented)
}else{
  library(speedglm, lib.loc = "/global/home/users/afbibaut/Rlibs")
  library(boot, lib.loc = "/global/home/users/afbibaut/Rlibs")
  library(segmented, lib.loc = "/global/home/users/afbibaut/Rlibs")
}

# Define inference functions ----------------------------------------------

# Wrapper for bootstrap of the variance of the IC of Psi_n(delta)
compute_variance <- function(observed_data, delta){
  TMLE_EY1.result <- TMLE_EY1_speedglm(observed_data, delta)
  c(TMLE_EY1.result$Psi_n, TMLE_EY1.result$var_IC)
}

# Main directions fits: perform moving windows regressions and score them
# This includes nice plotting

# Fit a broken line to the in-range area of the curve
fit_broken_line_gamma <- function(results_df, nb_breakpoints = 2, delta_min, delta_max, 
                                  gamma, plotting = F, var_IC.plot = NULL){
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
                          as.vector(segmented.out$psi[, 2]), 
                          log(delta_max) / log(10))
  broken_line_points_df <- matrix(0, ncol = 14, nrow = length(fitted_breakpoints) - 1)
  colnames(broken_line_points_df) <- c('x.start', 'x.end', 'y.start', 'y.end',
                                       'RSS', 'r_squared',
                                       'nb_points', 'leftmost', 'segment.squared_length',
                                       'gamma.slope', 'gamma.intercept',
                                       'exp.gamma.intercept',
                                       'loss', 'is_best')
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
    broken_line_points_df[i, 'nb_points'] <- nrow(segment.subset)
    broken_line_points_df[i, 'leftmost'] <- i
    broken_line_points_df[i, 'gamma.slope'] <- broken_line_df$slope[i]
    broken_line_points_df[i, 'gamma.intercept'] <- broken_line_df$intercept[i]
    broken_line_points_df[i, 'exp.gamma.intercept'] <- exp(broken_line_df$intercept[i])
  }
  broken_line_points_df <- as.data.frame(broken_line_points_df)
  
  var_IC.plot_with_broken_segments <- NULL
  if(plotting){
    var_IC.plot_with_broken_segments <- var_IC.plot + geom_segment(data = broken_line_points_df,
                                                                   mapping = aes(x = x.start, y = y.start,
                                                                                 xend = x.end, yend = y.end,
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
  # browser()
  broken_line_points_df$loss <- sapply(broken_line_points_df$gamma.slope, function(x) (abs(x) - 2 * gamma)^2 + 1e6 * as.numeric(x < -1 | x > 0))
  
  # Score the segments
  segments.squared_lengths <- (broken_line_points_df$x.end - broken_line_points_df$x.start)^2 +
    (broken_line_points_df$y.end - broken_line_points_df$y.start)^2
  
  broken_line_points_df$segment.squared_length <- segments.squared_lengths
  broken_line_points_df$relative_segment_squared_length <- segments.squared_lengths / max(segments.squared_lengths)
  broken_line_points_df$linearity <- -log(1 - broken_line_points_df$r_squared) / log(10)
  broken_line_points_df$relative_linearity <- broken_line_points_df$linearity / max(broken_line_points_df$linearity)
  
  # print(broken_line_points_df)
  
  if(plotting){
    var_IC.plot_with_broken_segments <- var_IC.plot_with_broken_segments + 
      geom_abline(intercept = as.numeric(broken_line_df$intercept[which.min(broken_line_points_df$loss)]),
                  slope = as.numeric(broken_line_df$slope[which.min(broken_line_points_df$loss)]),
                  colour = 'red', linetype = 'dotted')
    print(var_IC.plot_with_broken_segments)
  }
  
  cat('With ', nb_breakpoints, ', RSS = ', RSS, ', and BIC = ', BIC, '\n')
  list(var_IC.plot = var_IC.plot_with_broken_segments, 
       RSS = RSS, BIC = BIC,
       broken_line_points_df = broken_line_points_df)
}

# Compute variances
compute_variances <- function(observed_data, verbose = F){
  # Set up tasks
  deltas <- 10^seq(from = -5, to = -0.8, by = 0.05)
  
  results <- foreach(delta=deltas, .combine = rbind,
                     .packages = c('speedglm', 'boot'), .verbose = verbose,
                     .export = c('TMLE_EY1_speedglm', 'expit', 'logit', 'g_to_g_delta',
                                 'generate_data', 'compute_variance',
                                 'observed_data'),
                     .inorder = T) %dopar% {
                       
                       observed_data.df <- data.frame(L0 = observed_data$L0,
                                                      A0 = observed_data$A0,
                                                      L1 = observed_data$L1)
                       
                       var_IC <- compute_variance(observed_data.df, delta)[2]
                       
                       c(delta, var_IC)
                     }
  
  row.names(results) <- NULL
  colnames(results) <- c('delta', 'var_IC')
  
  results_df <- as.data.frame(results)
  cbind(results_df, log_delta = log(results_df$delta) / log(10),
        log_var_IC_delta = log(results_df$var_IC) / log(10))
}

# Narrow down the study domain
restrict_study_domain <- function(results_df){
  # Restrict the domain of study by cropping the curve
  # log_var_min <- 0.9 * min(results_df$log_var_IC_delta) + 0.1 * max(results_df$log_var_IC_delta)
  
  log_var_min <- min(results_df$log_var_IC_delta)
  log_var_max <- 0.1 * min(results_df$log_var_IC_delta) + 0.9 * max(results_df$log_var_IC_delta)
  delta_min <- max(results_df$delta[results_df$log_var_IC_delta >= log_var_max])
  delta_max <- max(results_df$delta[results_df$log_var_IC_delta >= log_var_min])
  
  # Is there a discontinuity ? (probably due to numerical issues)
  nb_points <- length(results_df$log_var_IC_delta)
  first_diff <- abs(results_df$log_var_IC_delta[2:nb_points] - results_df$log_var_IC_delta[1:(nb_points - 1)])
  jumps <- abs(first_diff) > 15 * mean(abs(first_diff))
  if(any(jumps)){
    delta_first_jump <- max(results_df$delta[which(jumps) + 1])
    if(sum(results_df$delta <= delta_max & results_df$delta >= delta_min) >= 10){
      delta_min <- max(delta_first_jump, delta_min)
      log_var_max <- results_df$log_var_IC_delta[which(results_df$delta == delta_min)]
    }
  }
  
  results_df <- cbind(results_df, in_range = results_df$delta <= delta_max & results_df$delta >= delta_min)
  list(var_IC_df = subset(results_df, in_range == TRUE), delta_min = delta_min, delta_max = delta_max)
}

# Plot log_var_IC against log(delta)
plot_log_var_IC <-  function(results_df){
  ggplot(results_df, aes(x = log_delta, y = log_var_IC_delta)) +  geom_line() + geom_point()
}

# Generate gamma datapoint
extract_gamma_features <- function(var_IC_df, gamma, plotting = F, var_IC.plot = NULL){
  restrict_study_domain.results <- restrict_study_domain(var_IC_df)
  var_IC_df <- restrict_study_domain.results$var_IC_df
  delta_min <- restrict_study_domain.results$delta_min
  delta_max <- restrict_study_domain.results$delta_max
  nrow.var_IC_df <- nrow(var_IC_df)
  
  broken_lines.1bp_fit <- cbind(fit_broken_line_gamma(var_IC_df, 1, delta_min, delta_max, gamma, 
                                                      plotting = plotting, var_IC.plot = var_IC.plot)$broken_line_points_df,
                                nrow.var_IC_df = nrow.var_IC_df, n_breakpoints = 1)
  broken_lines.1bp_fit <- cbind(broken_lines.1bp_fit,
                                slopes.ordering = paste(order(broken_lines.1bp_fit$gamma.slope), sep = '', collapse = ''))
  broken_lines.2bp_fit <- cbind(fit_broken_line_gamma(var_IC_df, 2, delta_min, delta_max, gamma, 
                                                      plotting = plotting, var_IC.plot = var_IC.plot)$broken_line_points_df,
                                nrow.var_IC_df = nrow.var_IC_df, n_breakpoints = 2)
  broken_lines.2bp_fit <- cbind(broken_lines.2bp_fit,
                                slopes.ordering = paste(order(broken_lines.2bp_fit$gamma.slope), sep = '', collapse = ''))
  
  broken_lines.results <- rbind(broken_lines.1bp_fit, broken_lines.2bp_fit)
  broken_lines.results$is_best <- rep(0, nrow(broken_lines.results))
  broken_lines.results$is_best[which.min(broken_lines.results$loss)] <- 1
  broken_lines.results$true_gamma <- rep(gamma, nrow(broken_lines.results))
  
  broken_lines.results <- cbind(id = 1, timevar = 1:5, as.data.frame(broken_lines.results))
  
  reshape(broken_lines.results, v.names = setdiff(colnames(broken_lines.results),
                                                  c("x.start", "y.start", "x.end", "y.end", "timevar", "id")),
          idvar = "id", timevar = "timevar", drop = c("x.start", "y.start", "x.end", "y.end"),
          direction = "wide")
}
