min_gap <- 5

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
                             y.start = regressions.results$lower_bound * regressions.results$slope + regressions.results$intercept,
                             y.end = regressions.results$upper_bound * regressions.results$slope + regressions.results$intercept,
                             fits.penalty = fits.penalty)
# Lines with true slope
intercepts.range <- range(c(regressions.results$y.end, regressions.results$y.start))
true_lines.df <- data.frame(intercept = seq(from = intercepts.range[1] - 2, 
                                            to = intercepts.range[2] + 2, 
                                            length = 15) + beta * min(regressions.results$lower_bound),
                            slope = rep(-beta, 15))

print(fin_diffs.plot + geom_segment(data = regressions.results, aes(x = lower_bound, y = y.start,
                                                                    xend = upper_bound, yend = y.end), colour = 'black'))
print(fin_diffs.plot)

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
