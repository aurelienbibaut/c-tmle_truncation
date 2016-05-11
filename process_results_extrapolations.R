library(plyr)
library(ggplot2)
library(gridExtra)
setwd('~/aurelien.bibaut@gmail.com/Data_PC/PhD Berkeley/TMLE_truncation/Taylor expansion based estimator/c-tmle_truncation/')
raw_results <- read.csv('TMLE_extrapolations_intermediate_results.csv')
parameters_grid <- read.csv('parameters_grid_extrapolations.csv')

raw_results <- cbind(raw_results, raw_results[, "bias"]^2 + raw_results[, "var"])
colnames(raw_results)[5] <- "MSE"
# A data generating process (dgp) is a distribution + a sample size. There will be one plot per data generating process
data_generating_process <- unique(parameters_grid[, c("type", "positivity_parameter", "alpha0", "beta0", "beta1", "beta2", "n")])

# Formatted_results matrix
formatted_results <- vector()
for(i in 1:nrow(data_generating_process)){
  dgp_job_indices <- as.numeric(rownames(subset(parameters_grid, type == data_generating_process[i, ]$type
                                & positivity_parameter == data_generating_process[i, ]$positivity_parameter
                                & alpha0 == data_generating_process[i, ]$alpha0
                                & beta0 == data_generating_process[i, ]$beta0
                                & beta1 == data_generating_process[i, ]$beta1
                                & beta2 == data_generating_process[i, ]$beta2
                                & n == data_generating_process[i, ]$n)))
  dgp_results <- subset(raw_results, parameters_tuple_id %in% dgp_job_indices)
  dgp_results <- cbind(dgp_results, order = parameters_grid[dgp_results$parameters_tuple_id, ]$order,
                       delta0 = parameters_grid[dgp_results$parameters_tuple_id, ]$delta0)
  dgp_results <- cbind(dgp_results, estimatorXorder = apply(dgp_results, 1, function(x) paste(c(x["estimator"], x["order"]), collapse = "")))
  
  formatted_results <- rbind(formatted_results,
                             cbind(dgp_id = i, dgp_results))
}

MSE_plots <- list()
for(i in 1:nrow(data_generating_process)){
  title <- paste(c(as.character(data_generating_process[i, ]$type), 
                   ", positivity param", data_generating_process[i, ]$positivity_parameter,
                   "\n and n = ", data_generating_process[i, ]$n), collapse = "")
  MSE_plots[[i]] <- ggplot(subset(formatted_results, dgp_id == i), aes(x = delta0, y = MSE, colour = estimatorXorder)) + 
    geom_line() + geom_point(data = subset(formatted_results, dgp_id == i & 
                                             estimatorXorder %in% c("TMLE.extr0", "TMLE.extr.bis0")), 
                               aes(x = delta0, y = MSE, colour = estimatorXorder)) + 
    geom_hline(yintercept = min(subset(formatted_results, dgp_id == i & 
                                     estimatorXorder %in% c("TMLE.extr0", "TMLE.extr.bis0"), select = MSE))) +
    ggtitle(title) + 
    coord_cartesian(ylim=c(0.9 * min(subset(formatted_results, dgp_id == i, select = MSE)), 
                           max(subset(formatted_results, dgp_id == i & estimatorXorder == "TMLE.extr0", select = MSE))))
  #try(print(MSE_plots[[i]]))
}

for(current_n in unique(parameters_grid$n)){
  id_dpg_with_current_n <- which(data_generating_process[, ]$n == current_n)
  grid.arrange(MSE_plots[[id_dpg_with_current_n[1]]], MSE_plots[[id_dpg_with_current_n[2]]], 
               MSE_plots[[id_dpg_with_current_n[3]]], MSE_plots[[id_dpg_with_current_n[4]]],
               nrow = 2, ncol = 2)
}