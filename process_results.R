library(plyr)
setwd('~/aurelien.bibaut@gmail.com/Data_PC/PhD Berkeley/TMLE_truncation/Taylor expansion based estimator/c-tmle_truncation/')
raw_results <- read.csv('C-TMLE_multi_orders_intermediate_results.csv')
parameters_grid <- read.csv('parameters_grid.csv')

EYd <- raw_results[1, "EYd"]

# Filter out absurd results like the ones for which there is a EYd out of [0,1]. Report them
absurd_results_ids <- which(raw_results[, "EYd"] < 0 | raw_results[, "EYd"] > 1)
cat("The following raw results had an absurd value for EYd (i.e. not in [0,1]):\n")
print(absurd_results_ids)
raw_results <- raw_results[-absurd_results_ids, ]


# Format results matrix: each row is the result of one estimation task (pair of a parameters tuple id and an estimator)
formatted_results <- rbind(cbind(raw_results[, -(4:6)], Estimate = raw_results[, "Utgtd.untr"], Estimator = "TMLE.EYd"), 
                           cbind(raw_results[, -(4:6)], Estimate = raw_results[, "Utgtd.extr"], Estimator = "C.TMLE-bis"), 
                           cbind(raw_results[, -(4:6)], Estimate = raw_results[, "C.TMLE"], Estimator = "C.TMLE"))

# An estimation task is the combination of a paramters tuple id and of an estimator
estimation_tasks <- expand.grid(parameters_tuple_id = 1:nrow(parameters_grid), Estimator = unique(formatted_results[, "Estimator"]))

bias_var_MSE_matrix <- matrix(nrow = nrow(estimation_tasks), ncol = 6)
colnames(bias_var_MSE_matrix) <- c("parameters_tuple_id", "Estimator", "bias", "var", "mse", "n")

for(i in 1:nrow(estimation_tasks)){
  task_results <- subset(formatted_results, select = c(Estimate, EYd), 
                         parameters_tuple_id == estimation_tasks[i,]$parameters_tuple_id &
                           Estimator == estimation_tasks[i,]$Estimator)
  bias_var_MSE_matrix[i, "n"] <- parameters_grid[estimation_tasks[i,]$parameters_tuple_id,]$n
  
  bias_var_MSE_matrix[i, "bias"] <- abs(mean(task_results$Estimate - task_results$EYd))
  
  bias_var_MSE_matrix[i, "var"] <- var(task_results$Estimate)
  
  if(!is.na(bias_var_MSE_matrix[i, "bias"]) & !is.nan(bias_var_MSE_matrix[i, "bias"])) 
    bias_var_MSE_matrix[i, "mse"] <- bias_var_MSE_matrix[i, "bias"]^2 + bias_var_MSE_matrix[i, "var"]
}

bias_var_MSE_matrix[, 1:2] <- as.matrix(estimation_tasks)
bias_var_MSE_df <- as.data.frame(bias_var_MSE_matrix)
bias_var_MSE_df <- transform(bias_var_MSE_matrix, parameters_tuple_id = as.numeric(as.character(parameters_tuple_id)),
                             bias = as.numeric(as.character(bias)),
                             var = as.numeric(as.character(var)),
                             mse = as.numeric(as.character(mse)),
                             Estimator = as.character(Estimator),
                             n = as.numeric(as.character(n)))

# TMLE indices counts for each n
C_TMLE_indices_counts <- list()
for(i in 1:nrow(parameters_grid))
#   C_TMLE_indices_counts[[parameters_tuple_id]] <- count(raw_results[raw_results[, "parameters_tuple_id"] == parameters_tuple_id, ], 
#                                                         c("order", "delta0"))
  C_TMLE_indices_counts[[i]] <- count(subset(raw_results, parameters_tuple_id == i), 
                                                        c("order", "delta0"))

print(C_TMLE_indices_counts)

# Plots
library(ggplot2); library(gridExtra)
# Plot MSE as a function of n. One plot per target parameter.
MSE_plots <- list(); delta0s_counts <- list(); delta0s_heatmaps <- list(); delta0s_frequencies <- list()
tp_parameters <- unique(parameters_grid[, c(1:6, 8)],)
for(i in 1:nrow(tp_parameters)){
  pg_row_numbers <- as.numeric(row.names(match_df(parameters_grid[, c(1:6, 8)], tp_parameters[i,])))
  
  title <- paste(c("L0 type", as.character(tp_parameters[i, "type"]), 
                   ", positivity parameter = ", tp_parameters[i, "positivity_parameter"],
                   "\norder = ", tp_parameters[i, "orders_set_id"] - 1),
                   collapse = " ")
  
  data_for_current_plot = na.omit(subset(bias_var_MSE_df, parameters_tuple_id %in% pg_row_numbers & n < 1000))
  MSE_plots[[i]] <- ggplot(data = data_for_current_plot,
                                                 aes(x = n, y = mse, colour = Estimator)) + geom_point() + geom_line() +
    ggtitle(title) + ylim(max(min(data_for_current_plot$mse), 0), 
                          min(max(data_for_current_plot$mse[data_for_current_plot$mse <= 1]), 1))
  
  data_for_delta0s_heatmap <- na.omit(subset(raw_results, parameters_tuple_id %in% pg_row_numbers))
  data_for_delta0s_heatmap = cbind(data_for_delta0s_heatmap,
                                    n = parameters_grid[as.numeric(as.character(data_for_delta0s_heatmap[, "parameters_tuple_id"])), "n"])

  delta0s_counts[[i]] <- count(data_for_delta0s_heatmap, c("n", "delta0"))
  delta0s_frequencies[[i]] <- delta0s_counts[[i]]
  for(j in 1:nrow(delta0s_counts[[i]]))
    delta0s_frequencies[[i]][j, "freq"] <- delta0s_counts[[i]][j, "freq"] / sum(subset(delta0s_counts[[i]], 
                                                n == delta0s_counts[[i]][j, "n"],
                                                select = freq))
  delta0s_heatmaps[[i]] <- ggplot(data = as.data.frame(delta0s_frequencies[[i]]), aes(x = factor(delta0), y = factor(n))) + 
    geom_tile(aes(fill = freq)) + ggtitle(title)
}

# grid.arrange(MSE_plots[[1]], MSE_plots[[2]], MSE_plots[[3]], nrow = 2)
for(i in 0:3)
  grid.arrange(MSE_plots[[4*i + 1]], MSE_plots[[4*i + 2]], MSE_plots[[4*i + 3]], MSE_plots[[4*i + 4]], nrow = 2)

for(i in 0:3)
  grid.arrange(delta0s_heatmaps[[4*i + 1]], delta0s_heatmaps[[4*i + 2]], 
               delta0s_heatmaps[[4*i + 3]], delta0s_heatmaps[[4*i + 4]], nrow = 2)
# Plot heatmap of frequencies of indices picked by C-TMLE
# ggplot(data = as.data.frame(C_TMLE_indices_counts[[2]]), aes(x = delta0, y = order)) + geom_tile(aes(fill = freq))
