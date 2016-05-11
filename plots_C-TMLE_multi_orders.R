setwd('~/aurelien.bibaut@gmail.com/Data_PC/PhD Berkeley/TMLE_truncation/Taylor expansion based estimator/c-tmle_truncation/')
results <- read.csv('C-TMLE_multi_orders_intermediate_results.csv')

Psi_d0 <- 0.2207294

# Compute the MSEs for each parameter tuple id
MSEs <- matrix(NA, nrow = nrow(parameters_grid), ncol = 4)
for(i in 1:nrow(parameters_grid)){
  MSE_utgtd <- mean((results[results[,"parameters_tuple_id"] == i, "Utgtd"] - Psi_d0)^2)
  MSE_C_TMLE <- mean((results[results[,"parameters_tuple_id"] == i, "C.TMLE"] - Psi_d0)^2)
  MSEs[i, ] <- c(i, parameters_grid[i,]$n, MSE_utgtd, MSE_C_TMLE)
}
colnames(MSEs) <- c("parameter_tuple_id", "n", "MSE_utgtd", "MSE_C-TMLE")

# Save results and plots
# pdf("plots.pdf")
plot(MSEs[, "n"],  MSEs[, "n"]* MSEs[,"MSE_utgtd"])
plot(MSEs[, "n"],  MSEs[, "n"]* MSEs[,"MSE_C-TMLE"])
# dev.off()

# save(full_results_matrix, MSEs, parameters_grid, file="C-TMLE_truncation-results")