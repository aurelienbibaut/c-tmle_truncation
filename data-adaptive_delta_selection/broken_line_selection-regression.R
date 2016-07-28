library(SuperLearner)
library(dummies)

broken_lines.results1 <- read.csv('broken_lines.results1.csv')
broken_lines.results2 <- read.csv('broken_lines.results2.csv')
broken_lines.results3 <- read.csv('broken_lines.results3.csv')
broken_lines.results4 <- read.csv('broken_lines.results4.csv')
broken_lines.results2$dataset_id <- broken_lines.results2$dataset_id + max(broken_lines.results1$dataset_id)
broken_lines.results3$dataset_id <- broken_lines.results3$dataset_id + max(broken_lines.results2$dataset_id)
broken_lines.results4$dataset_id <- broken_lines.results4$dataset_id + max(broken_lines.results3$dataset_id)
broken_lines.results <- rbind(broken_lines.results1, broken_lines.results2, broken_lines.results3, broken_lines.results4)

n_init <- nrow(broken_lines.results)
n_datasets <- n_init / 5
broken_lines.results <- cbind(broken_lines.results,
                              id = kronecker(1:n_datasets, rep(1, 5)))
# Filter out incomplete cases
incomplete_cases.dataset_ids <- unique(broken_lines.results$id[!complete.cases(broken_lines.results)])
complete_cases.dataset_ids <- which(! (1:n_datasets %in% incomplete_cases.dataset_ids))
broken_lines.results.cc <- subset(broken_lines.results, id %in% complete_cases.dataset_ids)
n_cc <- nrow(broken_lines.results.cc)

# Add a new feature: slopes ordering on the broken line
broken_lines.results.cc <- cbind(broken_lines.results.cc, slopes.ordering = NA)
for(current_dataset_id in broken_lines.results.cc$dataset_id){
  current_1bp_broken_line <- subset(broken_lines.results.cc, dataset_id == current_dataset_id & n_breakpoints == 1)
  ordering_1bp <- paste(current_1bp_broken_line$leftmost[order(-current_1bp_broken_line$slope)], collapse = '')
  broken_lines.results.cc[row.names(current_1bp_broken_line), 'slopes.ordering'] <- ordering_1bp
  
  current_2bp_broken_line <- subset(broken_lines.results.cc, dataset_id == current_dataset_id & n_breakpoints == 2)
  ordering_2bp <- paste(current_2bp_broken_line$leftmost[order(-current_2bp_broken_line$slope)], collapse = '')
  broken_lines.results.cc[row.names(current_2bp_broken_line), 'slopes.ordering'] <- ordering_2bp
}
broken_lines.results.cc <- dummy.data.frame(broken_lines.results.cc, names = "slopes.ordering")

# Define predictors
# predictors <- c('relative_linearity',
#   'linearity',
#   'relative_segment_squared_length',
#   'avg_log_p_value',
#   'nb_points',
#   'leftmost', 
#   'n', 'n_breakpoints',
#   colnames(broken_lines.results.cc)[which.dummy(broken_lines.results.cc)])
predictors <- c('relative_linearity',
                'linearity',
                'relative_segment_squared_length',
                'nb_points',
                'leftmost', 
                'n', 'n_breakpoints',
                colnames(broken_lines.results.cc)[which.dummy(broken_lines.results.cc)])


# Define test set and training set
training_set.dataset_ids <- sample(complete_cases.dataset_ids, floor(length(complete_cases.dataset_ids) * 4 / 5), F)
in_test_set.dataset_ids <- ! (complete_cases.dataset_ids %in% training_set.dataset_ids)
test_set.dataset_ids <- complete_cases.dataset_ids[in_test_set.dataset_ids]

training_set <- subset(broken_lines.results.cc, id %in% training_set.dataset_ids)
test_set <- subset(broken_lines.results.cc, id %in% test_set.dataset_ids)

# GLM
glm_formula <- as.formula(paste(c('is_best ~', paste(predictors, collapse = '+')), collapse = ''))

glm_fit <- glm(glm_formula, 
        data = subset(broken_lines.results.cc, id %in% training_set.dataset_ids),
        family = binomial)
# Predict on test set for GLM
test_set.prediction <- predict(glm_fit, newdata = subset(broken_lines.results.cc, id %in% test_set.dataset_ids), type = "response")


# Use SuperLearner
SL.library <- c("SL.glm", "SL.randomForest", "SL.gam",
                "SL.polymars", "SL.mean", "SL.nnet", "SL.bayesglm", "SL.step.interaction")
super_learner.result <- SuperLearner(Y = training_set$is_best, X = training_set[, predictors],
                                     newX = test_set[, predictors],
                                     family = binomial(), SL.library = SL.library)


#Subset of test set dataset ids, for which n > n0
n0 <- 1e4
test_set.dataset_ids.subset <- unique(subset(broken_lines.results.cc, dataset_id %in% test_set.dataset_ids & n > n0, select = dataset_id)[[1]])

glm.classification_error_rate <- 0
SL.classification_error_rate <- 0

MSE_best_choices <- rep(0, 5) # MSEs of the first best choice, of the 2 first best choices,..., of the 5 first best choices
MSE_glm <- 0
MSE_SL <- 0

truth_and_predictions.matrix <- vector()
for(dataset_id in test_set.dataset_ids.subset){
  cat('Dataset id ', dataset_id, ', rows\' numbers:\n')
  row_names <- (dataset_id - 1) * 5 + 1:5
  print(row_names)
  row_numbers_in_new_df <- which(as.numeric(row.names(broken_lines.results.cc)) == row_names)
  
  
  # print(broken_lines.results.cc$is_best[row_numbers])
  true_best <- which(broken_lines.results.cc$is_best[row_numbers_in_new_df] == 1)[1]
  cat("Is best: ", true_best, '\n')
  # cat("Predicted probabilites:\n")
  # print(test_set.prediction[as.character(row_numbers)])
  glm.predicted_best <- which.max(as.vector(test_set.prediction[as.character(row_names)]))
  SL.predicted_best <- which.max(super_learner.result$SL.predict[as.numeric(row.names(super_learner.result$SL.predict)) == row_names])
  cat("GLM Predicted_ best: ", glm.predicted_best, '\n')
  cat("SL Predicted_ best: ", SL.predicted_best, '\n')
  
  
  # Compute MSEs
  squared_residuals <- (-broken_lines.results.cc[row_numbers_in_new_df, 'slope'] - 
                          broken_lines.results.cc[row_numbers_in_new_df[1], 'true_gamma'])^2
  for(i in 1:5){
    MSE_best_choices[i] <- MSE_best_choices[i] + 1 / (length(test_set.dataset_ids.subset) * i) * sum(sort(squared_residuals)[1:i])
  }
  MSE_glm <- MSE_glm + 1 / length(test_set.dataset_ids.subset) * squared_residuals[glm.predicted_best]
  MSE_SL <- MSE_glm + 1 / length(test_set.dataset_ids.subset) * squared_residuals[SL.predicted_best]
  
  cat("Dataset id : ", dataset_id, ' and MSEs :\n')
  print(MSE_best_choices)
  
  # Add prediction and truth to results dataset
  truth_and_predictions.matrix <- rbind(truth_and_predictions.matrix,
                                        c(dataset_id, true_best, glm.predicted_best, SL.predicted_best))
  
  glm.classification_error_rate <- glm.classification_error_rate + as.numeric(true_best != glm.predicted_best) / length(test_set.dataset_ids.subset)
  SL.classification_error_rate <- SL.classification_error_rate + as.numeric(true_best != SL.predicted_best) / length(test_set.dataset_ids.subset)
}

truth_and_predictions.df <- as.data.frame(truth_and_predictions.matrix)

cat("Classification error rate:\n
    -GLM: ", glm.classification_error_rate, "\n
    -SL: ", SL.classification_error_rate, "\n")

cat('MSEs of first k best segments for each broken line, k from 1 to 5:\n')
print(MSE_best_choices)
cat('MSE of slope chosen by glm:', MSE_glm, '\n')
cat('MSE of slope chosen by SL:', MSE_SL, '\n')