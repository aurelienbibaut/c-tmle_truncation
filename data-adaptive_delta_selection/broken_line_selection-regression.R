library(SuperLearner)

broken_lines.results <- read.csv('broken_lines.results.new.csv')
n_init <- nrow(broken_lines.results)
n_datasets <- n_init / 5
broken_lines.results <- cbind(broken_lines.results,
                              id = kronecker(1:n_datasets, rep(1, 5)))
# Filter out incomplete cases
incomplete_cases.dataset_ids <- unique(broken_lines.results$id[!complete.cases(broken_lines.results)])
complete_cases.dataset_ids <- which(! (1:n_datasets %in% incomplete_cases.dataset_ids))
broken_lines.results.cc <- subset(broken_lines.results, id %in% complete_cases.dataset_ids)
n_cc <- nrow(broken_lines.results.cc)

# Define test set and training set
training_set.dataset_ids <- sample(complete_cases.dataset_ids, floor(length(complete_cases.dataset_ids) * 4 / 5), F)
in_test_set.dataset_ids <- ! (complete_cases.dataset_ids %in% training_set.dataset_ids)
test_set.dataset_ids <- complete_cases.dataset_ids[in_test_set.dataset_ids]

training_set <- subset(broken_lines.results.cc, id %in% training_set.dataset_ids)
test_set <- subset(broken_lines.results.cc, id %in% test_set.dataset_ids)

# Fit logistic regression model on training set
glm_fit <- glm(is_best ~ relative_linearity + linearity + relative_segment_squared_length +
                 avg_log_p_value + nb_points + leftmost + n + n_breakpoints, 
               data = subset(broken_lines.results.cc, id %in% training_set.dataset_ids),
               family = binomial)

# Predict on test set for GLM
test_set.prediction <- predict(glm_fit, newdata = subset(broken_lines.results.cc, id %in% test_set.dataset_ids), type = "response")

# Use SuperLearner
SL.library <- c("SL.glm", "SL.randomForest", "SL.gam",
                "SL.polymars", "SL.mean", "SL.nnet", "SL.bayesglm", "SL.step.interaction")
super_learner.result <- SuperLearner(Y = training_set$is_best, X = training_set[, c('relative_linearity',
                                                                                    'linearity',
                                                                                    'relative_segment_squared_length',
                                                                                    'avg_log_p_value',
                                                                                    'nb_points',
                                                                                    'leftmost', 
                                                                                    'n', 'n_breakpoints')],
                                     newX = test_set[, c('relative_linearity',
                                                         'linearity',
                                                         'relative_segment_squared_length',
                                                         'avg_log_p_value',
                                                         'nb_points',
                                                         'leftmost',
                                                         'n', 'n_breakpoints')],
                                     family = binomial(), SL.library = SL.library)


glm.classification_error_rate <- 0
SL.classification_error_rate <- 0
truth_and_predictions.matrix <- vector()
for(dataset_id in test_set.dataset_ids){
  cat('Dataset id ', dataset_id, ', rows\' numbers:\n')
  row_numbers <- (dataset_id - 1) * 5 + 1:5
  print(row_numbers)
  # print(broken_lines.results.cc$is_best[row_numbers])
  true_best <- which(broken_lines.results.cc$is_best[as.numeric(row.names(broken_lines.results.cc)) == row_numbers] == 1)[1]
  cat("Is best: ", true_best, '\n')
  # cat("Predicted probabilites:\n")
  # print(test_set.prediction[as.character(row_numbers)])
  glm.predicted_best <- which.max(as.vector(test_set.prediction[as.character(row_numbers)]))
  SL.predicted_best <- which.max(super_learner.result$SL.predict[as.numeric(row.names(super_learner.result$SL.predict)) == row_numbers])
  cat("GLM Predicted_ best: ", SL.predicted_best, '\n')
  cat("SL Predicted_ best: ", SL.predicted_best, '\n')
  
  # Add prediction and truth to results dataset
  truth_and_predictions.matrix <- rbind(truth_and_predictions.matrix,
                          c(dataset_id, true_best, glm.predicted_best, SL.predicted_best))
  
  glm.classification_error_rate <- glm.classification_error_rate + as.numeric(true_best == glm.predicted_best) / length(test_set.dataset_ids)
  SL.classification_error_rate <- SL.classification_error_rate + as.numeric(true_best == SL.predicted_best) / length(test_set.dataset_ids)
}

truth_and_predictions.df <- as.data.frame(truth_and_predictions.matrix)

cat("Classification error rate:\n
    -GLM: ", glm.classification_error_rate, "\n
    -SL: ", SL.classification_error_rate, "\n")
