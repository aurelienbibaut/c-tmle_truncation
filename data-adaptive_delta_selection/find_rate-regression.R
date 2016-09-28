library(h2o)
library(h2oEnsemble)
useEnsemble <- T

features_and_rates <- read.csv('merged_rate_inference.features.results.csv')[, -1]

if(!("true_rate" %in% colnames(features_and_rates))){
  features_and_rates <- cbind(features_and_rates, 
                              true_rate = 1 / (2 * (features_and_rates$true_gamma.1 + 1 - features_and_rates$true_beta)))
}
superfluous_columns <- colnames(features_and_rates)[grep(pattern = "Delta.delta_rate|leftmost|n_breakpoints|segment_id", colnames(features_and_rates))]
features_and_rates.cc <- features_and_rates[!is.na(features_and_rates$true_rate), setdiff(colnames(features_and_rates), superfluous_columns)]

# Define test set and training set
complete_cases.dataset_ids <- features_and_rates.cc$dataset_id
training_set.dataset_ids <- sample(complete_cases.dataset_ids, floor(length(complete_cases.dataset_ids) * 4 / 5), F)
in_test_set.dataset_ids <- ! (complete_cases.dataset_ids %in% training_set.dataset_ids)
test_set.dataset_ids <- complete_cases.dataset_ids[in_test_set.dataset_ids]

training_set <- subset(features_and_rates.cc, dataset_id %in% training_set.dataset_ids)
test_set <- subset(features_and_rates.cc, dataset_id %in% test_set.dataset_ids)

localH2O <- h2o.init(nthreads = -1)
training_set.h2o <- as.h2o(training_set)
test_set.h2o <- as.h2o(test_set)
full_dataset.h2o <- as.h2o(features_and_rates.cc)

# Train neural network on training set
cat("About to run the rate regression on training set\n")
# h2o.rate_regression_fit <- h2o.deeplearning(x = setdiff(colnames(features_and_rates.cc),
#                                                         c('dataset_id', 'id', 'true_rate',
#                                                           colnames(features_and_rates.cc)[grep(pattern = "loss|true|is_best|intercept", colnames(features_and_rates.cc))])),
#                                             y = 'true_rate', training_frame = training_set.h2o)

if(useEnsemble){
  learner <- c("h2o.glm.wrapper", "h2o.randomForest.wrapper", 
               "h2o.gbm.wrapper", "h2o.deeplearning.wrapper")
  metalearner <- "h2o.glm.wrapper"
  
  h2o.rate_regression_fit <- h2o.ensemble(x = setdiff(colnames(features_and_rates.cc),
                                                      c('dataset_id', 'id', 'true_rate',
                                                        colnames(features_and_rates.cc)[grep(pattern = "loss|true|is_best", colnames(features_and_rates.cc))])),
                                          y = 'true_rate', training_frame = training_set.h2o,
                                          family = "gaussian", 
                                          learner = learner, 
                                          metalearner = metalearner,
                                          cvControl = list(V = 5))
}else{
  h2o.rate_regression_fit <- h2o.deeplearning(x = setdiff(colnames(features_and_rates.cc),
                                                          c('dataset_id', 'id', 'true_rate',
                                                            colnames(features_and_rates.cc)[grep(pattern = "loss|true|is_best|intercept", colnames(features_and_rates.cc))])),
                                              y = 'true_rate', training_frame = training_set.h2o)
}
# cat("About to run the constant regression on training set\n")
# h2o.constant_regression_fit <- h2o.deeplearning(x = setdiff(colnames(features_and_rates.cc),
#                                                             c('dataset_id', 'id', 'true_rate',
#                                                               colnames(features_and_rates.cc)[grep(pattern = "loss|true|is_best", colnames(features_and_rates.cc))])),
#                                                 y = 'true_optimal_const', training_frame = training_set.h2o)

# Evaluate performance of the regression fit on the test set
true_rate <- test_set$true_rate
# true_optimal_const <- test_set$true_optimal_const
rate_regression.predictions <- as.vector(predict(h2o.rate_regression_fit, test_set.h2o)$pred)
# rate_regression.predictions <- as.vector(h2o.predict(h2o.rate_regression_fit, test_set.h2o))
# constant_regression.predictions <- as.vector(h2o.predict(h2o.constant_regression_fit, test_set.h2o))
rate_regression.MSE <- mean((rate_regression.predictions - true_rate)^2)
# constant_regression.MSE <- mean((constant_regression.predictions - true_optimal_const)^2)

cat("MSE for rate = ", rate_regression.MSE, "\n")
# cat("MSE for constant = ", constant_regression.MSE, "\n")

# Train neural network on full data
cat('Now training the neural net on the full dataset\n')
cat("About to run the rate regression on full dataset\n")
if(useEnsemble){
  h2o.rate_regression_fit.full_data <- h2o.ensemble(x = setdiff(colnames(features_and_rates.cc),
                                                                c('dataset_id', 'id', 'true_rate',
                                                                  colnames(features_and_rates.cc)[grep(pattern = "loss|true|is_best", colnames(features_and_rates.cc))])),
                                                    y = 'true_rate', training_frame = full_dataset.h2o,
                                                    family = "gaussian", 
                                                    learner = learner, 
                                                    metalearner = metalearner,
                                                    cvControl = list(V = 5))
}else{
  h2o.rate_regression_fit.full_data <- h2o.deeplearning(x = setdiff(colnames(features_and_rates.cc),
                                                                    c('dataset_id', 'id', 'true_rate',
                                                                      colnames(features_and_rates.cc)[grep(pattern = "loss|true|true|is_best", colnames(features_and_rates.cc))])),
                                                        y = 'true_rate', training_frame = full_dataset.h2o)
}

# cat("About to run the constant regression on full dataset\n")
# h2o.constant_regression_fit.full_data <- h2o.deeplearning(x = setdiff(colnames(features_and_rates.cc),
#                                                                       c('dataset_id', 'id', 'true_rate',
#                                                                         colnames(features_and_rates.cc)[grep(pattern = "loss|true|true|is_best", colnames(features_and_rates.cc))])),
#                                                           y = 'true_optimal_const', training_frame = full_dataset.h2o)
if(useEnsemble){
  rate_model.path <- h2o.save_ensemble(h2o.rate_regression_fit.full_data, path = './rate_regression.ensemble_fit', force = TRUE, export_levelone = FALSE)
}else{
  rate_model.path <- h2o.saveModel(h2o.rate_regression_fit.full_data, path = './rate_regression.deeplearning_fit', force = TRUE)
}
# constant_model.path <- h2o.saveModel(h2o.constant_regression_fit.full_data, path = './constant_regression.deeplearning_fit', force = TRUE)