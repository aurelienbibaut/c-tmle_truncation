library(h2o)
features_and_rates <- read.csv('merged_rate_inference.features.results.csv')[, -1]

features_and_rates <- cbind(features_and_rates, 
                            true_rate = 1 / (2 * (features_and_rates$true_gamma.1 + 1 - features_and_rates$true_beta)))

features_and_rates.cc <- features_and_rates[!is.na(features_and_rates$true_rate), ]

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

h2o.slope_regression_fit <- h2o.deeplearning(x = setdiff(colnames(features_and_rates.cc),
                                                         c('dataset_id', 'id', 'true_rate',
                                                           colnames(features_and_rates.cc)[grep(pattern = "loss|beta|true_gamma|is_best", colnames(features_and_rates.cc))])),
                                             y = 'true_rate', training_frame = training_set.h2o)

# Evaluate performance of the regression fit on the test set
true_rate <- test_set$true_rate
regression_predictions <- as.vector(h2o.predict(h2o.slope_regression_fit, test_set.h2o))
slope_regression.MSE <- mean((regression_predictions - true_rate)^2)

cat("MSE = ", slope_regression.MSE, "\n")

# Checking generalization
# newdatapoint <- generate_data_and_gamma_broken_line('L0_exp', 2, 4, -1, 1, 2, 1e4, 1/4, T, F)
# newdatapoint <- cbind(dataset_id = 1e4, newdatapoint)
# formatted_newdatapoint.h2o <- as.h2o(preprocess_dataset(newdatapoint)$dataset)
# h2o.predict(h2o.slope_regression_fit, formatted_newdatapoint.h2o)
