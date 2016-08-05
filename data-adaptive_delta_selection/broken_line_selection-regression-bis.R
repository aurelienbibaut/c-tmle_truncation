library(dummies)
library(h2o)
broken_lines.results <- read.csv('broken_lines.results.csv')

preprocess_dataset <- function(broken_lines.results){
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
  
  broken_lines.results.cc <- cbind(segment_id = kronecker(rep(1, nrow(broken_lines.results.cc) / 5), 1:5),
                                   broken_lines.results.cc)
  
  predictors <- c('relative_linearity',
                  'linearity',
                  'relative_segment_squared_length',
                  'nb_points',
                  'leftmost', 
                  'n', 'n_breakpoints',
                  'slopes.ordering')
  
  broken_lines.results.cc.colnames <- colnames(broken_lines.results.cc)
  broken_lines.results.cc <- as.data.frame(sapply(1:ncol(broken_lines.results.cc), 
                                                  function(j) as.numeric(broken_lines.results.cc[, j])))
  colnames(broken_lines.results.cc) <- broken_lines.results.cc.colnames
  
  
  wide_dataset <- reshape(data = broken_lines.results.cc[, c('dataset_id', 'segment_id', predictors, 'true_gamma', 'slope') ],
                          v.names = c(setdiff(predictors, 'n'), 'slope'), idvar = 'dataset_id', direction = 'wide',
                          timevar = 'segment_id')
  
  # Convert slopes orderings into factors
  for(j in 1:5){
    colname <- paste('slopes.ordering', j, sep = '.')
    wide_dataset[, colname] <- as.factor(wide_dataset[, colname])
  }
  
  # Add best segment feature
  wide_dataset <- cbind(wide_dataset, best_segment = NA)
  for(i in 1:nrow(wide_dataset)){
    current_dataset_id <- wide_dataset[i, 'dataset_id']
    wide_dataset[i, 'best_segment'] <- subset(broken_lines.results.cc, 
                                              dataset_id == current_dataset_id & is_best == 1, 
                                              select = segment_id)[[1]]
  }
  
  wide_dataset$best_segment <- as.factor(wide_dataset$best_segment)
  list(dataset = wide_dataset, complete_cases.dataset_ids = complete_cases.dataset_ids)
}

# Preprocess dataset
prepocessing_result <- preprocess_dataset(broken_lines.results)
wide_dataset <- prepocessing_result$dataset
complete_cases.dataset_ids <- prepocessing_result$complete_cases.dataset_ids

# Define test set and training set
training_set.dataset_ids <- sample(complete_cases.dataset_ids, floor(length(complete_cases.dataset_ids) * 4 / 5), F)
in_test_set.dataset_ids <- ! (complete_cases.dataset_ids %in% training_set.dataset_ids)
test_set.dataset_ids <- complete_cases.dataset_ids[in_test_set.dataset_ids]

training_set <- subset(wide_dataset, dataset_id %in% training_set.dataset_ids)
test_set <- subset(wide_dataset, dataset_id %in% test_set.dataset_ids)

localH2O <- h2o.init()
training_set.h2o <- as.h2o(training_set)
test_set.h2o <- as.h2o(test_set)
full_dataset.h2o <- as.h2o(wide_dataset)

h2o.segment_classification_fit <- h2o.deeplearning(x = setdiff(colnames(wide_dataset),
                                                               c('best_segment', 'dataset_id', 'slopes', 'true_gamma')), 
                                                   y = 'best_segment', training_frame = training_set.h2o)
h2o.slope_regression_fit <- h2o.deeplearning(x = setdiff(colnames(wide_dataset),
                                                         c('best_segment', 'dataset_id', 'true_gamma')), 
                                             y = 'true_gamma', training_frame = training_set.h2o)

# Evaluate performance of the classification fit on test set
classification_predictions <- as.numeric(as.vector(h2o.predict(h2o.segment_classification_fit, test_set.h2o)$predict))
slope.col_numbers <- which(colnames(test_set) %in% paste('slope', 1:5, sep = '.'))
predicted_slopes <- as.numeric(test_set[cbind(1:length(classification_predictions), slope.col_numbers[classification_predictions])])

# Classification error 
classification_matrix <- table(test_set$best,
                               classification_predictions,
                               dnn = c("actual", "predicted"))

classification_error_rate <- 1 - sum(diag(classification_matrix)) / sum(classification_matrix)

# MSE
true_gamma <- test_set$true_gamma
segment_classification.MSE <- mean((abs(predicted_slopes) - true_gamma)^2)

# Evaluate performance of the regression fit on the test set
regression_predictions <- as.vector(h2o.predict(h2o.slope_regression_fit, test_set.h2o))
slope_regression.MSE <- mean((regression_predictions - true_gamma)^2)

# Checking generalization
newdatapoint <- generate_data_and_gamma_broken_line('L0_exp', 2, 4, -1, 1, 2, 1e4, 1/4, T, F)
newdatapoint <- cbind(dataset_id = 1e4, newdatapoint)
formatted_newdatapoint.h2o <- as.h2o(preprocess_dataset(newdatapoint)$dataset)
h2o.predict(h2o.slope_regression_fit, formatted_newdatapoint.h2o)
