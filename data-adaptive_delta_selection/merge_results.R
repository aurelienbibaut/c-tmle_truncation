results_files <- list.files(pattern = "^broken_lines.results")
broken_lines.results <- vector()
for(results_file in results_files){
  current_results <- read.csv(results_file)
  if(!is.null(nrow(broken_lines.results)))
    current_results[, "dataset_id"] <- current_results[, "dataset_id"] + max(broken_lines.results[, "dataset_id"])
  broken_lines.results <- rbind(broken_lines.results, current_results)
}
write.csv(broken_lines.results, "merged_broken_lines.results.csv")
