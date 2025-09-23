#!/usr/bin/env Rscript
source('config.R')
source('R/utils.R')
source('R/single_response.R')

# Load test functions
source('tests/test_functions.R')
source('tests/test_pipeline.R')
source('tests/test_statistical_validation.R')

cat("=== Running Predictive Shrinkage Tests ===\n")

tryCatch({
  test_function_correctness()
  test_pipeline_integrity()
  test_statistical_validation()
  
  cat("All tests passed!\n")
}, error = function(e) {
  cat("Test failed:", e$message, "\n")
  quit(status = 1)
})