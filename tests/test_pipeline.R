test_pipeline_integrity <- function() {
  
  # Set test mode parameters
  old_test_mode <- TEST_MODE
  TEST_MODE <<- TRUE
  N_TRIALS <<- 2  # Minimal trials
  
  # Test configuration loading
  source('config.R')
  
  # Test required directories exist
  required_dirs <- c(RIDGE_OUTPUT_DIR, LASSO_OUTPUT_DIR, FIXED_PARAMS_DIR)
  for(dir in required_dirs) {
    if(!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    stopifnot(dir.exists(dir))
  }
  
  # Test ridge simulation runs and produces files
  ALPHA <<- 0
  source('scripts/run_simulations.R')
  ridge_files <- list.files(RIDGE_OUTPUT_DIR, pattern = "\\.rds$")
  stopifnot(length(ridge_files) > 0)
  
  # Test lasso simulation runs and produces files
  ALPHA <<- 1
  LASSO_SPARSITY_VALUES <<- c(950)  # Single value for test
  source('scripts/run_simulations.R')
  lasso_files <- list.files(LASSO_OUTPUT_DIR, pattern = "\\.rds$")
  stopifnot(length(lasso_files) > 0)
  
  # Test file content structure
  test_file <- paste0(LASSO_OUTPUT_DIR, lasso_files[1])
  results <- readRDS(test_file)
  required_fields <- c("ols_slope_vs", "lambda_only_slope_vs", "joint_slope_vs", 
                       "ols_mspe_vs", "lambda_only_mspe_vs", "joint_mspe_vs")
  stopifnot(all(required_fields %in% names(results)))
  stopifnot(length(results$ols_slope_vs) == N_TRIALS)
  
  # Restore original test mode
  TEST_MODE <<- old_test_mode
  
  cat("✓ Pipeline integrity tests passed\n")
}