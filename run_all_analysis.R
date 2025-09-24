#!/usr/bin/env Rscript
# Main analysis runner for Predictive Shrinkage project
cat("=== Predictive Shrinkage Analysis ===\n")

# Function to run ridge analysis
run_ridge_analysis <- function() {
  cat("\n" , rep("=", 50), "\n")
  cat("RUNNING RIDGE ANALYSIS\n")
  cat(rep("=", 50), "\n")
  
  # Set parameters for ridge
  assign("ALPHA", 0, envir = .GlobalEnv)
  
  # STEP 1: Run ridge simulations
  cat("Step 1: Running ridge simulations to generate data...\n")
  
  cat("- Joint cross-validation simulations...\n")
  source('scripts/run_simulations.R')
  
  cat("- Two-step procedure simulations...\n") 
  source('scripts/run_simulations_two_step.R')
  
  cat("- Copas procedure simulations...\n")
  source('scripts/run_simulations_copas.R')
  
  # Brief pause to ensure file system writes are complete
  Sys.sleep(1)
  
  # STEP 2: Run ridge analysis
  cat("Step 2: Running ridge analysis scripts...\n")
  cat("- Analyzing ridge regression results...\n")
  source('scripts/explore_ridge_regression_results.R')
  
  cat("RIDGE analysis complete!\n")
}

# Function to run lasso analysis (both comparison with sparsity and single sparsity)
run_lasso_analysis <- function() {
  cat("\n" , rep("=", 50), "\n")
  cat("RUNNING LASSO ANALYSIS (SPARSE & NON-SPARSE)\n")
  cat(rep("=", 50), "\n")
  
  # Set alpha for lasso
  assign("ALPHA", 1, envir = .GlobalEnv)
  
  # PHASE 1: Single sparsity value lasso
  cat("\n--- PHASE 1: Single Sparsity Value Lasso ---\n")
  cat("Using single sparsity value:", LASSO_SPARSITY_VALUE_SINGLE, "\n")
  
  # Temporarily set to single sparsity value
  original_sparsity <- LASSO_SPARSITY_VALUES
  LASSO_SPARSITY_VALUES <<- c(LASSO_SPARSITY_VALUE_SINGLE)
  
  cat("- Joint cross-validation simulations (non-sparse)...\n")
  source('scripts/run_simulations.R')
  
  cat("- Two-step procedure simulations (non-sparse)...\n") 
  source('scripts/run_simulations_two_step.R')
  
  cat("- Copas procedure simulations (non-sparse)...\n")
  source('scripts/run_simulations_copas.R')
  
  # PHASE 2: Multiple sparsity values lasso
  cat("\n--- PHASE 2: Multiple Sparsity Values Lasso ---\n")
  cat("Using multiple sparsity values:", paste(original_sparsity, collapse=", "), "\n")
  
  # Restore multiple sparsity values for sparse analysis
  LASSO_SPARSITY_VALUES <<- original_sparsity
  
  cat("- Joint cross-validation simulations (sparse)...\n")
  source('scripts/run_simulations.R')
  
  
  cat("- Two-step procedure simulations (sparse)...\n") 
  source('scripts/run_simulations_two_step.R')
  
  cat("- Copas procedure simulations (sparse)...\n")
  source('scripts/run_simulations_copas.R')
  
  # Brief pause to ensure file system writes are complete
  Sys.sleep(1)
  
  # STEP 2: Run lasso analysis scripts
  cat("\nStep 2: Running lasso analysis scripts...\n")
  
  cat("- Analyzing lasso regression results (sigma analysis)...\n")
  source('scripts/explore_lasso_regression_results.R')
  
  cat("- Analyzing lasso sparsity results...\n")
  source('scripts/explore_lasso_regression_results_vary_sparsity.R')
  
  cat("LASSO analysis complete!\n")
}

# Main execution
cat("Loading configurations...\n")
source('config.R')


cat("Loading functions...\n")
# Load all function libraries
source('R/utils.R')
source('R/single_response.R') 
source('R/multi_response.R')
source('R/multi_response_sparse.R')

cat("Setting up parameters...\n")
# Run setup scripts to generate fixed parameters
if (file.exists('scripts/setup/get_lasso_fixed_params.R')) {
  source('scripts/setup/get_lasso_fixed_params.R')
} else {
  cat("Warning: get_lasso_fixed_params.R not found\n")
}

if (file.exists('scripts/setup/get_ridge_fixed_params.R')) {
  source('scripts/setup/get_ridge_fixed_params.R')
} else {
  cat("Warning: get_ridge_fixed_params.R not found\n")
}

# Display configuration summary
cat("\nConfiguration Summary:\n")
cat("- Test mode:", TEST_MODE, "\n")
cat("- Ridge sigma values:", paste(RIDGE_SIGMA_VALUES, collapse=", "), "\n")
cat("- Lasso sigma values:", paste(LASSO_SIGMA_VALUES, collapse=", "), "\n")
cat("- Lasso single sparsity value:", LASSO_SPARSITY_VALUE_SINGLE, "\n")
cat("- Lasso multiple sparsity values:", paste(LASSO_SPARSITY_VALUES, collapse=", "), "\n")

# Run Ridge Analysis
run_ridge_analysis()

# Run Lasso Analysis (both sparse and non-sparse)
run_lasso_analysis()

cat("\nRunning additional simulations...\n")
# Multi-response simulations
if (file.exists('scripts/run_multi_response_sims(curds).R')) {
  cat("Running multi-response simulations...\n")
  source('scripts/run_multi_response_sims(curds).R')
} else {
  cat("Multi-response simulations script not found, skipping...\n")
}

if (file.exists('scripts/run_multi_response_sims(sparse_curds).R')) {
  cat("Running sparse multi-response simulations...\n")
  source('scripts/run_multi_response_sims(sparse_curds).R')
} else {
  cat("Sparse multi-response simulations script not found, skipping...\n")
}

cat("\n", rep("=", 60), "\n")
cat("ALL ANALYSIS COMPLETE!\n")
cat("Results generated:\n")
cat("- Ridge results: results/simulation_results/*ridge*\n")
cat("- Lasso non-sparse: results/simulation_results/*lasso* (main directories)\n") 
cat("- Lasso sparse: results/simulation_results/*/sparsity/\n")
cat("- Figures: results/figures/\n")
cat("- Multi-response: results/simulation_results/multi_response/\n")
cat(rep("=", 60), "\n")