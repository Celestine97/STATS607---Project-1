#!/usr/bin/env Rscript
# Main analysis runner for Predictive Shrinkage project

cat("=== Predictive Shrinkage Analysis ===\n")

cat("Loading configurations")
source('config.R')

cat("Loading functions...\n")
# Load all function libraries
source('R/utils.R')
source('R/single_response.R') 
source('R/multi_response.R')
source('R/multi_response_sparse.R')

cat("Setting up parameters...\n")
# Run setup scripts to generate fixed parameters
source('scripts/setup/get_lasso_fixed_params.R')
source('scripts/setup/get_ridge_fixed_params.R')  # if it exists

cat("Running single response simulations...\n")
# Run a subset of simulations (not full cluster scale)
cat("joint cross-validation")
source('scripts/run_simulations.R')
cat("two step")
source('scripts/run_simulations_two_step.R')
cat("copas")
source('scripts/run_simulations_copas.R')

# cat("Running multi-response simulations...\n") 
# source('scripts/run_multi_response_sims(curds).R')

# cat("Running sparse multi-response simulations...\n")
# source('scripts/run_multi_response_sims(sparse_curds).R')

cat("Running analysis")
source('scripts/explore_ridge_regression_results.R')
# source('scripts/explore_lasso_regression_results.R')
# source('scripts/explore_lasso_regression_results_vary_sparsity.R')

cat("Analysis complete! Check results/ directory for outputs.\n")