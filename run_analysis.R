#!/usr/bin/env Rscript
# Main analysis runner for Predictive Shrinkage project

cat("=== Predictive Shrinkage Analysis ===\n")
cat("Loading functions...\n")

source('config.R')

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
source('scripts/run_simulations.R')

cat("Running multi-response simulations...\n") 
source('scripts/run_multi_response_sims.R')

cat("Running sparse multi-response simulations...\n")
source('scripts/run_sparse_sims.R')

cat("Running fMRI analysis...\n")
source('scripts/run_fmri_analysis.R')

cat("Analysis complete! Check results/ directory for outputs.\n")