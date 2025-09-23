# config.R - Project configuration
TEST_MODE <- TRUE  # Set to FALSE for full analysis

RIDGE_METHOD <- 0      # Ridge regression
LASSO_METHOD <- 1      # Lasso regression

if(TEST_MODE) {
  # Test parameters - runs in minutes
  N_TRIALS <- 10
  K_RANGE_STEP <- 0.1       # Fewer K values to test
  CV_FOLDS <- 5              # Fewer folds
  
  # Ridge parameters
  RIDGE_SIGMA_VALUES <- c(5, 10)
  
  # Lasso parameters  
  LASSO_SIGMA_VALUES <- c(5, 10)
  LASSO_SPARSITY_VALUES <- c(900, 850)
  LASSO_SPARSITY_VALUE_SINGLE <- 950
  
  cat("=== TEST MODE: Reduced parameters for quick execution ===\n")
} else {
  # Full analysis parameters - runs in hours
  N_TRIALS <- 500
  K_RANGE_STEP <- 0.05
  CV_FOLDS <- 10
  
  # Ridge parameters (from shell scripts)
  RIDGE_SIGMA_VALUES <- c(5, 10, 20, 30, 40, 60, 70, 80, 90, 100)
  
  # Lasso parameters (from shell scripts)
  LASSO_SIGMA_VALUES <- c(5, 11, 12, 13, 14, 15)
  LASSO_SPARSITY_VALUES <- c(980, 950, 920, 890, 860, 830)
  LASSO_SPARSITY_VALUE_SINGLE <- 950
  
  cat("=== FULL MODE: Complete analysis ===\n")
}

# Output directories
RIDGE_OUTPUT_DIR <- 'results/simulation_results/joint_cv_ridge_results/'
RIDGE_TWO_STEP_OUTPUT_DIR <- 'results/simulation_results/two_step_ridge_results/'
RIDGE_COPAS_OUTPUT_DIR <- 'results/simulation_results/copas_ridge_results/'

LASSO_OUTPUT_DIR <- 'results/simulation_results/joint_cv_lasso_results/'
LASSO_TWO_STEP_OUTPUT_DIR <- 'results/simulation_results/two_step_lasso_results/'
LASSO_COPAS_OUTPUT_DIR <- 'results/simulation_results/copas_lasso_results/'

LASSO_CV_SPARSITY_OUTPUT_DIR <- 'results/simulation_results/joint_cv_lasso_results/sparsity/'
LASSO_TWO_STEP_SPARSITY_OUTPUT_DIR <- 'results/simulation_results/two_step_lasso_results/sparsity/'
LASSO_COPAS_SPARSITY_OUTPUT_DIR <- 'results/simulation_results/copas_lasso_results/sparsity/'

# Additional directories
FIXED_PARAMS_DIR <- 'results/simulation_results/fixed_params/'
FIGURES_DIR <- 'results/figures/'
RIDGE_FIGURES_DIR <- 'results/figures/ridge_figures/'
LASSO_FIGURES_DIR <- 'results/figures/lasso_figures/'

# for multi-response
MULTI_RESPONSE_DIR <- 'results/simulation_results/multi_response/'
MULTI_RESPONSE_FIGURE_DIR <- 'results/figures/multi_response/'

# Derived parameters
K_RANGE <- seq(0.5, 1.5, by = K_RANGE_STEP)

# Additional parameters
SAVE_FIG <- TRUE

# Function to get all parameter combinations
get_ridge_params <- function() {
  list(
    sigma_values = RIDGE_SIGMA_VALUES,
    sparsity = 0,
    method = RIDGE_METHOD,
    output_dirs = list(
      main = RIDGE_OUTPUT_DIR,
      two_step = RIDGE_TWO_STEP_OUTPUT_DIR,
      copas = RIDGE_COPAS_OUTPUT_DIR
    ),
    fixed_params_dir = FIXED_PARAMS_DIR,
    figures_dir = RIDGE_FIGURES_DIR
  )
}

get_lasso_params <- function() {
  # Use the current global value if it exists, otherwise use config default
  current_sparsity <- if(exists("LASSO_SPARSITY_VALUES") && 
                         !identical(LASSO_SPARSITY_VALUES, c(50, 100))) {
    LASSO_SPARSITY_VALUES
  } else {
    if(TEST_MODE) c(50, 100) else c(980, 950, 920, 890, 860, 830)
  }
  
  list(
    sigma_values = LASSO_SIGMA_VALUES,
    sparsity_values = current_sparsity,  # Use the current value
    method = LASSO_METHOD,
    output_dirs = list(
      main = LASSO_OUTPUT_DIR,
      two_step = LASSO_TWO_STEP_OUTPUT_DIR,
      copas = LASSO_COPAS_OUTPUT_DIR,
      cv_sparsity = LASSO_CV_SPARSITY_OUTPUT_DIR,
      two_step_sparsity = LASSO_TWO_STEP_SPARSITY_OUTPUT_DIR,
      copas_sparsity = LASSO_COPAS_SPARSITY_OUTPUT_DIR
    ),
    fixed_params_dir = FIXED_PARAMS_DIR,
    figures_dir = LASSO_FIGURES_DIR
  )
}

# Print configuration summary
cat("Configuration loaded:\n")
cat("- Test mode:", TEST_MODE, "\n")
cat("- N_TRIALS:", N_TRIALS, "\n")
cat("- Ridge sigma values:", paste(RIDGE_SIGMA_VALUES, collapse=", "), "\n")
cat("- Lasso sigma values:", paste(LASSO_SIGMA_VALUES, collapse=", "), "\n")
if(!TEST_MODE || length(LASSO_SPARSITY_VALUES) > 2) {
  cat("- Lasso sparsity values:", paste(LASSO_SPARSITY_VALUES, collapse=", "), "\n")
}