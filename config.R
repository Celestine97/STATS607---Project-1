# config.R - Project configuration
TEST_MODE <- TRUE  # Set to FALSE for full analysis

if(TEST_MODE) {
  # Test parameters - runs in minutes
  N_TRIALS <- 5
  K_RANGE_STEP <- 0.2        # Fewer K values to test
  CV_FOLDS <- 3              # Fewer folds
  SIGMA_VALUES <- c(5,10)       # Fewer noise level
  SNR_VALUES <- c(1.0,10.0)       # Test 2 SNR level
  cat("=== TEST MODE: Reduced parameters for quick execution ===\n")
} else {
  # Full analysis parameters - runs in hours
  N_TRIALS <- 500
  K_RANGE_STEP <- 0.05
  CV_FOLDS <- 10
  SIGMA_VALUES <- c(5, 10, 15, 20, 25)
  SNR_VALUES <- c(1.0, 10.0, 100.0)
  cat("=== FULL MODE: Complete analysis ===\n")
}

# Derived parameters
K_RANGE <- seq(0.5, 1.5, by = K_RANGE_STEP)

# 
ALPHA <- 1  # 1 for LASSO, 0 for ridge
SPARSITY <- c(0)

SAVE_FIG <- TRUE