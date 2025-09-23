test_function_correctness <- function() {
  # Test slope calculation
  set.seed(123)
  x <- c(1, 2, 3, 4, 5)
  y <- 2 * x + rnorm(5, 0, 0.1)  # y ≈ 2x
  slope <- get_slope(x, y)
  stopifnot(abs(slope - 2) < 0.5)  # Should be close to 2
  
  # Test MSPE calculation
  predicted <- c(1, 2, 3, 4, 5)
  observed <- c(1.1, 1.9, 3.2, 3.8, 5.1)
  mspe <- get_mspe(predicted, observed)
  expected_mspe <- mean((predicted - observed)^2)
  stopifnot(abs(mspe - expected_mspe) < 1e-10)
  
  # Test sparsity logic - use proper matrix dimensions
  n_obs <- 10
  n_dim <- 5
  beta_test <- c(1, 2, 3, 4, 5)  # Length matches dimensions
  x_test <- matrix(rnorm(n_obs * n_dim), nrow = n_obs, ncol = n_dim)  # 10x5 matrix
  
  data_sparse <- draw_cs_vs_sample(x_test, beta_test, s=2, sigma=1, n_vs=10)
  expected_beta <- c(0, 0, 3, 4, 5)  # First 2 should be zero
  stopifnot(all(data_sparse$true_beta == expected_beta))
  
  cat("✓ Function correctness tests passed\n")
}