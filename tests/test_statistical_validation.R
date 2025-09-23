test_statistical_validation <- function() {
  set.seed(42)
  
  cat("Testing statistical properties of shrinkage methods...\n")
  load(paste0(FIXED_PARAMS_DIR, 'lasso_fixed_params.RData'))
  
  n_sims <- 10
  
  # Test 1: Verify shrinkage actually shrinks coefficients
  cat("Testing coefficient shrinkage behavior...\n")
  shrinkage_test <- replicate(n_sims, {
    data <- draw_cs_vs_sample(x_cs[1:100, 1:200], beta_full[1:200], s=0, sigma=3, n_vs=100)
    
    # OLS coefficients
    ols_fit <- glmnet(data$x_cs, data$y_cs, lambda = 0)
    ols_coefs <- coef(ols_fit)[-1]  # Remove intercept
    
    # Ridge coefficients  
    ridge_fit <- cv.glmnet(data$x_cs, data$y_cs, alpha = 0)
    ridge_coefs <- coef(ridge_fit, s = ridge_fit$lambda.min)[-1]
    
    # Ridge coefficients should be smaller in magnitude
    mean(abs(ridge_coefs)) < mean(abs(ols_coefs))
  })
  
  shrinkage_rate <- mean(shrinkage_test, na.rm = TRUE)
  cat("Ridge shrinks coefficients in", round(shrinkage_rate * 100, 1), "% of cases\n")
  stopifnot(shrinkage_rate > 0.7)  # Should shrink most of the time
  
  # Test 2: Joint CV selects reasonable K values
  cat("Testing joint CV K selection...\n")
  k_test <- replicate(n_sims, {
    data <- draw_cs_vs_sample(x_cs[1:50, 1:100], beta_full[1:100], s=0, sigma=2, n_vs=50)
    
    tryCatch({
      cv_fit <- cv.glmnet(data$x_cs, data$y_cs, alpha = 1)
      joint_result <- choose_joint_lambda_k(data$x_cs, data$y_cs, 
                                            cv_fit$lambda, K_RANGE, 
                                            alpha = 1, nfolds = 3)
      
      # K should be positive and reasonable
      joint_result$k > 0 && joint_result$k < 2
    }, error = function(e) FALSE)
  })
  
  reasonable_k_rate <- mean(k_test, na.rm = TRUE)
  cat("Joint CV selects reasonable K in", round(reasonable_k_rate * 100, 1), "% of cases\n")
  stopifnot(reasonable_k_rate > 0.5)
  
  cat("✓ Statistical validation tests passed\n")
}