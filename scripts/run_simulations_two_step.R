# this script runs the simulations 
# for doing ridge (LASSO) then choosing K 
# with a second CV step. 
source('R/utils.R')
source('R/single_response.R')
source('R/multi_response.R')           # for multi-response scripts
source('R/multi_response_sparse.R')    # for sparse scripts
library(glmnet)
library(rlist)
set.seed(5636535)

# Set parameters from the configuration
alpha <- ALPHA       # 1 for LASSO, 0 for ridge

# Determine which parameters to use based on method
if (alpha == 0) {
  # Ridge regression
  sigma_values <- RIDGE_SIGMA_VALUES
  sparsity_values <- c(0)  # Ridge always uses 0 sparsity
  base_output_dir <- RIDGE_TWO_STEP_OUTPUT_DIR
  method <- 'ridge'
} else {
  # Lasso regression
  sigma_values <- LASSO_SIGMA_VALUES
  sparsity_values <- LASSO_SPARSITY_VALUES  # Uses modified values
  method <- 'lasso'
  
  # Determine output directory based on number of sparsity values
  if (length(sparsity_values) >= 2) {
    # Multiple sparsity values - use sparsity-specific directory
    base_output_dir <- LASSO_TWO_STEP_SPARSITY_OUTPUT_DIR
    cat("Using two-step sparsity-specific output directory for multiple sparsity values\n")
  } else {
    # Single sparsity value - use main two-step directory
    base_output_dir <- LASSO_TWO_STEP_OUTPUT_DIR
    cat("Using main two-step output directory for single sparsity value\n")
  }
}

# args <- commandArgs(TRUE)
# sigma <- as.double(args[1]) # std of error terms; we will vary this
# s <- as.double(args[2]) # sparsity
# alpha <- as.double(args[3]) # 0 if ridge, 1 if lasso
# out_dir <- as.character(args[4])
# 
# print(sigma)
# print(s)
# print(alpha)
# print(out_dir)

######################
# fixed parameters: 
######################
# Save the current value before loading
temp_sparsity <- LASSO_SPARSITY_VALUES

if(alpha == 0){
  load(paste0(FIXED_PARAMS_DIR, 'ridge_fixed_params.RData'))
}else if(alpha == 1){
  load(paste0(FIXED_PARAMS_DIR, 'lasso_fixed_params.RData'))
  # Restore the modified value after loading
  LASSO_SPARSITY_VALUES <<- temp_sparsity
}

# Use configuration values
n_trials <- N_TRIALS

# Loop through each sigma value
for (sigma_val in sigma_values) {
  sigma <- sigma_val
  
  # Loop through each sparsity value
  for (s_val in sparsity_values) {
    s <- s_val
    
    ######################
    # run simulations: 
    ######################
    
    # two rows, n_trials columns
    # first row stores the k's, second row stores the lambdas
    two_step_results_mat <- matrix(0, ncol = n_trials, nrow = 2)
    two_step_slope_cs <- rep(0, n_trials)
    two_step_slope_vs <- rep(0, n_trials)
    two_step_mspe_cs <- rep(0, n_trials)
    two_step_mspe_vs <- rep(0, n_trials)
    
    for(i in 1:n_trials){
      if(i %% 10 == 0){
        cat('trial', i, 'for sigma =', sigma, ', sparsity =', s, '\n')
      }
      data <- draw_cs_vs_sample(x_cs, beta_full, s, sigma, n_vs)  # draw data
      
      # K after CV - use configuration value for CV folds
      K <- get_two_step_cv_k(data$x_cs, data$y_cs, alpha = alpha, nfolds = CV_FOLDS)
      cv_fit <- cv.glmnet(data$x_cs, data$y_cs, alpha = alpha, nfolds = CV_FOLDS) 
      two_step_results_mat[1, i] <- K
      two_step_results_mat[2, i] <- cv_fit$lambda.min
      # get prediction
      two_step_cs_pred <- predict(cv_fit, newx = data$x_cs, 
                                  s = cv_fit$lambda.min) * K 
      two_step_vs_pred <- predict(cv_fit, newx = data$x_vs, 
                                  s = cv_fit$lambda.min) * K 
      
      two_step_slope_cs[i] <- drop(get_slope(two_step_cs_pred, data$y_cs))
      two_step_slope_vs[i] <- drop(get_slope(two_step_vs_pred, data$y_vs))
      two_step_mspe_cs[i] <- drop(get_mspe(two_step_cs_pred, data$y_cs))
      two_step_mspe_vs[i] <- drop(get_mspe(two_step_vs_pred, data$y_vs))
      
    }
    
    results_list <- 
      list(two_step_slope_cs = two_step_slope_cs,
           two_step_slope_vs = two_step_slope_vs,
           two_step_mspe_cs = two_step_mspe_cs,
           two_step_mspe_vs = two_step_mspe_vs, 
           two_step_results_mat = two_step_results_mat,
           sigma = sigma, 
           alpha = alpha, 
           s = s)
    
    # Create output directory if it doesn't exist
    if (!dir.exists(base_output_dir)) {
      dir.create(base_output_dir, recursive = TRUE)
      cat("Created output directory:", base_output_dir, "\n")
    }
    
    outfile <- paste0(base_output_dir,
                      method, '_two_step_sim_results_', 
                      'sigma', sigma, '_sparsity', s, '.rds')
    print(paste0('Done with sigma = ', sigma, ', sparsity = ', s, '. Saving results to ', outfile))
    list.save(results_list, outfile)
  }
}