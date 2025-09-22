# These are our original simulations 
# Run either lasso or ridge, get results 
# using OLS, traditional ridge (LASSO), 
# or joint cv-selected lambda and k

source('R/utils.R')
source('R/single_response.R')
source('R/multi_response.R')           # for multi-response scripts
source('R/multi_response_sparse.R')    # for sparse scripts
library(glmnet)
library(rlist)

set.seed(45345345)

# Set parameters from the configuration
alpha <- ALPHA       # 1 for LASSO, 0 for ridge

# Determine which parameters to use based on method
if (alpha == 0) {
  # Ridge regression
  sigma_values <- RIDGE_SIGMA_VALUES
  sparsity_values <- c(0)  # Ridge always uses 0 sparsity
  base_output_dir <- RIDGE_OUTPUT_DIR
  method <- 'ridge'
} else {
  # Lasso regression
  sigma_values <- LASSO_SIGMA_VALUES
  sparsity_values <- LASSO_SPARSITY_VALUES  # This will use the modified value
  method <- 'lasso'
  
  # Determine output directory based on number of sparsity values
  if (length(sparsity_values) >= 2) {
    # Multiple sparsity values - use sparsity-specific directory
    base_output_dir <- LASSO_CV_SPARSITY_OUTPUT_DIR
    cat("Using sparsity-specific output directory for multiple sparsity values\n")
  } else {
    # Single sparsity value - use main Copas directory
    base_output_dir <- LASSO_OUTPUT_DIR
    cat("Using main output directory for single sparsity value\n")
  }
}


# args <- commandArgs(TRUE)
# sigma <- as.double(args[1]) # std of error terms; we will vary this
# s <- as.double(args[2]) # sparsity
# alpha <- as.double(args[3]) # 0 if ridge, 1 if lasso
# out_dir <- as.character(args[4])

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
k_range <- K_RANGE

# Loop through each sigma value
for (sigma_val in sigma_values) {
  sigma <- sigma_val
  
  # Loop through each sparsity value
  for (s_val in sparsity_values) {
    s <- s_val
    
    ######################
    # run simulations: 
    # ######################
    
    # two rows, n_trials columns
    # first row stores the k's, second row stores the lambdas
    joint_results_mat <- matrix(0, ncol = n_trials, nrow = 2)
    rownames(joint_results_mat) <- c('k', 'lambda')
    cv_lambdas <- rep(0, n_trials)
    
    ols_slope_cs <- rep(0, n_trials)
    ols_slope_vs <- rep(0, n_trials)
    ols_mspe_cs <- rep(0, n_trials)
    ols_mspe_vs <- rep(0, n_trials)
    
    lambda_only_slope_cs <- rep(0, n_trials)
    lambda_only_slope_vs <- rep(0, n_trials)
    lambda_only_mspe_cs <- rep(0, n_trials)
    lambda_only_mspe_vs <- rep(0, n_trials)
    
    joint_slope_cs <- rep(0, n_trials)
    joint_slope_vs <- rep(0, n_trials)
    joint_mspe_cs <- rep(0, n_trials)
    joint_mspe_vs <- rep(0, n_trials)
    
    for(i in 1:n_trials){
      if(i %% 10 == 0){
        cat('trial', i, 'for sigma =', sigma, ', sparsity =', s, '\n')
      }
      data <- draw_cs_vs_sample(x_cs, beta_full, s, sigma, n_vs)  # draw data
      
      # if d < n compare with regression also
      if(dim(x_cs)[1] > dim(x_cs)[2]){
        ols_fit <- glmnet(data$x_cs, data$y_cs, lambda = 0)
        ols_cs_pred <- predict(ols_fit, newx = data$x_cs)
        ols_vs_pred <- predict(ols_fit, newx = data$x_vs)
        
        ols_slope_cs[i] <- drop(get_slope(ols_cs_pred, data$y_cs))
        ols_slope_vs[i] <- drop(get_slope(ols_vs_pred, data$y_vs))
        ols_mspe_cs[i] <- drop(get_mspe(ols_cs_pred, data$y_cs))
        ols_mspe_vs[i] <- drop(get_mspe(ols_vs_pred, data$y_vs))
      }
      
      # choose only lambda
      # Use configuration value for CV folds
      cv_fit <- cv.glmnet(data$x_cs, data$y_cs, alpha = alpha, nfolds = CV_FOLDS) 
      lambda_range <- cv_fit$lambda
      min_lambda <- cv_fit$lambda.min
      cv_lambdas[i] <- min_lambda 
      
      # get prediction
      lambda_only_cs_pred <- predict(cv_fit, newx = data$x_cs, s = min_lambda)
      lambda_only_vs_pred <- predict(cv_fit, newx = data$x_vs, s = min_lambda)
      
      # get results for fitting only lambda
      lambda_only_slope_cs[i] <- drop(get_slope(lambda_only_cs_pred, data$y_cs))
      lambda_only_slope_vs[i] <- drop(get_slope(lambda_only_vs_pred, data$y_vs))
      lambda_only_mspe_cs[i] <- drop(get_mspe(lambda_only_cs_pred, data$y_cs))
      lambda_only_mspe_vs[i] <- drop(get_mspe(lambda_only_vs_pred, data$y_vs))
      
      # choose joint lambda and k
      joint_results <- choose_joint_lambda_k(data$x_cs, data$y_cs, 
                                             lambda_range, k_range, 
                                             alpha = alpha, nfolds = CV_FOLDS)
      # if our default range for k isn't good enough ... 
      if(joint_results$k == max(k_range)){
        print('expanding k_range ... ')
        k_range_ <- seq(max(k_range), max(k_range) + 0.5, by = 0.05)
        joint_results <- 
          choose_joint_lambda_k(data$x_cs, data$y_cs, lambda_range, k_range_,
                                alpha = alpha, nfolds = CV_FOLDS)
      }else if(joint_results$k == min(k_range)){
        print('expanding k_range ... ')
        k_range_ <- seq(min(k_range) - 0.5, min(k_range), by = 0.05)
        joint_results <- choose_joint_lambda_k(data$x_cs, data$y_cs, 
                                               lambda_range, k_range_, 
                                               alpha = alpha, nfolds = CV_FOLDS)
      }
      
      joint_results_mat[1, i] <- joint_results$k
      joint_results_mat[2, i] <- joint_results$lambda
      
      # get prediction
      joint_cs_pred <- predict(cv_fit, newx = data$x_cs, 
                               s = joint_results$lambda) * joint_results$k 
      joint_vs_pred <- predict(cv_fit, newx = data$x_vs, 
                               s = joint_results$lambda) * joint_results$k 
      
      joint_slope_cs[i] <- drop(get_slope(joint_cs_pred, data$y_cs))
      joint_slope_vs[i] <- drop(get_slope(joint_vs_pred, data$y_vs))
      joint_mspe_cs[i] <- drop(get_mspe(joint_cs_pred, data$y_cs))
      joint_mspe_vs[i] <- drop(get_mspe(joint_vs_pred, data$y_vs))
      
    }
    
    results_list <- 
      list(ols_slope_cs = ols_slope_cs,
           ols_slope_vs = ols_slope_vs,
           ols_mspe_cs = ols_mspe_cs,
           ols_mspe_vs = ols_mspe_vs,
           lambda_only_slope_cs = lambda_only_slope_cs,
           lambda_only_slope_vs = lambda_only_slope_vs,
           lambda_only_mspe_cs = lambda_only_mspe_cs,
           lambda_only_mspe_vs = lambda_only_mspe_vs,
           joint_slope_cs = joint_slope_cs,
           joint_slope_vs = joint_slope_vs,
           joint_mspe_cs = joint_mspe_cs,
           joint_mspe_vs = joint_mspe_vs, 
           joint_results_mat = joint_results_mat, 
           cv_lambdas = cv_lambdas,
           sigma = sigma, 
           alpha = alpha,
           s = s)
    
    # Create output directory if it doesn't exist
    if (!dir.exists(base_output_dir)) {
      dir.create(base_output_dir, recursive = TRUE)
      cat("Created output directory:", base_output_dir, "\n")
    }
    
    outfile <- paste0(base_output_dir,
                      method, '_sim_results_', 
                      'sigma', sigma, '_sparsity', s, '.rds')
    print(paste0('Done with sigma = ', sigma, ', sparsity = ', s, '. Saving results to ', outfile))
    list.save(results_list, outfile)
  }
}
