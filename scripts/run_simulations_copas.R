# this script runs the simulations 
# for doing the copas procedure to choose K 
source('R/utils.R')
source('R/single_response.R')
source('R/multi_response.R')           # for multi-response scripts
source('R/multi_response_sparse.R')    # for sparse scripts
library(glmnet)
library(rlist)
set.seed(38989201)

# Set parameters from the configuration
alpha <- ALPHA       # 1 for LASSO, 0 for ridge

# Determine which parameters to use based on method
if (alpha == 0) {
  # Ridge regression
  sigma_values <- RIDGE_SIGMA_VALUES
  sparsity_values <- c(0)  # Ridge always uses 0 sparsity
  base_output_dir <- RIDGE_COPAS_OUTPUT_DIR
  method <- 'ridge'
} else {
  # Lasso regression
  sigma_values <- LASSO_SIGMA_VALUES
  sparsity_values <- LASSO_SPARSITY_VALUES  # This will use the modified value
  method <- 'lasso'
  
  # Determine output directory based on number of sparsity values
  if (length(sparsity_values) >= 2) {
    # Multiple sparsity values - use sparsity-specific directory
    base_output_dir <- LASSO_COPAS_SPARSITY_OUTPUT_DIR
    cat("Using Copas sparsity-specific output directory for multiple sparsity values\n")
  } else {
    # Single sparsity value - use main Copas directory
    base_output_dir <- LASSO_COPAS_OUTPUT_DIR
    cat("Using main Copas output directory for single sparsity value\n")
  }
}

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
    
    copas_k <- rep(0, n_trials)
    copas_slope_cs <- rep(0, n_trials)
    copas_slope_vs <- rep(0, n_trials)
    copas_mspe_cs <- rep(0, n_trials)
    copas_mspe_vs <- rep(0, n_trials)
    
    for(i in 1:n_trials){
      if(i %% 10 == 0){
        cat('trial', i, 'for sigma =', sigma, ', sparsity =', s, '\n')
      }
      # draw data
      data <- draw_cs_vs_sample(x_cs, beta_full, s, sigma, n_vs)  # draw data
      
      if(alpha == 0){
        # ridge case
        stopifnot(dim(data$x_cs)[1] > dim(data$x_cs)[2])
        x <- data$x_cs
        x_new <- data$x_vs
      }else if (alpha == 1){
        # lasso case with improved feature selection
        attempts <- 0
        nonzero_indx <- c()
        
        while(length(nonzero_indx) < 2 && attempts < 10){
          attempts <- attempts + 1
          if(attempts > 1) {
            # Redraw data if previous attempt failed
            data <- draw_cs_vs_sample(x_cs, beta_full, s, sigma, n_vs) 
          }
          
          lasso_fit <- cv.glmnet(data$x_cs, data$y_cs)
          nonzero_coefs <- coef(lasso_fit, s = lasso_fit$lambda.min)
          nonzero_indx <- which(nonzero_coefs[-1] != 0)  # Remove intercept, check for non-zero
          
          if(length(nonzero_indx) < 2 && attempts <= 3) {
            cat("Attempt", attempts, ": Only", length(nonzero_indx), "features selected, retrying...\n")
          }
        }
        
        # If still insufficient features, use fallback method
        if(length(nonzero_indx) < 2) {
          cat("Using fallback: selecting top features by coefficient magnitude\n")
          all_coefs <- abs(coef(lasso_fit, s = lasso_fit$lambda.min)[-1])  # Remove intercept
          nonzero_indx <- order(all_coefs, decreasing = TRUE)[1:max(2, min(10, length(all_coefs)))]
        }
        
        # Ensure we have valid indices
        nonzero_indx <- nonzero_indx[nonzero_indx > 0 & nonzero_indx <= ncol(data$x_cs)]
        
        x <- data$x_cs[, nonzero_indx, drop = FALSE]
        x_new <- data$x_vs[, nonzero_indx, drop = FALSE]
        
        # Final validation
        stopifnot(ncol(x) >= 2)
        stopifnot(ncol(x_new) >= 2)
        stopifnot(nrow(x) > ncol(x))
      }
      
      # run OLS 
      ols_fit <- glmnet(x, data$y_cs, lambda = 0)
      # get copas K
      sigma2_hat <- compute_sigma2_hat(data$y_cs, x, ols_fit$beta)
      K <- get_copas_k(x, ols_fit$beta, sigma2_hat)
      
      # get prediction
      copas_cs_pred <- predict(ols_fit, newx = x) * K 
      copas_vs_pred <- predict(ols_fit, newx = x_new) * K 
      
      copas_slope_cs[i] <- drop(get_slope(copas_cs_pred, data$y_cs))
      copas_slope_vs[i] <- drop(get_slope(copas_vs_pred, data$y_vs))
      copas_mspe_cs[i] <- drop(get_mspe(copas_cs_pred, data$y_cs))
      copas_mspe_vs[i] <- drop(get_mspe(copas_vs_pred, data$y_vs))
      
      copas_k[i] <- K
    }
    
    results_list <- 
      list(copas_slope_cs = copas_slope_cs,
           copas_slope_vs = copas_slope_vs,
           copas_mspe_cs = copas_mspe_cs,
           copas_mspe_vs = copas_mspe_vs, 
           copas_k = copas_k,
           sigma = sigma, 
           alpha = alpha, 
           s = s)
    
    # Create output directory if it doesn't exist
    if (!dir.exists(base_output_dir)) {
      dir.create(base_output_dir, recursive = TRUE)
      cat("Created output directory:", base_output_dir, "\n")
    }
    
    outfile <- paste0(base_output_dir,
                      method, '_copas_sim_results_', 
                      'sigma', sigma, '_sparsity', s, '.rds')
    print(paste0('Done with sigma = ', sigma, ', sparsity = ', s, '. Saving results to ', outfile))
    list.save(results_list, outfile)
  }
}