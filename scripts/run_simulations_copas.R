# this script runs the simulations 
# for doing the copas procedure to choose K 
source('R/utils.R')
source('R/single_response.R')
source('R/multi_response.R')           # for multi-response scripts
source('R/multi_response_sparse.R')    # for sparse scripts
library(glmnet)
library(rlist)
set.seed(38989201)

# Set default parameters for direct execution
alpha <- ALPHA       # 1 for LASSO, 0 for ridge

if (alpha == 0) {
  out_dir <- "results/simulation_results/copas_ridge_results/"
} else {
  out_dir <- "results/simulation_results/copas_lasso_results/"
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
if(alpha == 0){
  load('results/simulation_results/fixed_params/ridge_fixed_params.RData')
  method <- 'ridge'
}else if(alpha == 1){
  load('results/simulation_results/fixed_params/lasso_fixed_params.RData')
  method <- 'lasso'
}else{
  stop('alpha should be equal to 1 or 0 (values in between 
       not implemented/tested yet)')
}

# n_trials <- 40 
n_trials <- N_TRIALS

# Loop through each sigma value
for (sigma_val in SIGMA_VALUES) {
  sigma <- sigma_val
  
  # Loop through each sparsity value
  for (s_val in SPARSITY) {
    s <- s_val
    
    ######################
    # run simulations: 
    ######################
    
    # two rows, n_trials columns
    # first row stores the k's, second row stores the lambdas
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
        nonzero_indx <- 0
        while(sum(nonzero_indx) < 2){
          # lasso case
          # do feature selection
          data <- draw_cs_vs_sample(x_cs, beta_full, s, sigma, n_vs) 
          lasso_fit <- cv.glmnet(data$x_cs, data$y_cs)
          nonzero_indx <- which(coef(lasso_fit, s = lasso_fit$lambda.min) > 0)
          
          x <- data$x_cs[, nonzero_indx]
          x_new <- data$x_vs[, nonzero_indx]
          stopifnot(dim(x)[1] > dim(x)[2])
        }
        stopifnot(dim(x)[2] >= 2)
        stopifnot(dim(x_new)[2] >= 2)
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
    
    outfile <- paste0(out_dir,
                      method, '_copas_sim_results_', 
                      'sigma', sigma, '_sparsity', s, '.rds')
    print(paste0('Done with sigma = ', sigma, ', sparsity = ', s, '. Saving results to ', outfile))
    list.save(results_list, outfile)
  }
}