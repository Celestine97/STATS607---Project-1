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

args <- commandArgs(TRUE)
sigma <- as.double(args[1]) # std of error terms; we will vary this
s <- as.double(args[2]) # sparsity
alpha <- as.double(args[3]) # 0 if ridge, 1 if lasso
out_dir <- as.character(args[4])

print(sigma)
print(s)
print(alpha)
print(out_dir)

######################
# fixed parameters: 
######################
if(alpha == 0){
  load('../simulation_results/fixed_params/ridge_fixed_params.RData')
  method <- 'ridge'
}else if(alpha == 1){
  load('../simulation_results/fixed_params/lasso_fixed_params.RData')
  method <- 'lasso'
}else{
  stop('alpha should be equal to 1 or 0 (values in between 
       not implemented/tested yet)')
}

######################
# run simulations: 
######################
n_trials <- 100

# two rows, n_trials columns
# first row stores the k's, second row stores the lambdas
two_step_results_mat <- matrix(0, ncol = n_trials, nrow = 2)

two_step_slope_cs <- rep(0, n_trials)
two_step_slope_vs <- rep(0, n_trials)
two_step_mspe_cs <- rep(0, n_trials)
two_step_mspe_vs <- rep(0, n_trials)

for(i in 1:n_trials){
  if(i %% 10 == 0){
    cat('trial', i, '\n')
  }
  data <- draw_cs_vs_sample(x_cs, beta_full, s, sigma, n_vs)  # draw data
  
  # K after CV
  K <- get_two_step_cv_k(data$x_cs, data$y_cs, alpha = alpha, nfolds = 10)
  cv_fit <- cv.glmnet(data$x_cs, data$y_cs, alpha = alpha, nfolds = 10) 
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


outfile <- paste0(out_dir,
                  method, '_two_step_sim_results_', 
                  'sigma', sigma, '_sparsity', s, '.rds')
print(paste0('Done. Saving results to ', outfile))
list.save(results_list, outfile)
