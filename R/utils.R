require(MASS)
require(caret)

# Data generating functions and evaluation functions

###########################
# functions to draw data 
###########################
draw_linear_responses <- function(x, beta, sigma){
  # draws y ~ Normal(x\beta, sigma^2)
  # given covariate matrix x, regression coefficients beta
  # and standard errors sigma, 
  
  n_obs <- dim(x)[1]
  return(rnorm(n_obs, mean = drop(x %*% beta), sd = sigma))
}

draw_cs_vs_sample <- function(x_cs, beta_full, s, sigma, n_vs){
  # returns the x and y's for both the construction and validation sample
  
  # x_cs are the construction sample covariates
  # s is the sparsity: we set s entries of beta_full to be zero
  # sigma is standard error of observations about regression line
  # n_vs is the number of samples in the validation set
  
  # quick check that dimensions work
  stopifnot(dim(x_cs)[2] == length(beta_full)) 
  
  n_cs <- dim(x_cs)[1] # number of samples in construction sample
  
  d <- length(beta_full) # dimension of regression coefficients
  stopifnot(s <= d) # check sparsity is less than dimension 
  
  # the true beta
  beta_sparse <- beta_full
  
  if(s > 0){
    beta_sparse[1:s] <- 0
  }else{
    beta_sparse <- beta_full
  }
  
  
  y_cs <- draw_linear_responses(x_cs, beta_sparse, sigma)
  
  V <- t(x_cs) %*% x_cs / n_cs # empirical correlation
  
  # draw validation sample x and y
  x_vs <- mvrnorm(n_cs, rep(0, d), V) 
  y_vs <- draw_linear_responses(x_vs, beta_sparse, sigma)
  return(list(x_vs = x_vs, 
              y_vs = y_vs, 
              x_cs = x_cs, 
              y_cs = y_cs, 
              true_beta = beta_sparse))
}


###############################
# Functions to evaluate predictive accuracy 
###############################
get_slope <- function(predicted_y, obs_y){
  # returns the slope of observed y on predicted y
  
  # observed y should be a vector
  # predicted y can be a matrix, with rows corresponding to 
  # number of observations
  # and columns corresponding to different predictions 
  # (predictions from different lambdas, for example)
  
  predicted_y_centered <- scale(predicted_y, center = TRUE, scale = FALSE)
  obs_y_centered <- scale(obs_y, center = TRUE, scale = FALSE)
  
  return(t(obs_y_centered) %*% predicted_y / 
           diag(t(predicted_y_centered) %*% predicted_y_centered))
  
}

get_mspe <- function(predicted_y, obs_y){
  if(is.vector(predicted_y) && is.vector(obs_y)) {
    return(mean((predicted_y - obs_y)^2))
  } else {
    return(colMeans((predicted_y - obs_y)^2))
  }
}