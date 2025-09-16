require(MASS)
require(caret)

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
  return(colMeans((predicted_y - obs_y)**2))
}

#################################
# functions to choose lambda
##################################
choose_lambda_cv <- function(x, y, lambdas, alpha = 1, 
                             nfolds = 10, plot = FALSE){
  # given x, y observations, returns the lambda chosen by cross validation
  # alpha 1 corresponds to lasso, alpha 0 is ridge
  
  cv_fit <- cv.glmnet(x, y, lambda = lambdas, alpha = alpha, nfolds = nfolds)
  
  if(plot){
    plot(cv_fit)
  }
  
  return(cv_fit$lambda.min)
}

choose_lambda_ic <- function(x, y, lambdas){
  # returns the optimal lambdas corresponding to minimizing aic and bic
  
  lasso_fit <- glmnet(x, y, lambda = lambdas)
  
  # predict new y and get RSS
  pred_y <- predict(lasso_fit, x = data$x_cs, y = data$y_cs, 
                    newx = data$x_cs, s = lambdas, exact = TRUE)
  rss <- get_mspe(pred_y, y)
  
  # get degrees of freedom 
  # TODO: there is a more correct version of df for LASSO ... 
  k <- lasso_fit$df
  
  # get aic and bic
  n <- length(y)
  aic <- 2 * k + n * log(rss)
  bic <- log(n) * k + n * log(rss)
  
  # get optimal lambda
  lambda_aic <- lambdas[which.min(aic)]
  lambda_bic <- lambdas[which.min(bic)]
  return(list(lambda_aic = lambda_aic, 
              lambda_bic = lambda_bic, 
              rss = rss, 
              k = k, 
              aic = aic, 
              bic = bic))
}

choose_joint_lambda_k <- function(x, y, lambdas, k_range, 
                                  alpha = 1, nfolds = 10){
  # given features x and responses y chooses lambda and K by cross-validation
  # lambdas and k_range specify the values of K and lambda for which we will 
  # test 
  
  n_lambda <- length(lambdas)
  n_k <- length(k_range)
  n_obs <- dim(x)[1]
  
  folds <- createFolds(seq(1, n_obs), k = nfolds)
  
  # matrix of MSPEs: rows range over values of K, columns range over lambda
  mspe_k_lambda <- matrix(0, nrow = n_k, ncol = n_lambda)
  
  for(k in 1:n_k){
    print(k)
    mspe_k <-  rep(0, n_lambda) # vector in which we store the result for each k
    for(f in 1:nfolds){
      # split by fold
      x_train <- x[-folds[[f]], ]
      y_train <- y[-folds[[f]]]
      
      x_test <- x[folds[[f]], ]
      y_test <- y[folds[[f]]]
      
      # train and test
      lasso_fit <- glmnet(x_train, y_train, alpha = alpha)
      y_pred <- predict(lasso_fit, newx = x_test, s = lambdas)
      
      # update mspe
      mspe_k <- mspe_k + get_mspe(k_range[k] * y_pred, y_test)
    }
    mspe_k_lambda[k, ] <- mspe_k / nfolds
  }
  
  # find optimal lambda and K
  indx <- which(mspe_k_lambda == min(mspe_k_lambda), arr.ind = TRUE)
  return(list(lambda = lambdas[indx[2]], 
              k = k_range[indx[1]]))
}

get_two_step_cv_k <- function(x, y, alpha, nfolds = 10){
  n_obs <- dim(x)[1]
  folds <- createFolds(seq(1, n_obs), k = nfolds)
  
  K <- rep(0, nfolds)
  for(f in 1:nfolds){
    # split by fold
    x_train <- x[-folds[[f]], ]
    y_train <- y[-folds[[f]]]
    
    x_test <- x[folds[[f]], ]
    y_test <- y[folds[[f]]]
    
    # run ridge (lasso) 
    cv_fit <- cv.glmnet(x_train, y_train, alpha = alpha, nfolds = nfolds)
    
    # get prediction on held out set
    y_pred <- predict(cv_fit, newx = x_test, s = cv_fit$lambda.min)
    
    # get slope
    K[f] <- drop(get_slope(y_pred, y_test))
  }
  
  return(mean(K))
  
}

compute_sigma2_hat <- function(y, x, beta_hat){
  # mean squared error of residuals to estimate variance
  p <- length(beta_hat)
  n <- dim(x)[1]
  df <- n - p
  return(sum((y - x %*% beta_hat)^2) / df)
}

get_copas_k <- function(x, beta_hat, sigma2_hat){
  # returns copas shrinkage factor K
  # given x, the estimated regression coefficients beta_hat
  # and the estimated variance sigma2_hat
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  v <- t(x) %*% x / n
  k <- p - 2
  f <- drop(t(beta_hat) %*% v %*% beta_hat) * n / (p * sigma2_hat)
  
  return(1 - k / (f * p))
}
  
