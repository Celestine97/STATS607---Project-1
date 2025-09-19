require(MASS)
source('R/utils.R')
source('R/single_response.R')
source('R/multi_response.R')           # for multi-response scripts
source('R/multi_response_sparse.R')    # for sparse scripts
set.seed(46645)

d <- 1000 # total number of dimensions
n_cs <- 200 # number of observations in construction sample
n_vs <- 1000 # number of observations in validation sample

# draw X for construction sample. This will be fixed for our analysis 
covar_x <- diag(d) * 1. # true covariance of our draws of X
mu_x <- rep(0, d) # mean of our draws of X

x_cs_ <- mvrnorm(n_cs, mu_x, covar_x) 
x_cs <- scale(x_cs_, scale = FALSE, center = TRUE) # center the columns of X

# regression coefficients
beta_full <- rnorm(d, 0, 1)

outfile <- '../simulation_results/fixed_params/lasso_fixed_params.RData'
print(paste0('saving to: ', outfile))
save.image(outfile)

