library(latex2exp)
library(tidyverse)
library(glmnet)
set.seed(54235245)

# load the data
load("fmri_study/data/fMRIdata.RData")

n_samples <- dim(fit_feat)[1]
n_voxels <- dim(resp_dat)[2]

# split into training and test set
propn_test <- 0.3
test_indx <- sample.int(n_samples, size = round(propn_test * n_samples), 
                        replace = FALSE)
X_train <- fit_feat[-test_indx, ]
y_train <- resp_dat[-test_indx, ]

X_test <- fit_feat[test_indx, ]
y_test <- resp_dat[test_indx, ]
n_test <- dim(X_test)[1]


# define the training method
lasso_fit <- function(covariates, response, lambda = NULL) {
  glmnet(
    covariates,
    response,
    family = "gaussian",
    alpha = 1,
    lambda = lambda
  )
}

# fit a lasso model on each voxel
fit_per_voxel <- lapply(seq(20), function(voxel) {
  lasso_fit(X_train, y_train[, voxel])
})

# given:  a matrix of covariates, a vector of responses, 
#         a function for training a prediction model, 
#         a measure of performance on held out data, and 
#         an integer number of folds K,
# return: the kpi measure calculated on held out data, 
#         calculated on held out data & averaged over K folds
cv <- function(covariates, responses, fit_func, kpi = cor, K = 5) {
  # total amount of data handed to CV
  n <- nrow(covariates)
  
  # split indices {1,...,n} evenly into K groups
  unif <- array(1 / K, K)
  partition <- sample(cut(seq(n),
                          n * cumsum(c(0, unif)),
                          labels = seq(length(unif))))
  
  # function to hold out one of K groups, 
  # fit a model, and evaluate predictions
  hold_out_k <- function(k) {
    # split data to hold out group k from training set
    train_x <- covariates[partition != k, ]
    train_y <- responses[partition != k]
    held_x  <- covariates[partition == k, ]
    held_y  <- responses[partition == k]
    
    # fit the model on the training data except group k
    fit_func(train_x, train_y) %>% 
      # get predictions for held out group k
      predict(newx = held_x) %>%
      # calculate performance measure on predictions
      kpi(held_y) %>%
      # convert to single number (from list-like object)
      as.numeric
  }
  
  # hold out each group and average evaluation of predictions
  rowMeans(sapply(seq(K), hold_out_k))
}

cv_lasso <- Reduce(rbind, lapply(seq(20), function(voxel) {
  lambda <- fit_per_voxel[[voxel]]$lambda
  data.frame(voxel, lambda,
             correlation = cv(X_train, y_train[, voxel],
                              partial(lasso_fit, lambda = lambda)))
}))

cv_lasso_best <- cv_lasso %>%
  filter((!is.na(lambda)) & (!is.na(correlation))) %>%
  group_by(voxel) %>%
  filter(abs(correlation) == max(abs(correlation))) %>%
  dplyr::select(voxel, lambda) %>%
  as.data.frame()
save(cv_lasso, file = 'fmri_study/outputs/cv_lasso.RData')
save(cv_lasso_best, file = 'fmri_study/outputs/cv_lasso_best.RData')

load("fmri_study/outputs/cv_lasso.RData")
load("fmri_study/outputs/cv_lasso_best.RData")

lasso_fit_final <- lapply(seq(20), function(voxel) {
  lambda <- cv_lasso_best[voxel, 'lambda']
  lasso_fit(X_train, y_train[, voxel], lambda = lambda)
})

A_lasso <- sapply(seq(20), function(voxel){as.vector(lasso_fit_final[[voxel]]$beta)}) %>% t

## Estimate the best linear predictor via cv
BLP_CV <- function(x, y, A_lasso, cv_lasso_best, nfolds = 5){
  N <- dim(y)[1] 
  q <- dim(y)[2]
  p <- dim(x)[2]
  fold <- sort(rep(1:nfolds, ceiling(N/nfolds)))[1:N]
  
  Bs <- lapply(1:nfolds, function(hold){
    train_x <- x[fold != hold, ]
    train_y <- y[fold != hold, ]
    valid_x <- x[fold == hold, ]
    valid_y <- y[fold == hold, ]
    ## fit LASSO on train, predict on validate
    valid_y_hat <- sapply(seq(20), function(voxel) {
      lambda <- cv_lasso_best[voxel, 'lambda']
      fit <- lasso_fit(X_train, y_train[, voxel], lambda = lambda)
      predict(fit, newx = valid_x) %>% as.vector 
    })
    
    ## fit OLS to fitted values, store weights
    sapply(1:q, function(k){
      as.vector(coef(lm(valid_y[, k] ~ valid_y_hat + 0)))
    }) %>% t
  })
  ## average coefficients across folds
  B <- Reduce('+', Bs) / nfolds
  B %*% A_lasso
}

A_blpcv_5 <- BLP_CV(X_train, y_train, A_lasso, cv_lasso_best, nfolds = 5)

save(A_lasso, file = 'fmri_study/outputs/A_lasso.RData')
save(A_blpcv_5, file = 'fmri_study/outputs/A_blpcv_5_folds.RData')


y_hat_lasso <- X_test %*% t(A_lasso)
y_hat_blpcv_5 <- X_test %*% t(A_blpcv_5)
y_hat_blpcv_5t <- X_test %*% t(A_blpcv_5t)
y_hat_blpcv_10 <- X_test %*% t(A_blpcv_10)
y_hat_blpcv_20 <- X_test %*% t(A_blpcv_20)






cor_test <- rbind(
  data.frame(
    Method = "LASSO",
    voxel = seq(20),
    correlation = sapply(1:20, function(v){
      cor(y_hat_lasso[, v], y_test[, v])
    })),
  data.frame(
    Method = "BLP-CV",
    voxel = seq(20),
    correlation = sapply(1:20, function(v){
      cor(y_hat_blpcv_5[, v], y_test[, v])
    }))
)

p3 <- cor_test %>%
  ggplot(aes(x = Method, y = correlation, fill = Method, color = Method)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  labs(x = "",
       y = "Correlation Between
Predicted/Actual Response",
       title = "Test Set Comparison, all 20 Voxels") + 
  theme(axis.text.x = element_blank(),
        text = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, size = 26),
        axis.text = element_text(size = 10)) +
  facet_wrap( ~ voxel)

p3

