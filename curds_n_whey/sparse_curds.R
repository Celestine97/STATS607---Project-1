library(tidyverse)
library(dplyr)
library(MASS)
library(Matrix)
library(expm)
library(quadprog)
library(glmnet)
library(corrplot)
library(PMA)


## generate multi-response sparse linear model (x, y)  
## following the Design section of Breiman-Friedman (1997)
generate_sparse_data <- function(p, q, N, SNR, rho, s, Sigma = 'eye'){
  # generate covariates x from 
  # AR(1) Gaussian graphical model
  r <- 2*runif(1) - 1
  V <- r**abs(outer(1:p, 1:p, "-"))
  x <- mvrnorm(N, (1:p)*0, V) 
  
  # generate true coefficients A
  K <- 10
  C <- mvrnorm(K, (1:q)*0, rho**abs(outer(1:q, 1:q, "-"))) 
  j <- sample(seq(s), K, replace = TRUE)
  l <- sample(seq(6), K, replace = TRUE)
  G <- pmax(Matrix(rep(l, s), nrow = s, byrow = TRUE) - abs(outer(1:s, j, "-")), 0)^2
  G <- G / Matrix(rep(colSums(G), s), nrow = s, byrow = TRUE)
  A <- G %*% C
  zeros <- Matrix(rep(0, (p-s)*q), ncol = q)
  A <- abs(rbind(A, zeros))
  # "all coefficient values were normalized by the same scale factor so that 
  # the average ('signal') variance for each response was equal to 1.0"
  signal <- diag(t(A) %*% V %*% A)
  A <- (A / sqrt(mean(signal))) %>% t 
  
  # generate responses y = Ax + eps
  sigma2 = 1/SNR
  if(Sigma == 'eye'){
    Sigma = sigma2 * diag(rep(1, q))
  } else { 
    Sigma = sigma2 * diag(seq(q)^2)
  }
  eps <- mvrnorm(N, (1:q)*0, Sigma) 
  y <- x %*% t(A) + eps
  list(x = x, A = A, V = V, y = y)
}

## Run LASSO separately on each problem
lasso_separate <- function(x, y, lambda = NULL){
  q <- dim(y)[2]
  out <- lapply(1:q, function(k){
    fit <- cv.glmnet(x, y[, k], intercept = FALSE, nfolds = 5, lambda = lambda)
    beta <- coef(fit, s = "lambda.min")
    list(beta = beta[2:nrow(beta)],
         lambdas = fit$lambda)}) 
  A_lasso <- sapply(1:q, function(k){out[[k]]$beta}) %>% t
  lambdas <- lapply(1:q, function(k){out[[k]]$lambdas})
  A_lasso
}

## Estimate union of support based on LASSO 
## then run Curds & Whey on support
support_cw <- function(x, y, A_lasso){
  # how many coefficients are non-zero at least k out of q
  votes <- sapply(0:q, function(k){
    length(which(colSums(A_lasso != 0) > k))
  })
  # do any of the votes satisfy the two constraints for C&W?
  valid_cutoffs <- (which((q+1 < votes) & (votes < n)) - 1)
  # if not, take all possible
  if(length(valid_cutoffs) == 0){
    support_est <- which(colSums(A_lasso != 0) > 0) 
    if(length(support_est) < q + 2){
      support_est <- c(support_est, sample((1:p)[-support_est], q+2-length(support_est)))
    } else{
      support_est <- sample(support_est, q+2)
    }
  } else {
    cutoff <- median(valid_cutoffs) # super arbitrary
    support_est <- which(colSums(A_lasso != 0) > cutoff) 
  }
  
  A_cw <- curds_n_whey_gcv(x[, support_est], y, A_lasso[, support_est])
  A_est <- matrix(1:(p*q), nrow=q)*0
  A_est[, support_est] <- A_cw
  A_est
}


## Shrink multi-response predictions according to  
## "generalized CV" - selecting diagonal matrix w/ CV
sparse_cw_gcv <- function(x, y, A_lasso, nfolds = 5){
  N <- dim(y)[1] 
  q <- dim(y)[2]
  p <- dim(x)[2]
  fold <- sort(rep(1:nfolds, ceiling(N/nfolds)))[1:N]
  
  y_hat_lasso <- x %*% t(A_lasso)
  
  T_cv <- lapply(1:nfolds, function(hold){
    train_x <- x[fold != hold, ]
    train_y <- y[fold != hold, ]
    train_y_s <- scale(train_y)
    sqrtm(solve(t(train_y_s) %*% train_y_s)) %*% CCA(train_x, train_y, K = q, penaltyz = 1, trace = FALSE)$v
  })
  
  T_inv_cv <- Map(solve, T_cv)
  R_cv <- sapply(1:N, function(n){
    T_cv[[fold[n]]] %*% y_hat_lasso[n, ]
  })
  
  dvec <- sapply(1:N, function(n){
    (y[n, ] %*% T_inv_cv[[fold[n]]]) * R_cv[, n]
  }) %>% rowSums
  
  Dmat <- Reduce('+', 
                 lapply(1:N, function(n){
                   (t(T_inv_cv[[fold[n]]]) %*% T_inv_cv[[fold[n]]]) * outer(R_cv[, n], R_cv[, n])
                 }))
  
  D_ <- diag(pmax(solve(Dmat) %*% dvec, 0) %>% as.vector())
  
  y_s <- scale(y)
  T_ <- sqrtm(solve(t(y_s) %*% y_s)) %*% CCA(x, y, K = q, penaltyz = 1, trace = FALSE)$v
  B <- solve(T_) %*% D_ %*% T_
  
  B %*% A_lasso
}

## Estimate the best linear predictor via cv
BLP_CV <- function(x, y, A_lasso, nfolds = 5){
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
    valid_y_hat <- valid_x %*% t(lasso_separate(train_x, train_y))
    
    ## fit OLS to fitted values, store weights
    sapply(1:q, function(k){
      as.vector(coef(lm(valid_y[, k] ~ valid_y_hat + 0)))
    }) %>% t
  })
  ## average coefficients across folds
  B <- Reduce('+', Bs) / nfolds
  B %*% A_lasso
}

## estimate autocorrelation
autocor <- function(x){
  p <- dim(x)[1]
  mean(sapply(1:(p-1), function(j){cor(x[, j], x[, j+1])}))}

## tridiagonal matrix - from stackoverflow.com/questions/28974507
tridiag <- function(upper, lower, main){
  out <- matrix(0,length(main),length(main))
  diag(out) <- main
  indx <- seq.int(length(upper))
  out[cbind(indx+1,indx)] <- lower
  out[cbind(indx,indx+1)] <- upper
  return(out)
}

## estimate inverse covariance of covariates
## with knowledge of AR(1) structure
precision_est <- function(x, y, A_lasso){
  p <- dim(x)[2]
  n <- dim(x)[1]
  r_est <- autocor(x)
  
  # closed form precision estimate
  off_diag <- -r_est / (1 - r_est^2)
  on_diag <- c(1+r_est^2, rep(1+(1-r_est)^2, p-2), 1+r_est^2)
  prec_xx <- tridiag(rep(off_diag, (p-1)), rep(off_diag, (p-1)), on_diag)
  
  x_s <- scale(x, scale = FALSE)
  y_s <- scale(y, scale = FALSE)
  cov_xy <- t(x_s) %*% y_s
  prec_yy <- solve(t(y_s) %*% y_s)
  Q <- prec_yy %*% t(cov_xy) %*% prec_xx %*% cov_xy
  r <- (sum(A_lasso != 0) *log(p) / q) / n
  B <- solve((1-r)*diag(rep(1, q)) + r*t(solve(Q)))
  B %*% A_lasso
}

## Generate data, estimate coefficients, and calculate
## metrics used by Breiman & Friedman (1997)
sim <- function(p, q, N, SNR, rho, s){
  data <- generate_sparse_data(p, q, N, SNR, rho, s) 
  A <- data$A
  V <- data$V
  x <- data$x
  y <- data$y
  
  A_lasso <- lasso_separate(x, y)
  A_cw <- sparse_cw_gcv(x, y, A_lasso)
  A_sup <- support_cw(x, y, A_lasso)
  A_prec <- precision_est(x, y, A_lasso)
  A_blpcv <- BLP_CV(x, y, A_lasso)
  
  MSE_lasso <- diag((A - A_lasso) %*% V %*% t(A - A_lasso))
  MSE_cw <- diag((A - A_cw) %*% V %*% t(A - A_cw))
  MSE_sup <- diag((A - A_sup) %*% V %*% t(A - A_sup))
  MSE_prec <- diag((A - A_prec) %*% V %*% t(A - A_prec))
  MSE_blpcv <- diag((A - A_blpcv) %*% V %*% t(A - A_blpcv))
  
  cw <- data.frame(MSE_relative = c(sum(MSE_cw) / sum(MSE_lasso)), 
                   MSE_individual = c(mean(MSE_cw / MSE_lasso)), 
                   MSE_worst = c(max(MSE_cw / MSE_lasso)), 
                   Method = c('C&W (Sparse CCA)'), 
                   p = p, q = q, N = N, SNR = SNR, rho = rho, s = s)
  
  sup <- data.frame(MSE_relative = c(sum(MSE_sup) / sum(MSE_lasso)), 
                   MSE_individual = c(mean(MSE_sup / MSE_lasso)), 
                   MSE_worst = c(max(MSE_sup / MSE_lasso)), 
                   Method = c('C&W (Subset Selection)'), 
                   p = p, q = q, N = N, SNR = SNR, rho = rho, s = s)
  
  prec <- data.frame(MSE_relative = c(sum(MSE_prec) / sum(MSE_lasso)), 
                    MSE_individual = c(mean(MSE_prec / MSE_lasso)), 
                    MSE_worst = c(max(MSE_prec / MSE_lasso)), 
                    Method = c('C&W (Precision Estimation)'), 
                    p = p, q = q, N = N, SNR = SNR, rho = rho, s = s)
  
  blpcv <- data.frame(MSE_relative = c(sum(MSE_blpcv) / sum(MSE_lasso)), 
                      MSE_individual = c(mean(MSE_blpcv / MSE_lasso)), 
                      MSE_worst = c(max(MSE_blpcv / MSE_lasso)), 
                      Method = c('BLP - CV'), 
                      p = p, q = q, N = N, SNR = SNR, rho = rho, s = s)
  rbind(cw, sup, prec, blpcv)
}

p <- 200
q <- 5
N <- 50
SNRs <- c(1.0, 10.0, 100.0) 
SNR <- 10
rho <- .35
s <- 10
sim(p, q, N, SNR, rho, s) # one trial

qq <- c(5)
NN <- c(100, 200)

num_trials <- 250
sparse_results <- lapply(NN, function(N){
  lapply(SNRs, function(snr){
  lapply(rep(qq, num_trials), function(q){sim(p, q, N, snr, .35, 10)})})
})
sparse_results <- Reduce(rbind, rlang::squash(sparse_results), init = list())
save(sparse_results, file = 'sparse_results.rda')

p2 <- sparse_results %>% 
  mutate(N = paste0('n = ', N),
         s = paste0('p = ', p, ', s = ', s)) %>%
  dplyr::select(Method, MSE_relative, SNR, N, s) %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(SNR), y = MSE_relative, color = Method)) + 
  geom_abline(aes(slope = 0, intercept = 1), linetype = 2) + 
  labs(x = 'Signal-to-Noise Ratio',
       y = 'Total Squared Error, Relative to LASSO') + 
  ylim(c(0, 2)) + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  facet_grid(N ~ s)

p2
ggsave('TSE_sparse_curds.pdf', plot = p2, width = 14, height = 10)