library(tidyverse)
library(dplyr)
library(MASS)
library(Matrix)
library(expm)
library(quadprog)

## generate multi-response linear model (x, y)  
## following the Design section of Breiman-Friedman (1997)
generate_data <- function(p, q, N, SNR, rho, Sigma = 'eye'){
  # generate covariates x from 
  # AR(1) Gaussian graphical model
  r <- 2*runif(1) - 1
  V <- r**abs(outer(1:p, 1:p, "-"))
  x <- mvrnorm(N, (1:p)*0, V) 
  
  # generate true coefficients A
  K <- 10
  C <- mvrnorm(K, (1:q)*0, rho**abs(outer(1:q, 1:q, "-"))) 
  j <- sample(seq(p), K, replace = TRUE)
  l <- sample(seq(6), K, replace = TRUE)
  G <- pmax(Matrix(rep(l, p), nrow = p, byrow = TRUE) - abs(outer(1:p, j, "-")), 0)^2
  G <- G / Matrix(rep(colSums(G), p), nrow = p, byrow = TRUE)
  A <- G %*% C
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

## Run OLS separately on each problem
ols_separate <- function(x, y){
  q <- dim(y)[2]
  sapply(1:q, function(i){
    as.vector(lm(y[, i] ~ x + 0)$coefficients)}) %>% t
}

## Shrink multi-response predictions according to model-
## based closed form: "does not provide enough shrinkage"
curds_n_whey <- function(x, y, A_ols){
  can_cor <- cancor(x, y)
  T_ <- can_cor$ycoef
  C2 <- (can_cor$cor)^2
  r <- dim(x)[2] / dim(x)[1]
  D_ <- diag(C2 / (C2 + r*(1-C2)))
  B <- solve(T_) %*% D_ %*% T_
  q <- dim(y)[2]
  
  B %*% A_ols
}

## Shrink multi-response predictions according to  
## "generalized CV" - selecting diagonal matrix w/ CV
curds_n_whey_gcv <- function(x, y, A_ols, nfolds = 5){
  N <- dim(y)[1] 
  q <- dim(y)[2]
  p <- dim(x)[2]
  fold <- sort(rep(1:nfolds, ceiling(N/nfolds)))[1:N]
  
  y_hat_ols_cv <- Reduce(rbind,
                         lapply(1:nfolds, function(hold){
                           train_x <- x[fold != hold, ]
                           train_y <- y[fold != hold, ]
                           valid_x <- x[fold == hold, ]
                           A_ols_cv <- ols_separate(train_x, train_y)
                           valid_x %*% t(A_ols_cv)
                         }), NULL)
  
  T_cv <- lapply(1:nfolds, function(hold){
    train_x <- x[fold != hold, ]
    train_y <- y[fold != hold, ]
    cancor(train_x, train_y)$ycoef
  })
  
  T_inv_cv <- Map(solve, T_cv)
  R_cv <- sapply(1:N, function(n){
    T_cv[[fold[n]]] %*% y_hat_ols_cv[n, ]
  })
  
  dvec <- sapply(1:N, function(n){
    (y[n, ] %*% T_inv_cv[[fold[n]]]) * R_cv[, n]
  }) %>% rowSums
  
  Dmat <- Reduce('+', 
                 lapply(1:N, function(n){
                   (t(T_inv_cv[[fold[n]]]) %*% T_inv_cv[[fold[n]]]) * outer(R_cv[, n], R_cv[, n])
                 }))
  
  D_ <- diag(pmax(solve(Dmat) %*% dvec, 0) %>% as.vector())
  
  can_cor <- cancor(x, y)
  T_ <- can_cor$ycoef
  B <- solve(T_) %*% D_ %*% T_
  
  B %*% A_ols
}

## Generate data, estimate coefficients, and calculate
## metrics used by Breiman & Friedman (1997)
sim <- function(p, q, N, SNR, rho){
  data <- generate_data(p, q, N, SNR, rho) 
  A <- data$A
  V <- data$V
  x <- data$x
  y <- data$y
  
  A_ols <- ols_separate(x, y)
  A_cw <- curds_n_whey(x, y, A_ols)
  A_cw_cv <- curds_n_whey_gcv(x, y, A_ols)
  
  S <- (1/N) * t(x) %*% x
  sigma_hat_sq <- sum((y - x %*% t(A_ols))**2) / ((N-p)*q)
  nu <- N - p - 1
  D_copas <- diag(1 - ((p-2) * (sigma_hat_sq / N) * nu) / 
                    ((nu+2) * diag(A_ols %*% S %*% t(A_ols))))
  A_copas <- D_copas %*% A_ols
  
  MSE_ols <- diag((A - A_ols) %*% V %*% t(A - A_ols))
  MSE_cw  <- diag((A - A_cw) %*% V %*% t(A - A_cw))
  MSE_copas <- diag((A - A_copas) %*% V %*% t(A - A_copas))
  MSE_cw_cv <- diag((A - A_cw_cv) %*% V %*% t(A - A_cw_cv))
  
  cw <- data.frame(MSE_relative = c(sum(MSE_cw) / sum(MSE_ols)), 
                   MSE_individual = c(mean(MSE_cw / MSE_ols)), 
                   MSE_worst = c(max(MSE_cw / MSE_ols)), 
                   Method = c('Curds & Whey (Closed Form)'), 
                   p = p, q = q, N = N, SNR = SNR, rho = rho)
  
  copas <- data.frame(MSE_relative = c(sum(MSE_copas) / sum(MSE_ols)), 
                      MSE_individual = c(mean(MSE_copas / MSE_ols)), 
                      MSE_worst = c(max(MSE_copas / MSE_ols)), 
                      Method = c('Separate Copas Shrinkage'), 
                      p = p, q = q, N = N, SNR = SNR, rho = rho)
  
  cw_cv <- data.frame(MSE_relative = c(sum(MSE_cw_cv) / sum(MSE_ols)), 
                      MSE_individual = c(mean(MSE_cw_cv / MSE_ols)), 
                      MSE_worst = c(max(MSE_cw_cv / MSE_ols)), 
                      Method = c('Curds & Whey (GCV)'), 
                      p = p, q = q, N = N, SNR = SNR, rho = rho)
  rbind(cw, copas, cw_cv)
}