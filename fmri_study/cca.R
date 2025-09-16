library(tidyverse)
library(reshape2)
library(plotROC)
library(plot3D)
library(corrplot)
#library(CCA)
library(PMA)
library(glmnet)

load('outputs/lasso_fit_final.RData')


# load the data
load("data/fMRIdataWithSplits.RData")

feats_keep <- (apply(X_train, 2, sd) > 0)
X <- X_train[, feats_keep]
y <- y_train

out <- CCA(X, y, K = 20, penaltyz = 1, niter = 200) # vpos = TRUE, 

B <- function(V, cors, r){
  D_ <- diag(cors^2 / (cors^2 + r*(1-cors)^2))
  solve(V) %*% D_ %*% V
}

p <- dim(X)[2] 
n <- dim(X)[1]
r_ <- p/n
B_bf <- B(out$v, out$cors, p/n)
B_bf <- pmax(B_bf, 0)
B_bf <- B_bf / matrix(rep(rowSums(B_bf), 20), nrow = 20)

image(B_bf, axes=FALSE, zlim=c(0,1), 
      col = gray.colors(100, start = 1, end = 0, gamma = 1))




# how much to shrink B?
y_hat_test <- sapply(seq(20), 
  function(v){
    lasso_fit_final[[v]] %>% # final model
    predict(newx = X_test) %>% 
    as.vector()})

mse <- function(u, v){mean((u-v)**2)}

rs <- seq(-10,30,.1)
cvpmse <- sapply(rs, function(r){
  shrunk_preds <- t(B(out$v, out$cors, r) %*% t(y_hat_test))
  sapply(seq(20), function(v){mse(shrunk_preds[, v], y_test[, v])}) %>% mean
  })

cvcor <- sapply(rs, function(r){
  shrunk_preds <- t(B(out$v, out$cors, r) %*% t(y_hat_test))
  sapply(seq(20), function(v){cor(shrunk_preds[, v], y_test[, v])}) %>% mean
})


p1 <- ggplot() + 
  geom_line(aes(x = rs, y = log(1+cvpmse))) +
  theme_minimal() +
  labs(x = 'Shrinkage r', 
       y = 'Log of Prediction MSE', 
       title = 'MSE on Test Data vs Shrinkage Factor') +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot() + 
  geom_line(aes(x = rs, y = cvcor)) +
  theme_minimal() +
  labs(x = 'Shrinkage r', 
       y = 'Log of Prediction Correlation', 
       title = 'Cor on Test Data vs Shrinkage Factor') +
  theme(plot.title = element_text(hjust = 0.5))









