library(tidyverse)
library(glmnet)
source('R/utils.R')
source('R/single_response.R')
source('R/multi_response.R')           # for multi-response scripts
source('R/multi_response_sparse.R')    # for sparse scripts
# source('../simulation_studies/simulation_utils.R')

set.seed(54235245)
load('./fmri_data/fMRIdata.RData')
ls()

dim(fit_feat)
dim(resp_dat)

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

# save(X_train, y_train, X_test, y_test, file = './fmri_train_test_split.RData')

load_cv_fits <- TRUE
# get lasso results
if(load_cv_fits == FALSE){
  cv_fits <- list()
  for(i in 1:20){
    cat('lasso on voxel i: ', i, '\n')
    cv_fits[[i]] <- cv.glmnet(X_train, y_train[, i])
  }
  
  save(cv_fits, file = './fmri_results/lasso_cv_fits.RData')
}else{
  load('./fmri_results/lasso_cv_fits.RData')
}

# get cv predictions
cv_predictions <- matrix(0, ncol = 20, nrow = n_test)
cv_slopes <- rep(0, 20)
cv_mspes <- rep(0, 20)
cv_correlations <- rep(0, 20)
for(i in 1:20){
  cat('predicting on voxel i: ', i, '\n')
  cv_predictions[, i] <- predict(cv_fits[[i]], newx = X_test, 
                                 s = cv_fits[[i]]$lambda.min)
  
  cv_correlations[i] <- cor(cv_predictions[, i], y_test[, i])
  cv_slopes[i] <- drop(get_slope(cv_predictions[, i], y_test[, i]))
  cv_mspes[i] <- sqrt(sum((cv_predictions[, i] - y_test[, i])**2))
}

# get joint lambda k fits
load_joint_fits <- TRUE
if(load_joint_fits == FALSE){
  joint_results_mat <- matrix(0, nrow = 2, ncol = n_voxels)
  rownames(joint_results_mat) <- c('lambda', 'k')
  
  k_range <- seq(0.8, 2.0, by = 0.05)
  
  for(i in 1:n_voxels){
    cat('joint lambda k selection on voxel i: ', i, '\n')
    lambdas <- cv_fits[[i]]$lambda
  
    joint_fit <- 
      choose_joint_lambda_k(X_train, y_train[, i], lambdas, k_range, 
                            nfolds = 5)
    
    joint_results_mat[1, i] <- joint_fit$lambda
    joint_results_mat[2, i] <- joint_fit$k
  }
  
  cat('saving joint results mat')
  save(joint_results_mat, file = './fmri_results/joint_results2.RData')
}else(
  load('./fmri_results/joint_results2.RData')
)


# get joint predictions
joint_predictions <- matrix(0, ncol = 20, nrow = n_test)
joint_slopes <- rep(0, 20)
joint_mspes <- rep(0, 20)
joint_correlations <- rep(0, 20)
for(i in 1:n_voxels){
  cat('predicting on voxel i: ', i, '\n')
  joint_predictions[, i] <- predict(cv_fits[[i]], newx = X_test, 
                                 s = joint_results_mat['lambda', i]) * 
                              joint_results_mat['k', i]
  
  joint_correlations[i] <- cor(joint_predictions[, i], y_test[, i])
  joint_slopes[i] <- drop(get_slope(joint_predictions[, i], y_test[, i]))
  joint_mspes[i] <- sqrt(sum((joint_predictions[, i] - y_test[, i])**2))
}


# data.frame(lambda_only = cv_slopes,
#               joint = joint_slopes,
#            voxel = seq(1, 20, by = 1)) %>%
#   gather(key = 'method', value = 'slope', lambda_only, joint) %>%
#   ggplot() + geom_point(aes(x = voxel, y = slope, color = method)) +
#   geom_hline(yintercept = 1.0) +
#   scale_color_manual(values = c('green3', 'blue')) +
#   ggtitle('test set slopes for the twenty voxels') +
#   theme_bw() +
#   theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))
data.frame(lambda_only = cv_slopes,
           joint = joint_slopes,
           voxel = seq(1, 20, by = 1)) %>%
  gather(key = 'method', value = 'slope', lambda_only, joint) %>%
  filter(voxel < 20) %>%
  ggplot() + 
  geom_bar(aes(voxel, slope - 1, color = method, fill = method), stat = "identity", 
           position = "dodge", width = 0.5, alpha = 0.5) +
  scale_y_continuous(labels = function(x) {x + 1}) + 
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = c('green3', 'blue')) + 
  scale_fill_manual(values = c('green3', 'blue')) +
  ylab('slope') + ggtitle('Test set slopes each voxel') +
  theme_bw() +
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))

ggsave('../writing/single_response_fmri_figures/fmri_slopes.png', 
       width = 9, height = 6)

data.frame(lambda_only = cv_mspes,
           joint = joint_mspes,
           voxel = seq(1, 20, by = 1)) %>%
  gather(key = 'method', value = 'mspe', lambda_only, joint) %>%
  filter(voxel < 20) %>%
  ggplot() + 
  geom_bar(aes(voxel, mspe, color = method, fill = method), stat = "identity", 
           position = "dodge", width = 0.5, alpha = 0.5) +
  scale_color_manual(values = c('green3', 'blue')) + 
  scale_fill_manual(values = c('green3', 'blue')) +
  ylab('MSPE') + ggtitle('Test set MSPEs for each voxel') +
  theme_bw() +
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))
ggsave('../writing/single_response_fmri_figures/fmri_mspe.png', 
       width = 9, height = 6)


data.frame(lasso = cv_correlations,
           joint = joint_correlations,
           voxel = seq(1, 20, by = 1)) %>%
  gather(key = 'method', value = 'correlations', lasso, joint) %>%
  filter(voxel < 20) %>%
  ggplot() + 
  geom_bar(aes(voxel, correlations, color = method, fill = method), stat = "identity", 
           position = "dodge", width = 0.5, alpha = 0.8) +
  scale_color_manual(values = c('green3', 'blue')) + 
  scale_fill_manual(values = c('green3', 'blue')) +
  ylab('correlations') + ggtitle('Test set correlations for each voxel') +
  theme_bw() +
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))

ggsave('../writing/single_response_fmri_figures/fmri_correlations.png',
       width = 12, height = 6)

voxel <- 15
data.frame(y_obs = y_test[, voxel] - mean(y_test[, voxel]), 
           lambda_only = cv_predictions[, voxel] - mean(cv_predictions[, voxel]), 
           joint = joint_predictions[, voxel] - mean(joint_predictions[, voxel])) %>%
  gather(key = 'procedure', value = 'y_pred', lambda_only, joint) %>%
  ggplot() + geom_point(aes(x = y_pred, y = y_obs, colour = procedure), alpha = 0.2) + 
  stat_smooth(aes(x = y_pred, y = y_obs, colour = procedure), 
              linetype = 'solid', method = 'lm', se = FALSE, fullrange = TRUE) + 
  geom_abline(slope = 1) + scale_color_manual(values=c("green3", "blue")) + 
  xlab('Predicted y') + ylab('Observed y') + 
  theme_bw() + 
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))
ggsave('../writing/single_response_fmri_figures/obs_y_vs_pred.png', width = 9, height = 6)
