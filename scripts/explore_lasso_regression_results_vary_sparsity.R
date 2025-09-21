# script lasso varying sparsity

library(tidyverse)
library(MASS)
library(glmnet)
library(rlist)
library(superheat)

source('R/utils.R')

options(repr.plot.width=6, repr.plot.height=4) # plot sizes in this notebook

set.seed(5645654)
save_figs <- SAVE_FIG

results_dir <- 'results/simulation_results/joint_cv_lasso_results/'

file_list <- list.files(path = results_dir, pattern = '*.rds')

get_results_across_sparsity <- function(results_dir, name){
  file_list <- list.files(path = results_dir, pattern = '*.rds')
  
  result_vec <- c()
  sparsities <- c()
  for(i in 1:length(file_list)){
    result_list <- list.load(paste0(results_dir, file_list[i]))
    
    result_vec <- c(result_vec, result_list[[name]])
    # Extract sparsity from the saved results
    sparsities <- c(sparsities, rep(result_list$s, length(result_list[[name]])))
  }
  return(list(results = result_vec, 
              sparsities = sparsities))
}

slopes <- c()
sparsities <- c()
method <- c()

lambda_only_slopes <- get_results_across_sparsity(results_dir, 'lambda_only_slope_vs')
slopes <- c(slopes, lambda_only_slopes$results)
sparsities <- c(sparsities, lambda_only_slopes$sparsities)
method <- c(method, rep('lasso', length(lambda_only_slopes$sparsities)))

joint_slopes <- get_results_across_sparsity(results_dir, 'joint_slope_vs')
slopes <- c(slopes, joint_slopes$results)
sparsities <- c(sparsities, joint_slopes$sparsities)
method <- c(method, rep('joint', length(joint_slopes$sparsities)))

copas_result_dir <- 'results/simulation_results/copas_lasso_results/'
copas_slopes <- get_results_across_sparsity(copas_result_dir, 'copas_slope_vs')
slopes <- c(slopes, copas_slopes$results)
sparsities <- c(sparsities, copas_slopes$sparsities)
method <- c(method, rep('Copas', length(copas_slopes$sparsities)))

two_step_dir <- 'results/simulation_results/two_step_lasso_results/'
two_step_slopes <- get_results_across_sparsity(two_step_dir, 'two_step_slope_vs')
slopes <- c(slopes, two_step_slopes$results)
sparsities <- c(sparsities, two_step_slopes$sparsities)
method <- c(method, rep('two_step', length(two_step_slopes$sparsities)))

slopes_df <- data.frame(slopes = slopes, 
                        sparsities = sparsities, 
                        method = method)

slopes_df %>% 
  group_by(sparsities, method) %>% 
  summarize(median_slope = median(slopes, na.rm = TRUE)) %>%
  ggplot(aes(x = sparsities, y = median_slope)) + 
  geom_point(aes(color = method, shape = method), size = 4) + 
  geom_line(aes(color = method), linetype = 'dashed')  + 
  scale_color_manual(values = c('purple', 'green3', 'blue', 'orange')) + 
  geom_hline(yintercept = 1.0) + 
  ylab('median slope') + xlab('s') + theme_bw() + 
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))

if(save_figs){
  ggsave('results/figures/lasso_figures/slope_vs_sparsity.png', height = 6, width = 9)
}

mspes <- c()
sparsities <- c()
method <- c()

lambda_only_mspes <- get_results_across_sparsity(results_dir, 'lambda_only_mspe_vs')
mspes <- c(mspes, lambda_only_mspes$results)
sparsities <- c(sparsities, lambda_only_mspes$sparsities)
method <- c(method, rep('lasso', length(lambda_only_mspes$sparsities)))

joint_mspes <- get_results_across_sparsity(results_dir, 'joint_mspe_vs')
mspes <- c(mspes, joint_mspes$results)
sparsities <- c(sparsities, joint_mspes$sparsities)
method <- c(method, rep('joint', length(joint_mspes$sparsities))) 

copas_result_dir <- 'results/simulation_results/copas_lasso_results/'
copas_mspes <- get_results_across_sparsity(copas_result_dir, 'copas_mspe_vs')
mspes <- c(mspes, copas_mspes$results)
sparsities <- c(sparsities, copas_mspes$sparsities)
method <- c(method, rep('Copas', length(copas_mspes$sparsities)))

two_step_dir <- 'results/simulation_results/two_step_lasso_results/'
two_step_mspes <- get_results_across_sparsity(two_step_dir, 'two_step_mspe_vs')
mspes <- c(mspes, two_step_mspes$results)
sparsities <- c(sparsities, two_step_mspes$sparsities)
method <- c(method, rep('two_step', length(two_step_mspes$sparsities)))

mspes_df <- data.frame(mspe = mspes, 
                       sparsity = sparsities, 
                       method = method)

mspes_df %>% 
  group_by(sparsity, method) %>% 
  summarize(median_mspe = median(mspe, na.rm = TRUE)) %>% 
  spread(method, median_mspe) %>% 
  mutate(Copas = Copas / lasso, joint = joint / lasso, two_step = two_step / lasso)  %>% 
  gather(key = 'method', value = 'median_mspe', joint, lasso, two_step, Copas) %>% 
  filter(method != 'lasso')  %>% 
  
  filter(median_mspe < 3) %>%
  
  ggplot(aes(x = sparsity, y = median_mspe)) + 
  geom_point(aes(color = method, shape = method), size = 4) + 
  geom_line(aes(color = method), linetype = 'dashed')  + 
  scale_color_manual(values = c('purple', 'green3', 'orange')) + 
  geom_hline(yintercept = 1.0) + 
  ylab('MSPE relative to LASSO') + xlab('s') + theme_bw() + 
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))

if(save_figs){
  ggsave('results/figures/lasso_figures/mspe_vs_sparsity.png', height = 6, width = 9)    
}

