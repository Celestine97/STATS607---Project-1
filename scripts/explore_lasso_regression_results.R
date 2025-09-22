# script explore lasso regression results
library(tidyverse)
library(MASS)
library(glmnet)
library(rlist)
library(superheat)

source('R/utils.R')

options(repr.plot.width=6, repr.plot.height=4) # plot sizes in this notebook

set.seed(5645654)
save_figs <- SAVE_FIG

# Use config values for directories
results_dir <- LASSO_OUTPUT_DIR
file_list <- list.files(path = results_dir, pattern = '*.rds')

get_results_across_sigma <- function(results_dir, name){
  file_list <- list.files(path = results_dir, pattern = '*.rds')
  
  result_vec <- c()
  sigmas <- c()
  for(i in 1:length(file_list)){
    result_list <- list.load(paste0(results_dir, file_list[i]))
    
    result_vec <- c(result_vec, result_list[[name]])
    sigmas <- c(sigmas, rep(result_list$sigma, length(result_list[[name]])))
    
  }
  return(list(results = result_vec, 
              sigmas = sigmas))
}

slopes <- c()
sigmas <- c()
method <- c()

lambda_only_slopes <- get_results_across_sigma(results_dir, 'lambda_only_slope_vs')
slopes <- c(slopes, lambda_only_slopes$results)
sigmas <- c(sigmas, lambda_only_slopes$sigmas)
method <- c(method, rep('lasso', length(lambda_only_slopes$sigmas)))

joint_slopes <- get_results_across_sigma(results_dir, 'joint_slope_vs')
slopes <- c(slopes, joint_slopes$results)
sigmas <- c(sigmas, joint_slopes$sigmas)
method <- c(method, rep('joint', length(joint_slopes$sigmas)))  # Use joint_slopes length

# Use config values for Copas directory
copas_result_dir <- LASSO_COPAS_OUTPUT_DIR
copas_slopes <- get_results_across_sigma(copas_result_dir, 'copas_slope_vs')
slopes <- c(slopes, copas_slopes$results)
sigmas <- c(sigmas, copas_slopes$sigmas)
method <- c(method, rep('Copas', length(copas_slopes$sigmas)))

# Use config values for two-step directory
two_step_dir <- LASSO_TWO_STEP_OUTPUT_DIR
two_step_slopes <- get_results_across_sigma(two_step_dir, 'two_step_slope_vs')
slopes <- c(slopes, two_step_slopes$results)
sigmas <- c(sigmas, two_step_slopes$sigmas)
method <- c(method, rep('two_step', length(two_step_slopes$sigmas)))

slopes_df <- data.frame(slopes = slopes, 
                        sigma = sigmas, 
                        method = method)

slopes_df %>% 
  filter(sigma < 7, slopes < 2) %>%
  ggplot() + 
  geom_boxplot(aes(x = factor(sigma), y = slopes, colour = method), position=position_dodge(width=0.8))

slopes_df %>% 
  filter(sigma < 7, slopes < 2) %>%
  group_by(sigma, method) %>% 
  summarize(median_slope = median(slopes)) %>%
  ggplot(aes(x = sigma, y = median_slope)) + 
  geom_point(aes(color = method, shape = method), size = 4) + 
  geom_line(aes(color = method), linetype = 'dashed')  + 
  scale_color_manual(values = c('purple', 'green3', 'blue', 'orange')) + 
  geom_hline(yintercept = 1.0) + 
  ylab('median slope') + theme_bw() + 
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))

if(save_figs){
  ggsave(paste0(LASSO_FIGURES_DIR, 'slope_vs_sigma.png'), height = 6, width = 9)
}

# Lets look at one particular sigma

slopes_df %>% 
  filter(sigma == 5, slopes < 2) %>%
  ggplot() + 
  geom_boxplot(aes(x = sigma, y = slopes, colour = method, fill = method), alpha = 0.3, 
               position=position_dodge(width=0.8), width = 0.5) + 
  scale_color_manual(values = c('purple', 'green3', 'blue', 'orange')) + 
  scale_fill_manual(values = c('purple', 'green3', 'blue', 'orange')) + 
  geom_hline(yintercept = 1.0, color = 'red', linetype = 'dashed') + 
  ylab('slope') + 
  theme_bw() + 
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + scale_x_continuous(breaks = NULL)

if(save_figs){
  ggsave(paste0(LASSO_FIGURES_DIR, 'slopes_sigma5.png'), height = 6, width = 9)    
}

slopes_df %>% 
  filter(sigma == 5, method == 'lasso') %>%
  ggplot() + 
  geom_histogram(aes(x = slopes), colour = 'blue', fill = 'blue', alpha = 0.3) + 
  ggtitle('Slopes when running LASSO') + 
  geom_vline(xintercept = 1.0, linetype = 'dashed') + 
  theme_bw() + 
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5)) 

if(save_figs){
  ggsave(paste0(LASSO_FIGURES_DIR, 'hist_lasso_sigma5.png'), height = 6, width = 9)    
}

# Look at MSPEs
mspes <- c()
sigmas <- c()
method <- c()

lambda_only_mspes <- get_results_across_sigma(results_dir, 'lambda_only_mspe_vs')
mspes <- c(mspes, lambda_only_mspes$results)
sigmas <- c(sigmas, lambda_only_mspes$sigmas)
method <- c(method, rep('lasso', length(lambda_only_mspes$sigmas)))

joint_mspes <- get_results_across_sigma(results_dir, 'joint_mspe_vs')
mspes <- c(mspes, joint_mspes$results)
sigmas <- c(sigmas, joint_mspes$sigmas)
method <- c(method, rep('joint', length(joint_mspes$sigmas)))  # <- Change this line

# Use config values for Copas directory
copas_result_dir <- LASSO_COPAS_OUTPUT_DIR
copas_mspes <- get_results_across_sigma(copas_result_dir, 'copas_mspe_vs')
mspes <- c(mspes, copas_mspes$results)
sigmas <- c(sigmas, copas_mspes$sigmas)
method <- c(method, rep('Copas', length(copas_mspes$sigmas)))

# Use config values for two-step directory
two_step_dir <- LASSO_TWO_STEP_OUTPUT_DIR
two_step_mspes <- get_results_across_sigma(two_step_dir, 'two_step_mspe_vs')
mspes <- c(mspes, two_step_mspes$results)
sigmas <- c(sigmas, two_step_mspes$sigmas)
method <- c(method, rep('two_step', length(two_step_mspes$sigmas)))

mspes_df <- data.frame(mspe = mspes, 
                       sigma = sigmas, 
                       method = method)

mspes_df %>% 
  filter(sigma < 7) %>% 
  group_by(sigma, method) %>% 
  summarize(median_mspe = median(mspe, na.rm = TRUE)) %>% 
  spread(method, median_mspe) %>% 
  mutate(Copas = Copas / lasso, joint = joint / lasso, two_step = two_step / lasso) %>% 
  gather(key = 'method', value = 'median_mspe', joint, lasso, two_step, Copas) %>% 
  filter(method != 'lasso')  %>% 
  filter(median_mspe < 2.5) %>%
  ggplot(aes(x = sigma, y = median_mspe)) + 
  geom_point(aes(color = method, shape = method), size = 4) + 
  geom_line(aes(color = method), linetype = 'dashed')  + 
  scale_color_manual(values = c('purple', 'green3', 'orange')) + 
  geom_hline(yintercept = 1.0) + 
  ylab('MSPE relative to LASSO') + theme_bw() + 
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5))
if(save_figs){
  ggsave(paste0(LASSO_FIGURES_DIR, 'mspe_vs_sigma.png'), height = 6, width = 9)    
}

mspes_df %>% 
  filter(sigma == 5) %>%
  ggplot() + 
  geom_boxplot(aes(x = sigma, y = mspe, colour = method, fill = method), alpha = 0.3, 
               position=position_dodge(width=0.8), width = 0.5) + 
  scale_color_manual(values = c('purple', 'green3', 'blue', 'orange')) + 
  scale_fill_manual(values = c('purple', 'green3', 'blue', 'orange')) + 
  ylab('MSPE') + 
  theme_bw() + 
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + scale_x_continuous(breaks = NULL)

if(save_figs){
  ggsave(paste0(LASSO_FIGURES_DIR, 'mspes_sigma5.png'), height = 6, width = 9)    
}