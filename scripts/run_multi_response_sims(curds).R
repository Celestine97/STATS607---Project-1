# curds - script

source('R/utils.R')
source('R/single_response.R')
source('R/multi_response.R')           # for multi-response scripts

set.seed(1234567891)

# Use test mode parameters
if(TEST_MODE) {
  pp <- c(50)           # Test with one dimension only
  NN <- c(100)          # Test with one sample size
  qq <- c(10, 15)       # Test with fewer response variables
  SNRs <- c(1.0)        # Test with one SNR
  num_trials <- 5       # Much fewer trials
} else { # all setting used in Breiman & Friedman (1997) simulations
  pp <- c(50, 75)
  NN <- c(100, 200)
  qq <- c(10, 15, 20, 25, 30)
  SNRs <- c(1.0, 3.0)
  num_trials <- 250
}
rhos <- seq(-.7, .7, .35) # not same as average corr

results <- lapply(pp, function(p){
  lapply(NN, function(N){
    lapply(rep(qq, num_trials), function(q){sim(p, q, N, 1.0, .35)}) 
  })}) 
results <- Reduce(rbind, rlang::squash(results), init = list())
# save(results, file = 'results.rda')
# Fix output paths
save(results, file = 'results/simulation_results/multi_response/multi_response_results(curds).rda')
library(ggplot2)

p1 <- results %>% 
  mutate(N = paste0('n = ', N),
         p = paste0('p = ', p)) %>%
  dplyr::select(Method, MSE_relative, q, N, p) %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(q), y = MSE_relative, color = Method)) + 
  geom_abline(aes(slope = 0, intercept = 1), linetype = 2) + 
  labs(x = 'Number of Response Variables',
       y = 'Total Squared Error, Relative to OLS') + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  facet_grid(N ~ p)

# p1 
ggsave('results/figures/TSE_curds.pdf', plot = p1, width = 14, height = 10)