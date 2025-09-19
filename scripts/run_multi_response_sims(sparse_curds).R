# sparse curds - script

source('R/utils.R')
source('R/single_response.R')
source('R/multi_response.R')           # for multi-response scripts
source('R/multi_response_sparse.R')    # for sparse scripts
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