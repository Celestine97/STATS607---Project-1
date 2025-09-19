#! /bin/bash

for s in 980 ; do
  echo $s
  Rscript run_simulations_copas.R 5 $s 1 '../simulation_results/copas_lasso_results/varying_sparsity/'
done;

