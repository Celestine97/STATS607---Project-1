#! /bin/bash

for sigma in 10 20 30 40 50 60 70 80 90 100; do
  echo $sigma
  Rscript run_simulations_copas.R $sigma 0 0 '../simulation_results/copas_ridge_results/'
done;

