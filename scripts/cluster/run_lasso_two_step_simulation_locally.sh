#! /bin/bash

for s in 980 950 920 890 860 830; do
  echo $s
  Rscript run_simulations_two_step.R 5 $s 1 '../simulation_results/two_step_lasso_results/varying_sparsity/'
done;

