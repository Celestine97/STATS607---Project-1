#! /bin/bash

for sparsity in 920; do
  sbatch \
  --export=sigma=5,s=${sparsity}\
  run_lasso_simulations.sh
done;

