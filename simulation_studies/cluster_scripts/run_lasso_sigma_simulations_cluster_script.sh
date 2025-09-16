#! /bin/bash

for sigma in 11 12 13 14 15 ; do
  sbatch \
  --export=sigma=${sigma},s=950\
  run_lasso_simulations.sh
done;

