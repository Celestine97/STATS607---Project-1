#! /bin/bash

for sigma in 10 20 30 40 60 70 80 90 100; do
  sbatch \
  --export=sigma=${sigma} \
  run_ridge_simulations.sh
done;

