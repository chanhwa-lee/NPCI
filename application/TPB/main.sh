#!/bin/bash

# module add r/4.0.1

mkdir -p log
mkdir -p result

### Estimator Simulation ###
sbatch -t 2- --mem=32000 --output=log/estimator_rep-%a.out --array=1-50 --wrap="Rscript Estimation.R"