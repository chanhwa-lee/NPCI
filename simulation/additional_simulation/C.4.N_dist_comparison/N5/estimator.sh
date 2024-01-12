#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --output=log/estimator_rep-%a.out

m=$1 # Number of clusters per simulation
r=$2 # Number of sample splitting repetition

Rscript "../../estimator.R" \
  -m $m \
  -r $r