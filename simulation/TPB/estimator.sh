#!/bin/bash
#SBATCH -t 40:00:00
#SBATCH --output=log/estimator_rep-%a.out

n=$1 # Number of clusters per simulation

Rscript "../../estimator.R" \
  -n $n