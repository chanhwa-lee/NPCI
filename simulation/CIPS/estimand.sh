#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH --output=estimand/estimand_rep-%a.out

m=10000 # Number of clusters in one simulation
r=100 # Number of binary vector sampling

Rscript "estimand.R" \
  -m $m \
  -r $r