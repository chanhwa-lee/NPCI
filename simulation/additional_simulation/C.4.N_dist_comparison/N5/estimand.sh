#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --output=log/estimand_rep-%a.out

m=100000 # Number of clusters in one simulation
r=100 # Number of binary vector sampling

Rscript "../../estimand.R" \
  -m $m \
  -r $r