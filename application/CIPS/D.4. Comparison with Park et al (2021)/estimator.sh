#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mem=32000
#SBATCH --output=log/estimator_rep-%a.out

# module add r/4.0.1
Rscript "estimator.R"
