#!/bin/bash

# module add r/4.0.1

mkdir -p log
mkdir -p Rdata

### Estimator Simulation ###
jid_estimator=$( sbatch --array=1-45 estimator.sh | awk '{print $4}' )
