#!/bin/bash

# module add r/4.0.1

m=$1 # Number of clusters per each simulation
r=$2 # Number of subsampling approximation

d=data/m${m}_r${r}

echo $d
mkdir -p $d
mkdir -p $d/log
mkdir -p $d/Rdata
cd $d

### Estimand Simulation ###
jid_estimand=$( sbatch --array=1-100 ../../estimand.sh | awk '{print $4}' )

### Estimator Simulation ###
jid_estimator=$( sbatch --array=1-1000 ../../estimator.sh $m $r | awk '{print $4}' )
