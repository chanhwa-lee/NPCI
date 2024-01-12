#!/bin/bash

### Estimator Simulation ###
m=$1 # Number of clusters per each simulation
r=$2 # Number of subsampling approximation

d=data/m${m}_r${r}

echo $d
mkdir -p $d
mkdir -p $d/log
mkdir -p $d/Rdata
cd $d

jid_estimator=$( sbatch --array=1-2000 ../../estimator.sh $m $r | awk '{print $4}' )
