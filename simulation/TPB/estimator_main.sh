#!/bin/bash

n=$1 # Number of clusters per each simulation

d=data/n${n}

echo $d
mkdir -p $d
mkdir -p $d/log
mkdir -p $d/Rdata
cd $d

### Estimator Simulation ###
jid_estimator=$( sbatch --array=1-1000 ../../estimator.sh $n | awk '{print $4}' )