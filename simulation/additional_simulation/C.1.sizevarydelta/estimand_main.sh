#!/bin/bash

# module add r/4.0.1

### Estimand Simulation ###
mkdir -p estimand
jid_estimand=$( sbatch --array=1-100 estimand.sh | awk '{print $4}' )