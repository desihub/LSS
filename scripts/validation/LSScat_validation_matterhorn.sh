#!/bin/bash

#pipeline for all tracers starting from comb files; 16.3 hours total

set -e

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#export LSSCODE=$HOME ; do this before script, e.g., export LSSCODE=$HOME/LSScode for desica
#PYTHONPATH=$PYTHONPATH:$LSS/LSS/py

verspec=matterhorn-v2
survey=DA3
#LSS=$HOME/LSScode
#PYTHONPATH=$LSS/LSS/py:$PYTHONPATH
scriptdir=$HOME/LSScode/LSS/scripts
srun -N 1 -C cpu -t 02:0:00 --qos interactive --account desi python $scriptdir/validation/validation_sky.py --tracers all --version $1 --weight_col WEIGHT_SYS --survey $survey --verspec $verspec