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
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $scriptdir/validation/validation_skyclus.py --tracers all --version $1 --weight_col WEIGHT_SYS --survey $survey --verspec $verspec

srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $scriptdir/validation/validation_focal.py --survey $survey --verspec $verspec --version $1 

srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $scriptdir/validation/validation_improp_clus.py --survey $survey --verspec $verspec --version $1 

srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $scriptdir/validation/validation_tsnr_zbin.py --survey $survey --verspec $verspec --version $1 