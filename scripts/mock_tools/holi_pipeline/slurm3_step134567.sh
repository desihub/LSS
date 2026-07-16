#!/usr/bin/env bash

# Description:
# ===========
# slurm adaptor for step134567_one.sh

#LSS_DIR=$1
#DS_DIR=$2     # root directory of mock
#OFFSET=$3     # offset id simu

SIM_ID=$((SLURM_PROCID + $3))

time ./step134567_one.sh $1 $2 ${SIM_ID} step134567_${SIM_ID}.log  $2>&1
