#!/usr/bin/env bash

# Description:
# ===========
# slurm adaptor for step134567_one.sh

#LSS_DIR=$1
#DS_DIR=$2     # root directory of mock
#ID_FIRST=$3     # first id mock to process

PROCID=${SLURM_PROCID:-0}
# ID mock to process in task rank
SIM_ID=$((PROCID + $3))

CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}

if [[ "$CPUS_PER_TASK" =~ ^[0-9]+$ ]] && (( CPUS_PER_TASK < 4 )); then
	time ./step134567_one.sh $1 $2 ${SIM_ID} > $2/step134567_seed${SIM_ID}.log 2>&1
else
    echo "**Used step134567_oneNcpu.sh**"
	time ./step134567_oneNcpu.sh $1 $2 ${SIM_ID} > $2/step134567_seed${SIM_ID}.log 2>&1
fi
