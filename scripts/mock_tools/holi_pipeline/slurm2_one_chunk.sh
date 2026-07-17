#!/usr/bin/env bash

# Description:
# ===========
# Holi pipeline for N mock used by fiber assignement
# Phase 1: N mock simu (step134567) in parallel with 1 CPU each
# Phase 2: fiber assgnment with N CPU
# Phase 3: clean

set -e

LSS_DIR=$1 
DS_DIR=$2

# bash syntax : default value is 4 with syntax -4 ...
NCPU=${SLURM_CPUS_PER_TASK:-4}
ID_ARRAY=${SLURM_ARRAY_TASK_ID:-14}

# fix 1 CPU by simu
SIZE_CHUNK=$NCPU
OFFSET=$((SIZE_CHUNK*ID_ARRAY))

echo "Phase 1 : $NCPU simulations in parallel, first ID ${OFFSET}"
# launch N tasks in parallel
srun -n $NCPU -c1 ./slurm3_step134567.sh $LSS_DIR $DS_DIR $OFFSET

echo "Phase 2 : Fiber assign with $NCPU CPUs"
srun -n 1 -c $NCPU step_fa.sh $LSS_DIR $DS_DIR $OFFSET $NCPU > $DS_DIR/fa_chunk_${OFFSET}.log  2>&1

echo "Phase 3 : clean holi pipeline"
#TODO