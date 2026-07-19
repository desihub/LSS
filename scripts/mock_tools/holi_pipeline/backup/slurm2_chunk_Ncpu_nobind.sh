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
NTASKS=${SLURM_NTASKS:-2}
NCPU_PT=${SLURM_CPUS_PER_TASK:-4}
ID_ARRAY=${SLURM_ARRAY_TASK_ID:-0}
ALL_CPU=$((NTASKS*NCPU_PT))
SZ_CHUNK=$NTASKS
ID_FIRST=$((SZ_CHUNK*ID_ARRAY))

# Ensure nested sruns do not inherit strict binding masks.
export SLURM_CPU_BIND=none
unset SLURM_CPU_BIND_LIST

echo "Phase 1 : $NTASKS simulations in parallel with $NCPU_PT CPU, first ID ${ID_FIRST}"
# launch one process by chunk
srun --cpu-bind=none --export=ALL,SLURM_CPU_BIND=none -n $SZ_CHUNK -c $NCPU_PT ./slurm3_step134567.sh $LSS_DIR $DS_DIR $ID_FIRST

echo "Phase 2 : Fiber assign with $ALL_CPU CPUs"
srun --cpu-bind=none --export=ALL,SLURM_CPU_BIND=none -n 1 -c $ALL_CPU step8_chunk.sh $LSS_DIR $DS_DIR $ID_FIRST $SZ_CHUNK $NCPU_PT > $DS_DIR/fa_chunk_${ID_FIRST}.log  2>&1

echo "Phase 3 : clean holi pipeline"
#TODO