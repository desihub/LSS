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

NCPU=${SLURM_CPUS_PER_TASK:10}
ID_ARRAY=${SLURM_ARRAY_TASK_ID:0}

# fix 1 CPU by simu
SIZE_CHUNK=$NCPU

OFFSET=$((SIZE_CHUNK*ID_ARRAY))

echo "Phase 1 : $NCPU simulations in parallel, first ID ${OFFSET}"
srun -n "$NCPU" -c1 ./slurm_step134567.sh $LSS_DIR $DS_DIR $OFFSET

echo "Phase 2 : Fiber assign with $NCPU CPUs"
srun -n1 -c"$N" step_fa.sh ./fa_chunk_${OFFSET}.log  $2>&1

echo "Phase 3 : clean holi pipeline"
#TODO