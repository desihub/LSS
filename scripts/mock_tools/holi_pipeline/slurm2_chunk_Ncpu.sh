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
ID0_ARRAY=$3

cd $LSS_DIR/scripts/mock_tools/holi_pipeline

# bash syntax : default value is 4 with syntax -4 ...
NTASKS=${SLURM_NTASKS:-2}
NCPU_PT=${SLURM_CPUS_PER_TASK:-4}
ID_ARRAY=${SLURM_ARRAY_TASK_ID:-0}
#ID0_ARRAY=${SLURM_ARRAY_TASK_MIN:-100}
SZ_CHUNK=$NTASKS
ID_FIRST=$((SZ_CHUNK*ID_ARRAY + ID0_ARRAY))

echo "Phase 1 : $NTASKS simulations in parallel with $NCPU_PT CPU, first ID ${ID_FIRST}"
# launch one process by chunk without --exclusive
srun -n $SZ_CHUNK -c $NCPU_PT ./slurm3_step134567.sh $LSS_DIR $DS_DIR $ID_FIRST


#echo "Phase 2 : Fiber assign with  CPUs"
#srun --exclusive -n $SZ_CHUNK -c $NCPU_PT ./step8_chunk.sh $LSS_DIR $DS_DIR $ID_FIRST $SZ_CHUNK $NCPU_PT 
#
# Remove this job, because runAltMTLRealizations 
#   * can't used all CPU allocated (use only 1 by seed , OpenMP not used ?)
#   * very long processing ~30 h (to be confirmed), limitation 48h can be exceded
# 
# => replace by sbatch with same number of task with 1 CPU
#

sbatch -n $NTASKS ./slurm3_step8.sh  $LSS_DIR $DS_DIR $ID_FIRST 
#echo "Phase 3 : clean holi pipeline"
#TODO
