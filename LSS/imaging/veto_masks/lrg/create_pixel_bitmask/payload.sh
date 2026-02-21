#!/bin/bash

module load parallel
source /project/projectdirs/desi/software/desi_environment.sh 20.7

if [[ -z "${SLURM_NODEID}" ]]; then
    echo "need \$SLURM_NODEID set"
    exit
fi
if [[ -z "${SLURM_NNODES}" ]]; then
    echo "need \$SLURM_NNODES set"
    exit
fi

# cat $1 |                                               \
# awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" \
# 'NR % NNODE == NODEID' |                               \
# parallel --jobs 1 python lrg_gaiamask.py {}

# cat $1 |                                               \
# awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" \
# 'NR % NNODE == NODEID' |                               \
# parallel --jobs 1 python lrg_wisemask.py {}

cat $1 |                                               \
awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" \
'NR % NNODE == NODEID' |                               \
parallel --jobs 1 python combine_lrg_masks.py {}
