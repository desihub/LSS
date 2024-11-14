#!/bin/bash

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK



num=$(($1 + 1)) 

PROG=$2

if [ "$PROG" == "DARK" ]; then
	python $LSSCODE/scripts/mock_tools/prepare_mocks_Y3_dark.py --mockver ab_secondgen --mockpath /global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky_v4_1 --realmin $1 --realmax $num --isProduction y --split_snapshot y --new_version AbacusSummit_v4_1 --prog dark

elif [ "$PROG" == "BRIGHT" ]; then
	python $LSSCODE/scripts/mock_tools/prepare_mocks_Y1_bright.py  --mockver ab_secondgen_cosmosim --realmin $1 --realmax $num --prog bright --isProduction y --rbandcut 19.5
    
else
    echo "Need to define DARK or BRIGHT"
fi

