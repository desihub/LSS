#!/bin/bash




export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

first_rea=0
last_rea=25

cpus_needed=$((last_rea - first_rea))

#num=$(($1 + 1)) 

PROG=$1
echo $PROG

if [ "$PROG" == "DARK" ]; then
	srun -n 1 -C cpu -t 00:40:00 --cpus-per-task $cpus_needed --qos interactive --account desi python $LSSCODE/scripts/mock_tools/prepare_mocks_Y3_dark.py --mockver ab_secondgen --mockpath /global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky_v4_1 --realmin $first_rea --realmax $last_rea --isProduction y --split_snapshot y --new_version AbacusSummit_v4_1 --prog dark

elif [ "$PROG" == "BRIGHT" ]; then
	srun -n 1 -C cpu -t 00:40:00 --cpus-per-task $cpus_needed --qos interactive --account desi python $LSSCODE/scripts/mock_tools/prepare_mocks_Y3_bright.py  --mockver ab_secondgen_cosmosim --realmin $1 --realmax $num --prog bright --isProduction y --rbandcut 19.5
    
else
    echo "Need to define DARK or BRIGHT"
fi

