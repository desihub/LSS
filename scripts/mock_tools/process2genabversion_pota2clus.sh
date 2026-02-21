#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

for (( i=$1;i<=$2;i++ ))
do
 srun -N 1 -C cpu -t 00:45:00 --qos interactive --account desi python scripts/mock_tools/pota2clus_fast.py --realization $i --mockver AbacusSummit_$3
 mv $SCRATCH/AbacusSummit_$3/mock$i/*GC* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_$3/mock$i/
 mv $SCRATCH/AbacusSummit_$3/mock$i/*nz* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_$3/mock$i/
 rm $SCRATCH/AbacusSummit_$3/mock$i/*
 chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_$3/mock$i/*clustering*
 chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_$3/mock$i/*nz*
done