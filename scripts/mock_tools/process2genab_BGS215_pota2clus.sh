#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
 srun -N 1 -C cpu -t 00:45:00 --qos interactive --account desi python scripts/mock_tools/pota2clus_fast.py --mockver AbacusSummitBGS --tracer BGS_BRIGHT-21.5 --realization $i --prog BRIGHT
 mv $SCRATCH/AbacusSummitBGS/mock$i/*GC* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/mock$i/
 mv $SCRATCH/AbacusSummitBGS/mock$i/*nz* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/mock$i/
 rm $SCRATCH/AbacusSummitBGS/mock$i/*
 chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/mock$i/*clustering*
 chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/mock$i/*nz*
done