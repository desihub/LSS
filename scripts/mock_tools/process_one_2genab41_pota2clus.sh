#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py

ver=v4_1

srun -N 1 -C cpu -t 00:45:00 --qos interactive --account desi python scripts/mock_tools/pota2clus_fast.py --realization $1 --mockver AbacusSummit_$ver
mv $SCRATCH/AbacusSummit_$ver/mock$1/*GC* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_$ver/mock$1/
mv $SCRATCH/AbacusSummit_$ver/mock$1/*nz* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_$ver/mock$1/
rm $SCRATCH/AbacusSummit_$ver/mock$1/*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_$ver/mock$1/*clustering*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_$ver/mock$1/*nz*
