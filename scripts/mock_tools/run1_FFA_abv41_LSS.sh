#!/bin/bash
#python scripts/readwrite_pixel_bitmask.py --tracer lrg --input $1 --cat_type Y1EZmock
python scripts/mock_tools/ffa2clus_fast.py --mockver AbacusSummit_v4_1 --realization $1
mv $SCRATCH/AbacusSummit_v4_1/mock$1/*GC* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/mock$1/
mv $SCRATCH/AbacusSummit_v4_1/mock$1/*nz* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/mock$1/
rm $SCRATCH/AbacusSummit_v4_1/mock$1/*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/mock$1/*clustering*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/mock$1/*nz*