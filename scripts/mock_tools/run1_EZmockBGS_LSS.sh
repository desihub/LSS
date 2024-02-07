#!/bin/bash
python scripts/mock_tools/ffa2clus_fast.py --mockver EZmock/FFA_BGS --tracer BGS --realization $1
mv $SCRATCH/EZmock/FFA_BGS/mock$1/*GC* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/mock$1/
mv $SCRATCH/EZmock/FFA_BGS/mock$1/*nz* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/mock$1/
rm $SCRATCH/EZmock/FFA_BGS/mock$1/*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/mock$1/*clustering*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/mock$1/*nz*