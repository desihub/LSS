#!/bin/bash
python scripts/readwrite_pixel_bitmask.py --tracer lrg --input $1 --cat_type Y1EZmock
python scripts/mock_tools/ffa2clus_fast.py --mockver EZmock/FFA --realization $1
cp $SCRATCH/EZmock/FFA/mock$1/* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/mock$1/
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/mock$1/*clustering*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/mock$1/*nz*