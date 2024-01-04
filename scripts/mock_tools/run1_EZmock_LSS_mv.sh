#!/bin/bash

set -e

mv $SCRATCH/EZmock/FFA/mock$1/*GC* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/mock$1/
mv $SCRATCH/EZmock/FFA/mock$1/*nz* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/mock$1/
rm $SCRATCH/EZmock/FFA/mock$1/*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/mock$1/*clustering*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/mock$1/*nz*