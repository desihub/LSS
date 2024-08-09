#!/bin/bash
python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/ALTMTL_EZmock/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum $1 --survey Y1 --add_gtl y --specdata iron --tracer dark --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/ALTMTL_EZmock --combd y --joindspec y

#mv $SCRATCH/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$1/mock$1/* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$1/mock$1/
#mv $SCRATCH/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$1/fba$1/* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$1/fba$1/
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/ALTMTL_EZmock/altmtl$1/mock$1/*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/ALTMTL_EZmock/altmtl$1/fba$1/*
