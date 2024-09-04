#!/bin/bash
python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum $1 --survey DA2 --add_gtl y --specdata jura-v1 --tracer dark --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1 --combd y --joindspec y

#mv $SCRATCH/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$1/mock$1/* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$1/mock$1/
#mv $SCRATCH/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$1/fba$1/* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$1/fba$1/
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1/mock$1/*
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1/fba$1/*
