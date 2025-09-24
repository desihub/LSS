#!/bin/bash

#combine information from the assignment files (--combd y) and real data spec files (--joindspec y)
python $LSSCODE/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1 --mockver ab_secondgen --mocknum $1 --survey DA2 --add_gtl y --specdata kibo-v1 --tracer dark --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --combd y --joindspec y --par y --minr 0 --maxr 18 --outmd 'noscratch'

#get the dark_*full_noveto.ran.fits files
python $LSSCODE/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1 --mockver ab_secondgen --mocknum $1 --survey DA2 --add_gtl y --specdata kibo-v1 --tracer dark --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --fullr y --par y --minr 0 --maxr 18 --outmd 'noscratch'


