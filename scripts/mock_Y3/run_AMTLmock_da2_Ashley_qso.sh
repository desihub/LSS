#!/bin/bash

#get QSO fulld (--fulld y), and masked data (--apply_veto y) and randoms (--apply_veto_ran y), and add tileloc info to randoms (--add_tlcomp y)
python $LSSCODE/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1 --mockver ab_secondgen --mocknum $1 --survey DA2 --add_gtl n --specdata kibo-v1 --tracer QSO --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --fulld y --apply_veto y --apply_veto_ran y --add_tlcomp y --par y --minr 0 --maxr 18 --outmd 'noscratch'


#get QSO clustering catalogs (--mkclusdat y)
python $LSSCODE/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1 --mockver ab_secondgen --mocknum $1 --survey DA2 --add_gtl n --specdata kibo-v1 --tracer QSO --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --mkclusdat y --par y --minr 0 --maxr 18 --outmd 'noscratch'


#get QSO clustering catalogs (--mkclusran y)
python $LSSCODE/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1 --mockver ab_secondgen --mocknum $1 --survey DA2 --add_gtl n --specdata kibo-v1 --tracer QSO --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --mkclusran y --par y --minr 0 --maxr 18 --outmd 'noscratch'

#split them NGC/SGC (--splitGC y)
python $LSSCODE/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1 --mockver ab_secondgen --mocknum $1 --survey DA2 --add_gtl n --specdata kibo-v1 --tracer QSO --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --splitGC y --par y --minr 0 --maxr 18 --outmd 'noscratch'

#refactor/add FKP weights (--nz y)
python $LSSCODE/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$1 --mockver ab_secondgen --mocknum $1 --survey DA2 --add_gtl n --specdata kibo-v1 --tracer QSO --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --nz y --par y --minr 0 --maxr 18 --outmd 'noscratch'
