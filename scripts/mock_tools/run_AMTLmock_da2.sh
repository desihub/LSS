#!/bin/bash

#python mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum $1 --survey DA2 --specdata jura-v1 --add_gtl n --tracer QSO --notqso n --minr 0 --maxr 18 --outmd 'noscratch' --par y --fulld n --fullr n --apply_veto n --apply_veto_ran y --add_tl y --use_map_veto _HPmapcut --mkclusran y --mkclusdat y --splitGC y --nz y --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1

python mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum $1 --survey DA2 --specdata jura-v1 --add_gtl n --tracer ELG_LOP --notqso y --minr 0 --maxr 18 --outmd 'noscratch' --par y --fulld y --fullr y --apply_veto y --apply_veto_ran y --add_tl y --use_map_veto _HPmapcut --mkclusran y --mkclusdat y --splitGC y --nz y --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1



#python mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum $1 --survey Y1 --specdata iron --add_gtl n --tracer QSO --notqso n --minr 0 --maxr 18 --outmd 'scratch' --par y --fulld y --fullr y --apply_veto y --apply_veto_ran y --add_tl y --use_map_veto _HPmapcut --mkclusran y --mkclusdat y --splitGC y --nz y --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1

#python mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum $1 --survey Y1 --specdata iron --add_gtl n --tracer LRG --notqso n --minr 0 --maxr 18 --outmd 'scratch' --par y --fulld y --fullr y --apply_veto y --apply_veto_ran y --add_tl y --use_map_veto _HPmapcut --mkclusran y --mkclusdat y --splitGC y --nz y --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1


#python mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum 11  --survey Y1 --add_gtl n --specdata iron --tracer ELG_LOP --notqso y --fulld n --fullr n --apply_veto n --mkclusran y  --nz y --mkclusdat y --splitGC y --outmd 'notscratch' --use_map_veto _HPmapcut --compmd altmtl --add_bitweights /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl11_R128_v41/DARK_bitweights.fits --add_weight_ntile y --par y --maxr 18 --minr 0
