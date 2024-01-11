#!/bin/bash

OUTBASE=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl{MOCKNUM}

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --survey Y1 --add_gtl y --specdata iron --tracer ELG_LOP --notqso y --minr 0 --maxr 18 --fulld y --fullr y --apply_veto y --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --targDir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --survey Y1 --add_gtl y --specdata iron --tracer LRG --notqso n --minr 0 --maxr 18 --fulld y --fullr y --apply_veto y --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --targDir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --survey Y1 --add_gtl y --specdata iron --tracer QSO --notqso n --minr 0 --maxr 18 --fulld y --fullr y --apply_veto y --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --targDir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1

mv $SCRATCH/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl$1/mock$1/LSScats/* /global/cfs/cdirs/desi/survey/catalogs//Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl$1/mock$1/LSScats/
chmod 775 /global/cfs/cdirs/desi/survey/catalogs//Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl$1/mock$1/LSScats/*