#!/bin/bash
OUTBASE=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl{MOCKNUM}

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --add_gtl y --apply_veto y --tracer ELG_LOP --notqso y --minr 0 --maxr 18 --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --outmd 'notscratch'

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --add_gtl y --apply_veto y  --tracer LRG --notqso n --minr 0 --maxr 18  --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --outmd 'notscratch'

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --add_gtl y --apply_veto y --tracer QSO --notqso n --minr 0 --maxr 18  --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --outmd 'notscratch'

