#!/bin/bash
OUTBASE='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1_window/altmtl{MOCKNUM}'

cp /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl$1/mock$1/LSScats/*_HPmapcut*.ran.fits /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1_window/altmtl$1/mock$1/LSScats/

cp /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl$1/mock$1/LSScats/*frac_tlobs.fits /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1_window/altmtl$1/mock$1/LSScats/

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --tracer ELG_LOP --notqso y --minr 0 --maxr 18 --use_map_veto '_HPmapcut' --mkclusran y  --nz y --mkclusdat y --splitGC y --outmd 'notscratch' 

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --tracer LRG --notqso n --minr 0 --maxr 18  --use_map_veto '_HPmapcut' --mkclusran y  --nz y --mkclusdat y --splitGC y --outmd 'notscratch' 

python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output $OUTBASE --mockver ab_secondgen --mocknum $1  --tracer QSO --notqso n --minr 0 --maxr 18  --use_map_veto '_HPmapcut' --mkclusran y  --nz y --mkclusdat y --splitGC y --outmd 'notscratch' 

