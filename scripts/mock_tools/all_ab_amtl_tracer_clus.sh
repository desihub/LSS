#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for ((i=$1;i<=$2;i++ ))
do
 srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum $i  --survey Y1 --add_gtl n --specdata iron --tracer $3 --notqso $4 --minr 0 --maxr 18 --fulld n --fullr n --apply_veto n --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --targDir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit
done

