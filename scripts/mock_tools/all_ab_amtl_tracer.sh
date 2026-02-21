#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for ((i=$1;i<=$2;i++ ))
do
 srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /pscratch/sd/a/acarnero/SecondGenMocks/AbacusSummit/altmtl{MOCKNUM} --mockver ab_secondgen --mocknum $i  --survey Y1 --add_gtl y --specdata iron --tracer $3 --notqso $4 --minr 0 --maxr 18 --fulld y --fullr y --apply_veto y --use_map_veto _HPmapcut --mkclusran y  --nz y --mkclusdat y --splitGC y --targDir /pscratch/sd/a/acarnero/SecondGenMocks/AbacusSummit
done

