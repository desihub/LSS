#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
 j=$i+1
 srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/mkCat_SecondGen.py --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/altmtl{MOCKNUM} --mockver ab_secondgen --mockmin $i --mockmax $j --survey Y1 --add_gtl y --specdata iron --tracer $3 --notqso y --minr 0 --maxr 18 --fulld y --fullr y --apply_map_veto y --apply_veto y --use_map_veto _HPmapcut --mkclusran y --resamp y --nz y --mkclusdat y --targDir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit
done

