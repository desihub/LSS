#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

scriptdir=/global/homes/d/desica/LSScode/LSS/scripts

for ((i=$1;i<=$2;i++ ))
do
  echo $i
  python $scriptdir/mock_tools/mkCat_SecondGen_amtl.py --mocknum $i --tracer BGS_BRIGHT-02 --absmagmd redshiftdep --simName SecondGenMocks/AbacusSummitBGS_v2 --outmd cfs --specdata kibo-v1 --mkclusran y --mkclusdat y --nz y --splitGC y --par y
done

