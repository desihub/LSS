#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
#module load LSS/main
LSSCODE=$HOME/LSScode
PYTHONPATH=$LSSCODE/LSS/py:$PYTHONPATH
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

scriptdir=/global/homes/d/desica/LSScode/LSS/scripts

for ((i=$1;i<=$2;i++ ))
do
  echo $i
  python $scriptdir/mock_tools/mkCat_amtl.py --mocknum $i --tracer BGS_BRIGHT --absmagmd redshiftdep --simName SecondGenMocks/AbacusSummitBGS_v2  --specdata loa-v1 --par y --base_altmtl_dir '/global/cfs/cdirs/desi/survey/catalogs/' --combd y --usepota y --joindspec y --fulld y --apply_veto y
  python $scriptdir/mock_tools/mkCat_amtl.py --mocknum $i --tracer BGS_BRIGHT-02 --absmagmd redshiftdep --simName SecondGenMocks/AbacusSummitBGS_v2  --specdata loa-v1 --par y --base_altmtl_dir '/global/cfs/cdirs/desi/survey/catalogs/'  --mkclusran y --mkclusdat y --nz y --splitGC y --par y --apply_oldfoot y --doimlin y --replace_syscol --transfer_cfs
done

