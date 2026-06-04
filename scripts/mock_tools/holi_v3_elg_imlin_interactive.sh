#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$1
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_v3
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

for ((i=$1;i<=$2;i++ )); do
  python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir /global/cfs/cdirs/desi/mocks/cai/LSS/ --outmd cfs --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y  --doimlin y --par y --imsys_zbin split
done
