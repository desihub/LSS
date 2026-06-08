#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_v3

for ((i=$1;i<=$2;i++ )); do
 python $scriptdir/mock_tools/mkCat_amtl.py --mocknum $i --tracer LRG --simName $sim --transfer_cfs
done