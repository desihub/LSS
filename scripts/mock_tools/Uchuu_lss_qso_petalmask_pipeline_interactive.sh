#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
#module load LSS/main
#DR2-mocks-v1
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#export LSSCODE=$LSS
#export PYTHONPATH=$PYTHONPATH:$LSSCODE/py
#
mocknum=$1
scriptdir=$HOME/LSScode/LSS/scripts
#scriptdir=$LSSCODE/scripts
#/global/homes/d/desica/LSScode/LSS/scripts
sim=GLAM-Uchuu_v2/
bdir=/global/cfs/cdirs/desi/mocks/cai/LSS/
edir=petalmask
export LSSCODE=/global/homes/d/desica/LSScode/LSS
PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH


srun --exclusive -N 1 -n 1 --cpus-per-task=128 python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $bdir --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer QSO  --mkclusdat y --mkclusran y --splitGC y --nz y --par y   --extra_clusdir $edir --outmd cfs  --mask_petalnight_tids --doimlin y --replace_syscol --par y --imsys_zbin split --transfer_cfs


