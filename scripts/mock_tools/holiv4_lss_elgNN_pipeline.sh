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
sim=holi_v4/altmtl
bdir=/global/cfs/cdirs/desi/mocks/cai/LSS/
edir=NN
export LSSCODE=/global/homes/d/desica/LSScode/LSS
PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH


srun --exclusive -N 1 -n 1 --cpus-per-task=128 python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $bdir --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y --redo_fracz y --compmd data --nearestneighbor y  --extra_clusdir $edir --outmd cfs --prep4sysnet y --nran4imsys 18


#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $bdir --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --prep4sysnet y --nran4imsys 18 --extra_clusdir $edir --outmd cfs

$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELG_LOPnotqso $bdir/DA2/mocks/$sim/altmtl$mocknum/loa-v1/mock$mocknum/LSScats/$edir/

srun --exclusive -N 1 -n 1 --cpus-per-task=128  python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $bdir --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --addsysnet y --replace_syscol --extra_clusdir $edir --outmd cfs

