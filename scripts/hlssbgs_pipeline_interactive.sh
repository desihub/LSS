#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$1
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_bgs
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --add_gtl y --specdata loa-v1 --tracer bright --targDir $SCRATCH/DA2/mocks/$sim --combd y --joindspec y --par y --usepota y

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer BGS_BRIGHT --fulld y --apply_veto y --par y 

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer BGS_BRIGHT-21.35  --mkclusdat y --mkclusran y --splitGC y --nz y --par y 

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer BGS_BRIGHT-21.35 --doimlin y --replace_syscol --par y --imsys_zbin split

python $scriptdir/mock_tools/mkCat_amtl.py --mocknum $mocknum --tracer LRG --simName $sim --transfer_cfs