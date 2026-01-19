#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
#DR2-mocks-v1
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#export LSSCODE=$HOME/LSScode/LSS
#export PYTHONPATH=$PYTHONPATH:$LSSCODE/py
#
mocknum=102
scriptdir=$HOME/LSScode/LSS/scripts
#scriptdir=$LSSCODE/scripts
#/global/homes/d/desica/LSScode/LSS/scripts
sim=GLAM-Uchuu_v1
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --add_gtl y --specdata loa-v1 --tracer dark --targDir $SCRATCH/DA2/mocks/$sim --combd y --joindspec y --par y --usepota y
#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer LRG --fulld y --apply_veto y --par y
#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer QSO --fulld y --apply_veto y --par y
#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer ELG_LOP --notqso y --fulld y --apply_veto y --par y
#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer ELG --notqso y --fulld y --apply_veto y --par y


python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer QSO  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer LRG  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y


srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer QSO  --doimlin y --replace_syscol --par y --imsys_zbin split
srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer LRG  --doimlin y --replace_syscol --par y --imsys_zbin fine
srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --prep4sysnet y --nran4imsys 18
#$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELG_LOPnotqso $SCRATCH/DA2/mocks/$sim/altmtl$mocknum/loa-v1/mock$mocknum/LSScats/
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --addsysnet y --replace_syscol
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG --notqso y --par y --prep4sysnet y --nran4imsys 18
#$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELGnotqso $SCRATCH/DA2/mocks/$sim/altmtl$mocknum/loa-v1/mock$mocknum/LSScats/
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG --notqso y --par y --addsysnet y --replace_syscol
