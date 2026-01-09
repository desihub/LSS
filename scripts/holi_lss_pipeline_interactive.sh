#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
LSSCODE=$HOME/LSScode
PYTHONPATH=$PYTHONPATH:$LSSCODE/py
#module load LSS/DR2-mocks-v0
mocknum=201
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_v1
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --add_gtl y --specdata loa-v1 --tracer dark --targDir $SCRATCH/DA2/mocks/$sim --combd y --joindspec y --par y --usepota y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer LRG --fulld y --apply_veto y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer QSO --fulld y --apply_veto y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer ELG_LOP --notqso y --fulld y --apply_veto y --par y

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer QSO  --mkclusdat y #--mkclusran y --splitGC y --nz y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer LRG  --mkclusdat y #--mkclusran y --splitGC y --nz y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y  --mkclusdat y #--mkclusran y --splitGC y --nz y --par y


#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer QSO  --doimlin y --replace_syscol --par y --imsys_zbin split
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer LRG  --doimlin y --replace_syscol --par y --imsys_zbin fine
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --prep4sysnet y --nran4imsys 18
#$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELG_LOPnotqso $SCRATCH/DA2/mocks/$sim/altmtl$mocknum/loa-v1/mock$mocknum/LSScats/
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --addsysnet y --replace_syscol
