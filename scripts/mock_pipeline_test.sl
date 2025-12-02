#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q debug
#SBATCH -t 00:30:00

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
mocknum=201
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --add_gtl y --specdata loa-v1 --tracer dark --targDir $SCRATCH/DA2/mocks/Holi/v3.00 --combd y --joindspec y --par y --usepota y
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/Holi/v3.00 --tracer LRG --fulld y --apply_veto y --par y
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/Holi/v3.00 --tracer QSO --fulld y --apply_veto y --par y
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/Holi/v3.00 --tracer ELG_LOP --notqso y --fulld y --apply_veto y --par y

srun python $LSSCODE/LSS/scripts/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer LRG  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
srun python $LSSCODE/LSS/scripts/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test

srun python $LSSCODE/LSS/scripts/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer QSO  --mkclusdat y --mkclusran y --splitGC y --nz y --par y



#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer LRG  --doimlin y --replace_syscol --par y --imsys_zbin fine
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --prep4sysnet y --nran4imsys 18
#$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELG_LOPnotqso $SCRATCH/DA2/mocks/Holi/v3.00/altmtl$mocknum/loa-v1/mock$mocknum/LSScats/
#srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName Holi/v3.00 --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --addsysnet y --replace_syscol
