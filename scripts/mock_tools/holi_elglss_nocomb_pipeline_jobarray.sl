#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 02:30:00
#SBATCH --array=211-350
#test
source /global/common/software/desi/desi_environment.sh main
module load LSS/DR2-mocks-v0
mocknum=$SLURM_ARRAY_TASK_ID
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_v1
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer ELG --notqso y --fulld y --apply_veto y --par y


source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y


srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG --notqso y --par y --prep4sysnet y --nran4imsys 18
$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELGnotqso $SCRATCH/DA2/mocks/$sim/altmtl$mocknum/loa-v1/mock$mocknum/LSScats/
srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG --notqso y --par y --addsysnet y --replace_syscol

srun python scripts/mock_tools/mkCat_amtl.py --mocknum $mocknum --tracer LRG --simName $sim --transfer_cfs