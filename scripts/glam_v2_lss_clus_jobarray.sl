#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 00:25:00
#SBATCH --array=200-249,550-599,850-899,1100-1149
#SBATCH --dependency=afterany:51314869
#test
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$SLURM_ARRAY_TASK_ID
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=GLAM-Uchuu_v2
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer QSO  --mkclusdat y --mkclusran y --splitGC y --nz y --par y 
srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer LRG  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y

