#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH --array=559-648
#SBATCH --dependency=afterany:51250457
#test
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$SLURM_ARRAY_TASK_ID
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=GLAM-Uchuu_v2
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer QSO  --doimlin y --replace_syscol --par y --imsys_zbin split &
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer LRG  --doimlin y --replace_syscol --par y --imsys_zbin fine &
wait

