#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 00:15:00
#SBATCH --array=559-561
#test
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$SLURM_ARRAY_TASK_ID
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=GLAM-Uchuu_v2
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --add_gtl y --specdata loa-v1 --tracer dark --targDir $SCRATCH/DA2/mocks/$sim --combd y --joindspec y --par y --usepota y

