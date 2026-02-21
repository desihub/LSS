#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --qos=debug
#SBATCH --constraint=cpu
#SBATCH --array=0-4
#SBATCH --account=desi
#SBATCH --mem=5G


source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun /pscratch/sd/a/acarnero/codes/LSS/bin/run_altMTL_batch.sh $SLURM_ARRAY_TASK_ID
