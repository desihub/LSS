#!/bin/bash
#SBATCH --time=00:50:00
#SBATCH --qos=shared
#SBATCH --cpus-per-task=128
#SBATCH --constraint=cpu
#SBATCH --array=1
#SBATCH --account=desi
#SBATCH --mem=40G

source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py

#Execute it as: >> sbatch abamtl_applyphotmask_sbatch_da2.sh DARK or abamtl_applyphotmask_sbatch_da2.sh BRIGHT
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system
srun python $LSSCODE/scripts/mock_tools/apply_phot_mask.py $SLURM_CPUS_PER_TASK $SLURM_ARRAY_TASK_ID $1
