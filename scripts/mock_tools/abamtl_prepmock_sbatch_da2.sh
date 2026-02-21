#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --qos=shared
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --constraint=cpu
#SBATCH --array=0-24
#SBATCH --account=desi

source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

#Execute it as: >> sbatch abamtl_prepmock_sbatch_da2.sh DARK or sbatch abamtl_prepmock_sbatch_da2.sh BRIGHT
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system
srun $LSSCODE/scripts/mock_tools/run1_prepmock_LSS_DA2.sh $SLURM_ARRAY_TASK_ID $1
