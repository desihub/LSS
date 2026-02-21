#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=1-24
#SBATCH --account=desi

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun scripts/mock_tools/process_one_DR2_2genab41_pota2clus.sh $SLURM_ARRAY_TASK_ID