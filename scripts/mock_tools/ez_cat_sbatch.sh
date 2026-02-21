#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=875-876
#SBATCH --account=desi

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun scripts/mock_tools/run1_EZmock_LSS.sh $SLURM_ARRAY_TASK_ID