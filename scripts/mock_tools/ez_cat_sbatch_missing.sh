#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=104,108,113,114,978,979,994
#SBATCH --account=desi

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun scripts/mock_tools/run1_EZmock_LSS_scratch.sh $SLURM_ARRAY_TASK_ID
srun scripts/mock_tools/run1_EZmock_LSS_mv.sh $SLURM_ARRAY_TASK_ID