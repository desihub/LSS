#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=1-24
#SBATCH --account=desi

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun scripts/mock_tools/run1_AMTLmock_LSS_3_1fix.sh $SLURM_ARRAY_TASK_ID