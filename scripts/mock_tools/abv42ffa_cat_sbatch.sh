#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=2-24
#SBATCH --account=desi

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun scripts/mock_tools/run1_FFA_abv42_LSS.sh $SLURM_ARRAY_TASK_ID