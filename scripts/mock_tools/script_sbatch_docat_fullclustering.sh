#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --qos=regular
##SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=1-2
#SBATCH --account=desi
#SBATCH --mem=384G
#SBATCH --cpus-per-task=9
#SBATCH --ntasks=1 
##SBATCH --exclusive

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#PYTHONPATH=$PYTHONPATH  #:$HOME/LSS

srun /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/run_AMTLmock_v4_1.sh $SLURM_ARRAY_TASK_ID
