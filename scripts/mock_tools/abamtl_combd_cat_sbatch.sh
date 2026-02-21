#!/bin/bash
#SBATCH --time=04:30:00
#SBATCH --qos=regular
##SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=1,2
#SBATCH --account=desi
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=460G

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:/pscratch/sd/a/acarnero/codes/LSS

srun /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/run1_AMTLmock_combd_LSS.sh $SLURM_ARRAY_TASK_ID
