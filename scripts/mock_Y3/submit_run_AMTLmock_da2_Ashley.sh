#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=11-24
#SBATCH --account=desi
#SBATCH --output=creation_LSS_mock.%j.%N.out

set -e
export LSSCODE=$HOME/project/DESI/LSS

source /global/common/software/desi/desi_environment.sh main
srun /global/homes/z/zxzhai/project/DESI/LSS_Y3_v3/run_AMTLmock_da2_Ashley_dark.sh $SLURM_ARRAY_TASK_ID

srun /global/homes/z/zxzhai/project/DESI/LSS_Y3_v3/run_AMTLmock_da2_Ashley_lrg.sh $SLURM_ARRAY_TASK_ID

srun /global/homes/z/zxzhai/project/DESI/LSS_Y3_v3/run_AMTLmock_da2_Ashley_elg.sh $SLURM_ARRAY_TASK_ID


source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
srun /global/homes/z/zxzhai/project/DESI/LSS_Y3_v3/run_AMTLmock_da2_Ashley_qso.sh $SLURM_ARRAY_TASK_ID
