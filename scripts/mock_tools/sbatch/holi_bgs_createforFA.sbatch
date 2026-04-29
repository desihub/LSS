#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 01:00:00
#SBATCH --array=150-199,201-999
#SBATCH --reservation=finish_dr2_mocks

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

scriptdir=/global/homes/d/desica/LSScode/LSS/

python $scriptdir/scripts/mock_tools/prepare_holi_bgs_seed.py $SLURM_ARRAY_TASK_ID

