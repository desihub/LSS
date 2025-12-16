#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH --array=0-48
#SBATCH -t 01:00:00

source /global/common/software/desi/desi_environment.sh main
module load LSS/DR2-mocks-v0
mockver=holi_v1
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts


source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test

srun python $scriptdir/mock_tools/pota2clus_fast.py --realization $SLURM_ARRAY_TASK_ID --base_dir /global/cfs/cdirs/desi/mocks/cai/LSS/ --data_dir /global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/ --mockver $mockver   --specrel loa-v1