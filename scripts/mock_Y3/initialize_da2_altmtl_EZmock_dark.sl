#!/bin/bash
#SBATCH --time=30:00
#SBATCH --qos=shared
#SBATCH --cpus-per-task=80
#SBATCH --constraint=cpu
##SBATCH --array=2-50
#SBATCH --array=80-80
#SBATCH --account=desi
#SBATCH --mem=160G
#SBATCH --ntasks=1
#SBATCH -J prepare_mtl
#SBATCH -o ./stdout/%x_%a.o%j
#SBATCH -e ./stdout/%x_%a.e%j


source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:/global/u2/j/jerryou/installed_packages/desihub/LSS

seed=$SLURM_ARRAY_TASK_ID

#printf -v realization "%04d" $seed

codepath="/global/u2/j/jerryou/installed_packages/desihub/LSS/scripts/mock_tools/"

srun -n 1 --cpu-bind=none python ${codepath}/initialize_amtl_mocks_da2.py /pscratch/sd/j/jerryou/DESI_Y3/SecondGenMocks/EZmock/forFA/forFA${seed}.fits /pscratch/sd/j/jerryou/DESI_Y3/SecondGenMocks/EZmock/altmtl${seed} DARK

