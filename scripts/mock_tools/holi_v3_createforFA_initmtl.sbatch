#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 01:00:00
#SBATCH --array=450-999

bash scripts/mock_tools/run1_456.sh $SLURM_ARRAY_TASK_ID webjax_v4.81 webjax_v4.80 webjax_v4.80 /pscratch/sd/d/desica/DA2/mocks/holi_v3 holi
source /global/common/software/desi/desi_environment.sh main
module load desitarget/3.0.0
module load LSS/main

srun --cpus-per-task=64 --ntasks=1 python $LSS/scripts/mock_tools/initialize_amtl_mocks_da2.py /pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA$SLURM_ARRAY_TASK_ID.fits /pscratch/sd/d/desica/DA2/mocks/holi_v3/altmtl$SLURM_ARRAY_TASK_ID/ DARK
