#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 01:00:00
#SBATCH --array=6-9
#SBATCH --reservation=finish_dr2_mocks
##SBATCH --dependency=afterany:51952224

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
module load desitarget/3.0.0

srun --cpus-per-task=32 --ntasks=1 python $LSS/scripts/mock_tools/initialize_amtl_mocks_da2.py /pscratch/sd/d/desica/DA2/mocks/holi_bgs/forFA$SLURM_ARRAY_TASK_ID.fits /pscratch/sd/d/desica/DA2/mocks/holi_bgs/altmtl$SLURM_ARRAY_TASK_ID/ BRIGHT

