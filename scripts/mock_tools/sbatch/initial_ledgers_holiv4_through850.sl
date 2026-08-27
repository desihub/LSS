#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 01:30:00
#SBATCH --array=601-850
#SBATCH --reservation=_CAP_holi_dr3_altmtl
#SBATCH --dependency=afterany:55383732
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
module load desitarget/3.0.0

mocknum=$SLURM_ARRAY_TASK_ID
scriptdir=/global/u2/d/desica/LSScode/LSS/scripts/mock_tools/

srun python $scriptdir/initialize_amtl_mocks_da3.py /pscratch/sd/d/desica/DA3/mocks/holi_v4/altmtl/forFA"$mocknum".fits /pscratch/sd/d/desica/DA3/mocks/holi_v4/altmtl/altmtl"$mocknum"/ DARK
