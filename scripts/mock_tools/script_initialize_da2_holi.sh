#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --qos=shared
##SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --constraint=cpu
#SBATCH --array=0-14
#SBATCH --account=desi
#SBATCH --mem=200G
#SBATCH --ntasks=1

source /global/common/software/desi/desi_environment.sh main
#source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
module load desitarget/3.0.0

#PYTHONPATH=$PYTHONPATH:$HOME/LSS/py

srun python initialize_amtl_mocks_da2.py /global/cfs/cdirs/desi/mocks/cai/holi/altMTL/forFA$SLURM_ARRAY_TASK_ID.fits /global/cfs/cdirs/desi/mocks/cai/holi/altMTL/altmtl$SLURM_ARRAY_TASK_ID DARK
