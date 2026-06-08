#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --qos=regular
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
#SBATCH --array=435-999
#SBATCH --account=desi
#SBATCH --nodes=1
#SBATCH -J lrgmask_holi
#SBATCH --reservation=altmtl

source /global/common/software/desi/desi_environment.sh main
#PYTHONPATH=$PYTHONPATH:$HOME/LSS/py

#Execute it as: >> sbatch abamtl_lrgmask_sbatch_da2.sh
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system
srun python readwrite_pixel_bitmask_da2.py --tracer lrg -i $SLURM_ARRAY_TASK_ID --cat_type Generic --nproc $SLURM_CPUS_PER_TASK --input /pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA$SLURM_ARRAY_TASK_ID.fits --output /pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA"$SLURM_ARRAY_TASK_ID"_matched_input_full_lrg_imask.fits
