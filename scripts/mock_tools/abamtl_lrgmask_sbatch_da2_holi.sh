#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --qos=regular
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
#SBATCH --array=500,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,648,649,650
#SBATCH --account=desi
#SBATCH --nodes=1

source /global/common/software/desi/desi_environment.sh main
#PYTHONPATH=$PYTHONPATH:$HOME/LSS

#Execute it as: >> sbatch abamtl_lrgmask_sbatch_da2.sh
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system
srun python readwrite_pixel_bitmask_da2.py --tracer lrg -i $SLURM_ARRAY_TASK_ID --cat_type Generic --nproc $SLURM_CPUS_PER_TASK --input /pscratch/sd/d/desica/DA2/mocks/holi_v1/forFA$SLURM_ARRAY_TASK_ID.fits --output /pscratch/sd/d/desica/DA2/mocks/holi_v1/forFA"$SLURM_ARRAY_TASK_ID"_matched_input_full_lrg_imask.fits
