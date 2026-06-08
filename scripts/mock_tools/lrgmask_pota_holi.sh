#!/bin/bash
#SBATCH --time=02:20:00
#SBATCH --qos=regular
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
##SBATCH --array=451-470
##SBATCH --array=0,1,2,5,7,8,9,10,11,12,13,14,23,24,25,28,29,30,31,33,34,35,38,39,43,44,45,48,49,52,53,54,56,57,58,59,60,61,62,64,65,67,68,69,91,95,100,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223
#SBATCH --array=224,229,230,242,243,253,255,256,257,258,264,266,268,272,276,277,282,283,284,287,289,291,294,296,297,308,312,313
#SBATCH --account=desi
#SBATCH --nodes=1

source /global/common/software/desi/desi_environment.sh main
#PYTHONPATH=$PYTHONPATH:$HOME/LSS/py


#Execute it as: >> sbatch abamtl_lrgmask_sbatch_da2.sh
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system

srun python readwrite_pixel_bitmask_da2.py --tracer lrg -i $SLURM_ARRAY_TASK_ID --cat_type Generic --nproc $SLURM_CPUS_PER_TASK --input /pscratch/sd/d/desica/DA2/mocks/holi_v2/forFA$SLURM_ARRAY_TASK_ID.fits --output /pscratch/sd/d/desica/DA2/mocks/holi_v2/forFA"$SLURM_ARRAY_TASK_ID"_matched_input_full_lrg_imask.fits

srun python getpotaDA2_mock.py --realization $SLURM_ARRAY_TASK_ID --mock holi --base_output /pscratch/sd/d/desica/DA2/mocks/holi_v2/ --base_input /pscratch/sd/d/desica/DA2/mocks/holi_v2/forFA$SLURM_ARRAY_TASK_ID.fits


#printf -v padded_i "%04d" $SLURM_ARRAY_TASK_ID

#srun python readwrite_pixel_bitmask_da2.py --tracer lrg -i $SLURM_ARRAY_TASK_ID --cat_type Generic --nproc $SLURM_CPUS_PER_TASK --input /pscratch/sd/d/desica/DA2/mocks/holi_v2/forFA$padded_i.fits --output /pscratch/sd/d/desica/DA2/mocks/holi_v2/forFA"$SLURM_ARRAY_TASK_ID"_matched_input_full_lrg_imask.fits

#srun python getpotaDA2_mock.py --realization $SLURM_ARRAY_TASK_ID --mock holi --base_output /pscratch/sd/d/desica/DA2/mocks/holi_v2/ --base_input /pscratch/sd/d/desica/DA2/mocks/holi_v2/forFA$padded_i.fits
