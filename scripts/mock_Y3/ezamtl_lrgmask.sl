#!/bin/bash
#SBATCH --account=desi
#SBATCH --ntasks=1
#SBATCH --time=40:00
#SBATCH --qos=shared
#SBATCH --cpus-per-task=16
#SBATCH --constraint=cpu
#SBATCH --array=15
##SBATCH --array=51-100
#SBATCH -J lrg_mask
#SBATCH -o ./stdout/%x_%a.o%j
#SBATCH -e ./stdout/%x_%a.e%j
#SBATCH --dependency=afterok:38715455_15

source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/installed_packages/desihub/LSS/py

LSSCODE="$HOME/installed_packages/desihub/LSS/"

#same as LSS/scripts/mock_tools/abamtl_lrgmask_sbatch_da2.sh

realization=$SLURM_ARRAY_TASK_ID
echo "realization=${realization}"

input="/pscratch/sd/j/jerryou/DA2/mocks/EZmock/forFA${realization}.fits"
output="/pscratch/sd/j/jerryou/DA2/mocks/EZmock/forFA${realization}_matched_input_full_lrg_imask.fits"

time python $LSSCODE/scripts/mock_tools/readwrite_pixel_bitmask_da2.py --tracer lrg --input $input --output $output --cat_type Generic --nproc $SLURM_CPUS_PER_TASK
