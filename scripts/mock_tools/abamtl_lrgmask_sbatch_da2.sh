#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --qos=regular
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
#SBATCH --array=0-24
#SBATCH --account=desi
#SBATCH --nodes=1

source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

#Execute it as: >> sbatch abamtl_lrgmask_sbatch_da2.sh
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system
srun python $LSSCODE/scripts/mock_tools/readwrite_pixel_bitmask_da2.py --tracer lrg -i $SLURM_ARRAY_TASK_ID --cat_type Ab2ndgen --secgen_ver AbacusSummit_v4_1 --nproc $SLURM_CPUS_PER_TASK
