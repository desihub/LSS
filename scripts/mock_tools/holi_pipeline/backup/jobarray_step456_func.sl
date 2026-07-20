#!/bin/bash
#SBATCH --ntasks=1 # force srun with option -n 1 , ie 1 task
#SBATCH --cpus-per-task=3
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -J FA456
#SBATCH -t 0:20:00
#SBATCH --array=12-80
#SBATCH --output=holi_%j.out
#SBATCH --error=holi_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jcolley@lpnhe.in2p3.fr

#
# Environment
#
#source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

# module load LSS/main/global/cfs/cdirs/desi/users/colley/fa_snake/scripts/concatenate_tracers_to_fba_snake.py
LSS_DIR=/global/cfs/cdirs/desi/users/colley/LSS
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH

PROC_DIR=$LSS_DIR/scripts/mock_tools/pipe_holi
DS_DIR=/global/cfs/cdirs/desi/mocks/fa4acm/holi/webjax_v4.80


#
# process one seed
#
$PROC_DIR/step_456_seed.sh  $LSS_DIR  $SLURM_ARRAY_TASK_ID  $DS_DIR 