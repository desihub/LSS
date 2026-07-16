#!/usr/bin/env bash

# Description:
# ===========
# top level of Holi pipeline with slurm
# job array to process chunk of siluation = cpus-per-task

#SBATCH --ntasks=1 # force srun with option -n 1 , ie 1 task
#SBATCH --cpus-per-task=32
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -J HoliArray
#SBATCH -t 06:00:00
#SBATCH --array=0-2
#SBATCH --output=holi_%j.out
#SBATCH --error=holi_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jcolley@lpnhe.in2p3.fr

#
# Environment
#
#source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

LSS_DIR=/global/cfs/cdirs/desi/users/colley/LSS
DS_DIR=/global/cfs/cdirs/desi/mocks/fa4acm/holi/webjax_v4.80

# 
HOLI_DIR=$LSS_DIR/scripts/mock_tools/pipe_holi
cd $HOLI_DIR

# run chunk
# chunk size is the number of CPU : cpus-per-task
./slurm2_chunk_all_step.sh $LSS_DIR $DS_DIR
